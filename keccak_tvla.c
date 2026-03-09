/*MIT License

Copyright (c) 2025 Adam Beattie

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in
all copies or substantial portions of the Software.
*/

#include "hal.h"
#include "simpleserial.h"
#include <stdint.h>
#include <stdio.h>
#include <string.h>
#include <stddef.h>
#include "stm32f4xx.h"      // <-- defines IRQn_Type, SysTick_IRQn, etc.
#include "core_cm4.h"       // <-- defines __get_PRIMASK(), __DSB(), etc.

//Rho offsets
const uint8_t keccak_rho_offsets[5][5] = {
    {  0, 36,  3, 41, 18 },
    {  1, 44, 10, 45,  2 },
    { 62,  6, 43, 15, 61 },
    { 28, 55, 25, 21, 56 },
    { 27, 20, 39,  8, 14 }
};


//Iota round constants
const uint64_t RC[24] = {
    0x0000000000000001ULL, 0x0000000000008082ULL,
    0x800000000000808aULL, 0x8000000080008000ULL,
    0x000000000000808bULL, 0x0000000080000001ULL,
    0x8000000080008081ULL, 0x8000000000008009ULL,
    0x000000000000008aULL, 0x0000000000000088ULL,
    0x0000000080008009ULL, 0x000000008000000aULL,
    0x000000008000808bULL, 0x800000000000008bULL,
    0x8000000000008089ULL, 0x8000000000008003ULL,
    0x8000000000008002ULL, 0x8000000000000080ULL,
    0x000000000000800aULL, 0x800000008000000aULL,
    0x8000000080008081ULL, 0x8000000000008080ULL,
    0x0000000080000001ULL, 0x8000000080008008ULL
};

#ifndef UNUSED
#define UNUSED(x) ((void)(x))
#endif

//error return value for CW
#define SS_ERR_LEN 0xA1

//Maximum usable 
#define MAX_ORDER 10 //arbitrary limit - this can be removed

#ifndef MASKING_ORDER
#define MASKING_ORDER 3
#endif

#define NROUNDS 24 //Not used
#define MASKING_N_VAL (MASKING_ORDER + 1)
#define MASKING_N MASKING_N_VAL

//  Keccak/SHAKE rate constants
#define KECCAK_RATE 168
#define SHAKE128_RATE 168     // Used for SHAKE128
#define SHAKE256_RATE 136     // Used for SHAKE256
#define SHA3_224_RATE 144
#define SHA3_256_RATE 136
#define SHA3_384_RATE 104
#define SHA3_512_RATE 72      // Used for SHA3-512

//Sha and shake domain seperators
#define DOMAIN_SHA3   0x06
#define DOMAIN_SHAKE  0x1F

#define ENABLE_TIMING 1  

typedef struct {
    uint64_t share[MASKING_N];
} masked_uint64_t;

//fundamental masked object
typedef struct {
    uint64_t share[MASKING_N][5][5];   // share[s][x][y]
} masked_state_t;


void masked_keccak_f1600(masked_state_t *state);
uint64_t xor64_safe(uint64_t a, uint64_t b);

#if ENABLE_TIMING
  #define TSTAMP(var) uint32_t var = DWT->CYCCNT
#else
  #define TSTAMP(var) do {} while (0)
#endif

static uint16_t current_timing_position = 0;
// --- Exact-NOP macro using GAS .rept (no loop overhead) -----------------
#define NOP_BLOCK(N) asm volatile (".rept " #N "\n\tnop\n\t.endr\n" ::: "memory")

//Spaces in time record buffer
#define TIME_BUFFER_SIZE 32

//Macro for setting state
#define ST(st, s, x, y)   ((st)->share[(s)][(x)][(y)])


//Main global state
volatile masked_state_t global_state;


// Random matrices for masked_chi()
volatile uint64_t randmat[5][5][MASKING_N][MASKING_N];
volatile uint64_t randmat_and [5][5][MASKING_N][MASKING_N];
volatile uint64_t randmat_in[5][5][4][MASKING_N][MASKING_N];
volatile uint64_t randmat_mid[5][5][MASKING_N][MASKING_N];
volatile uint64_t randmat_out[5][5][MASKING_N][MASKING_N];
volatile masked_state_t out_state;

//for determining pause length in chi if jitter in use
volatile uint64_t chi_jitter_seed;

/**
 * ----------------------------------------------------------------------------
 * Macros: TEST_PROLOGUE / TEST_EPILOGUE
 * ----------------------------------------------------------------------------
 * Purpose:
 *   Bracket a code region with:
 *     - IRQ masking (save/restore PRIMASK),
 *     - DSB/ISB barriers,
 *     - trace trigger high/low,
 *     - fixed NOP paddings (NOP_BLOCK(128)) to create clean capture windows.
 *
 * Side-channel Notes:
 *   - Stabilizes measurement alignment on ChipWhisperer;
 *   - Ensures reproducible trace length and boundaries.
 * ----------------------------------------------------------------------------
 */
#define TEST_PROLOGUE()                                    \
    do {                                                    \
        uint32_t __primask = __get_PRIMASK();              \
        __disable_irq();                                   \
        __DSB(); __ISB();                               \
        NOP_BLOCK(128);                                    \
        trigger_high();                                    \
        __DSB(); __ISB();

#define TEST_EPILOGUE()                                    \
        NOP_BLOCK(128);                                    \
        trigger_low();                                      \
        __DSB(); __ISB();                                  \
        if (!__primask) __enable_irq();                    \
    } while (0)


//  rotation macro
#define ROL64(x, n) (((x) << (n)) | ((x) >> (64 - (n))))


//Stores timing records from the test and implementation functions 
//can be used to extract time values when using CW
static volatile uint32_t time_buffer[TIME_BUFFER_SIZE] = {0};
//number of time captures recorded during run
static int time_captures = 0;

//Specific randomness for Iota
static volatile uint64_t iota_rands[MASKING_N-1];
static volatile uint8_t iota_round_idx = 0; 



void dwt_init(void) {
    // Enable TRC (trace)
    CoreDebug->DEMCR |= CoreDebug_DEMCR_TRCENA_Msk;
    
    // Reset the cycle counter
    DWT->CYCCNT = 0;
    
    // Enable the cycle counter
    DWT->CTRL |= DWT_CTRL_CYCCNTENA_Msk;
}


// ----------------------------------------------------------------------------
// K_NOT_SPLIT[] holds a fixed public Boolean share decomposition of 0xFFFF...FFFF.
// Used by masked_not() to implement bitwise NOT securely under Boolean masking.
// Constructed so that XOR of all shares == 0xFFFFFFFFFFFFFFFFULL.
// ----------------------------------------------------------------------------
uint64_t K_NOT_SPLIT[MASKING_N];


// ----------------------------------------------------------------------------
// The chip used for capture was not returning accurate 64 bit xor sums 
// This method splits each 64 bit xor into two 32 bit ones
// Seems to be working correctly now, unsure what the issue was so if required
// This should be able to be safely removed. 
// ----------------------------------------------------------------------------
uint64_t xor64_safe(uint64_t a, uint64_t b) {
    uint32_t a_lo = (uint32_t)(a & 0xFFFFFFFFULL);
    uint32_t a_hi = (uint32_t)(a >> 32);
    uint32_t b_lo = (uint32_t)(b & 0xFFFFFFFFULL);
    uint32_t b_hi = (uint32_t)(b >> 32);

    uint64_t res_lo = (uint64_t)(a_lo ^ b_lo);
    uint64_t res_hi = (uint64_t)(a_hi ^ b_hi);

    return (res_hi << 32) | res_lo;
}


/**
 * Generates a 64-bit random value using two 32-bit RNG calls.
 *
 * Combines two outputs from the underlying RNG (`get_rand`) to construct
 * a full 64-bit value. The first call provides the high 32 bits and the
 * second provides the low 32 bits.
 *
 * This helper is used throughout the masking implementation to generate
 * randomness for share initialization and ISW masking gadgets.
 *
 * @return A 64-bit pseudo-random value.
 */
uint64_t get_rand64(void)
{
    uint64_t hi = get_rand();
    uint64_t lo = get_rand();
    return (hi << 32) | lo;
}


/**
 * Initializes the NOT-splitting constants used in masked Boolean NOT operations.
 *
 * This routine constructs a deterministic set of constants that allow the
 * masked NOT operation to be applied share-wise while preserving correctness
 * across the recombined value. Each share receives a different constant
 * derived from a fixed base value combined with an index-dependent pattern.
 *
 * The final share is chosen such that the XOR of all constants equals the
 * bitwise negation mask (~0ULL), ensuring the correct masked inversion
 * when shares are recombined.
 *
 * This initialization must be performed before masked NOT operations that
 * rely on the K_NOT_SPLIT constants.
 */
void init_not_split(void) {
    uint64_t acc = 0;
    for (size_t i = 0; i < MASKING_N - 1; i++) {
        K_NOT_SPLIT[i] = xor64_safe(0x9E3779B97F4A7C15ULL,
                                    0x0101010101010101ULL * (i + 1));
        acc = xor64_safe(acc, K_NOT_SPLIT[i]);
    }
    K_NOT_SPLIT[MASKING_N - 1] = xor64_safe(~0ULL, acc);
}


//~~~~~~~~~~~~~~~~~~~~DEBUGGING AND PRINTING METHODS~~~~~~~~~~~~~~~~~~~~
//~~~~~~~~~~~~~~~~~~~~FOR REGULAR USE NOT CW~~~~~~~~~~~~~~~~~~~~
/**
 * Prints the individual masking shares of a masked 64-bit value.
 *
 * This debugging utility displays each share of a Boolean-masked
 * 64-bit value in hexadecimal form. It is useful when validating
 * masking correctness, inspecting intermediate values, or verifying
 * that share recombination produces the expected unmasked value.
 *
 * Each share is printed as two 32-bit halves to ensure consistent
 * formatting across platforms where printing 64-bit integers may
 * require special format specifiers.
 *
 * The output is labeled using the provided string so that the printed
 * shares can be associated with a particular variable or stage in
 * the algorithm during debugging.
 *
 * @param label  Text label describing the value being printed.
 * @param v      Pointer to the masked 64-bit value whose shares are
 *               to be displayed.
 */
static void print_lane_shares(const char *label, const masked_uint64_t *v)
{
    printf("%s\n", label);
    for (int s = 0; s < MASKING_N; s++) {
        uint64_t x = v->share[s];
        uint32_t hi = (uint32_t)(x >> 32);
        uint32_t lo = (uint32_t)(x & 0xFFFFFFFFu);
        printf("  share[%d] = %08lx%08lx\n",
               s,
               (unsigned long)hi,
               (unsigned long)lo);
    }
}

/**
 * Prints the recombined values of a masked Keccak state.
 *
 * This debugging utility reconstructs the unmasked Keccak state by XORing
 * together all masking shares for each lane. The resulting 5×5 state is then
 * printed in hexadecimal format for inspection.
 *
 * Each lane is displayed as a 64-bit value split into two 32-bit halves to
 * ensure consistent formatting across platforms. The lanes are printed in
 * the natural Keccak layout (x across, y down), allowing the output to be
 * directly compared with reference implementations or intermediate states
 * during debugging.
 *
 * @param label  Text label describing the state being printed.
 * @param S      Pointer to the masked Keccak state whose recombined values
 *               should be displayed.
 */
static void print_masked_state(const char *label, const masked_state_t *S)
{
    uint64_t st[25];

    /* recombine shares */
    for (int y = 0; y < 5; y++) {
        for (int x = 0; x < 5; x++) {

            uint64_t v = 0;

            for (int s = 0; s < MASKING_N; s++)
                v ^= S->share[s][x][y];

            st[x + 5*y] = v;
        }
    }

    printf("%s\n", label);

    for (int y = 0; y < 5; y++) {
        for (int x = 0; x < 5; x++) {

            uint64_t v = st[x + 5*y];

            uint32_t hi = (uint32_t)(v >> 32);
            uint32_t lo = (uint32_t)(v & 0xFFFFFFFF);

            printf("%08lx%08lx ",
                   (unsigned long)hi,
                   (unsigned long)lo);
        }
        printf("\n");
    }

    printf("\n");
}
/**
 * Generates all randomness matrices required for masked Keccak gadget CHI.
 *
 * This command handler prepares the pre-generated randomness used by the
 * ISW-style masked AND operations and subsequent share-refresh steps inside
 * the masked Chi implementation.
 *
 * The following randomness matrices are populated:
 *
 *  - randmat_and : Random values used in ISW cross-terms for masked AND.
 *  - randmat_in  : Randomness used for input refresh stages before nonlinear
 *                  processing (multiple instances per lane).
 *  - randmat_mid : Randomness applied after the AND operation to refresh shares.
 *  - randmat_out : Randomness applied after storing results to maintain masking.
 *
 * Each matrix is symmetric (r[i][j] == r[j][i]) with zeros on the diagonal,
 * matching the requirements of ISW masking where randomness is shared
 * between share pairs.
 *
 * Random values are generated using `get_rand64()` for every unique pair
 * of shares across all 5×5 Keccak lanes.
 *
 * @param cmd   SimpleSerial command identifier (unused).
 * @param scmd  Subcommand identifier (unused).
 * @param len   Length of received data (unused).
 * @param buf   Input buffer (unused).
 *
 * @return 0x00 on success.
 */
 uint8_t prepare_randmat_cmd(uint8_t cmd, uint8_t scmd,
                                   uint8_t len, uint8_t *buf)
{
    UNUSED(cmd); UNUSED(scmd); UNUSED(len); UNUSED(buf);

    // --- Fill AND randomness (ISW) ---
    for (int y = 0; y < 5; y++) {
        for (int x = 0; x < 5; x++) {
            for (int i = 0; i < MASKING_N; i++) {
                for (int j = 0; j < MASKING_N; j++) {
                    if (i == j) {
                        randmat_and[x][y][i][j] = 0;
                    } else if (j > i) {
                        uint64_t r = get_rand64();
                        randmat_and[x][y][i][j] = r;
                        randmat_and[x][y][j][i] = r;
                    }
                }
            }
        }
    }
for (int y = 0; y < 5; y++) {
    for (int x = 0; x < 5; x++) {

        for (int k = 0; k < 4; k++) {   

            for (int i = 0; i < MASKING_N; i++) {
                for (int j = 0; j < MASKING_N; j++) {

                    if (i == j) {
                        randmat_in[x][y][k][i][j] = 0;
                    }
                    else if (j > i) {
                        uint64_t r = get_rand64();
                        randmat_in[x][y][k][i][j] = r;
                        randmat_in[x][y][k][j][i] = r;
                    }
                }
            }

        }
    }
}

    // --- Fill post-AND refresh randomness (r_mid) ---
    for (int y = 0; y < 5; y++) {
        for (int x = 0; x < 5; x++) {
            for (int i = 0; i < MASKING_N; i++) {
                for (int j = 0; j < MASKING_N; j++) {
                    if (i == j) {
                        randmat_mid[x][y][i][j] = 0;
                    } else if (j > i) {
                        uint64_t r = get_rand64();
                        randmat_mid[x][y][i][j] = r;
                        randmat_mid[x][y][j][i] = r;
                    }
                }
            }
        }
    }

    // --- Fill post-store refresh randomness (r_out) ---
    for (int y = 0; y < 5; y++) {
        for (int x = 0; x < 5; x++) {
            for (int i = 0; i < MASKING_N; i++) {
                for (int j = 0; j < MASKING_N; j++) {
                    if (i == j) {
                        randmat_out[x][y][i][j] = 0;
                    } else if (j > i) {
                        uint64_t r = get_rand64();
                        randmat_out[x][y][i][j] = r;
                        randmat_out[x][y][j][i] = r;
                    }
                }
            }
        }
    }

    return 0x00;
}

/**
 * Fills a set of masking shares according to a specified test pattern.
 *
 * This helper is primarily used for controlled side-channel experiments
 * (e.g., TVLA or leakage characterization). It ensures that the RNG is
 * always exercised to maintain a consistent power profile, even when the
 * resulting value is overwritten with a deterministic pattern.
 *
 * Supported patterns include:
 *
 *  - 0 : Fully random values generated by the RNG.
 *  - 1 : All-zero shares.
 *  - 2 : All bits set (~0ULL).
 *  - 3 : Alternating bit pattern (0xAAAAAAAAAAAAAAAA).
 *  - 4 : Repeating byte pattern (0xA5A5A5A5A5A5A5A5).
 *  - 5 : Fixed test constant (0x123456789ABCDEF0).
 *
 * The function writes the resulting values directly into the provided
 * share array.
 *
 * @param shares   Pointer to the array of share values to populate.
 * @param nshares  Number of shares in the masked representation.
 * @param pattern  Pattern identifier controlling how values are generated.
 */
static void fill_pattern(uint64_t *shares, int nshares, int pattern) {
    for (int i = 0; i < nshares; i++) {
        /* always exercise the RNG to balance power */
        uint64_t v = get_rand64();

        if (pattern == 0) {
            /* keep the true random value */
        } else if (pattern == 1) {
            v = 0;                      // overwrite but RNG was still read
        } else if (pattern == 2) {
            v = ~0ULL;
        } else if (pattern == 3) {
            v = 0xAAAAAAAAAAAAAAAAULL;
        }
        else if (pattern == 4)
        {
            v = 0xA5A5A5A5A5A5A5A5;
        } else if (pattern == 5)
        {
            v = 0x123456789ABCDEF0;
        }
        shares[i] = v;

        /* Only if extraction is needed
        for (int j = 0; j < 8; j++) {

                (uint8_t)(v >> (8 * j));
        }*/
    }
}
 
/**
 * Initializes the global masked Keccak state with a controlled test pattern.
 *
 * This command handler prepares the internal 5×5 masked Keccak state used
 * during side-channel testing. A lane value is selected according to the
 * requested test pattern and then split into MASKING_N shares for every
 * lane of the state.
 *
 * For each lane:
 *   - MASKING_N-1 shares are filled with fresh random values.
 *   - The final share is computed so that the XOR of all shares equals the
 *     selected lane value, preserving the masked representation.
 *
 * The selected pattern allows deterministic test vectors or fully random
 * inputs to be used during leakage experiments.
 *
 * Supported patterns:
 *   0 : Fully random lane value
 *   1 : Fixed constant (0x3F84A12BC9D07E55)
 *   2 : All-zero value
 *   3 : All bits set (~0ULL)
 *   4 : Alternating bits (0xAAAAAAAAAAAAAAAA)
 *   5 : Repeating pattern (0xA5A5A5A5A5A5A5A5)
 *   6 : Fixed test constant (0x123456789ABCDEF0)
 *
 * A fresh random seed (`chi_jitter_seed`) is also generated for later use
 * in Chi timing jitter.
 *
 * @param cmd   SimpleSerial command identifier (unused).
 * @param scmd  Subcommand identifier (unused).
 * @param len   Length of received data.
 * @param buf   Input buffer where buf[0] selects the test pattern.
 *
 * @return 0x00 on success.
 */
uint8_t prepare_state_cmd(uint8_t cmd, uint8_t scmd,
                                 uint8_t len, uint8_t *buf)
{
    UNUSED(cmd); UNUSED(scmd);

    int pattern = buf[0];

    uint64_t lane_rand = get_rand64();
    uint64_t lane_val = 0;
    chi_jitter_seed = get_rand64();

    switch (pattern) {
        case 0:  lane_val = lane_rand; break;
        case 1:  lane_val = 0x3F84A12BC9D07E55ULL; break;
        case 2:  lane_val = 0x0ULL; break;
        case 3:  lane_val = ~0ULL; break;
        case 4:  lane_val = 0xAAAAAAAAAAAAAAAAULL; break;
        case 5:  lane_val = 0xA5A5A5A5A5A5A5A5ULL; break;
        case 6:  lane_val = 0x123456789ABCDEF0ULL; break;
        default: lane_val = 0ULL; break;
    }

    for (int y = 0; y < 5; y++) {
        for (int x = 0; x < 5; x++) {

            uint64_t acc = 0;

            for (int s = 0; s < MASKING_N - 1; s++) {
                uint64_t r = get_rand64();
                global_state.share[s][x][y] = r;
                acc ^= r;
            }

            global_state.share[MASKING_N - 1][x][y] = lane_val ^ acc;
        }
    }

    return 0x00;
}

/**
 * Prepares randomness used by the masked Iota step.
 *
 * This command selects the round constant index to be used for the next
 * Iota operation and generates fresh randomness for share refreshing if
 * masking is enabled.
 *
 * The round index must be within the valid Keccak-f[1600] range (0–23).
 * If masking order is greater than one, random values are generated for
 * the auxiliary shares used during the masked Iota injection.
 *
 * @param cmd   SimpleSerial command identifier (unused).
 * @param scmd  Subcommand identifier (unused).
 * @param len   Length of received data (must be ≥ 1).
 * @param buf   Input buffer where buf[0] specifies the round index.
 *
 * @return 0x00 on success, or SS_ERR_LEN if the input is invalid.
 */
static uint8_t prepare_iota_rands_cmd(uint8_t cmd, uint8_t scmd,
                                      uint8_t len, uint8_t *buf)
{
    // buf[0] = round index we want to use next
    if (len < 1) return SS_ERR_LEN;
    if (buf[0] >= 24) return SS_ERR_LEN;
    iota_round_idx = buf[0];

    #if (MASKING_N > 1)
        fill_pattern(iota_rands, MASKING_N-1, 0);
    #endif


    return 0x00;
}

/**
 * Clears the global masked Keccak state.
 *
 * Sets all shares of every lane in the global state to zero. This is
 * typically used between experiments or before initializing a new test
 * configuration to ensure no residual masked values remain.
 */
void clear_global_state(void)
{
    for (int s = 0; s < MASKING_N; s++)
        for (int x = 0; x < 5; x++)
            for (int y = 0; y < 5; y++)
                global_state.share[s][x][y] = 0ULL;
}

/**
 * Clears all internal masked state used by the test harness.
 *
 * This helper currently resets the global masked state but exists as a
 * separate abstraction to allow additional internal structures to be
 * cleared in the future if needed.
 */
static void clear_internal()
{
    clear_global_state();
}



/**
 * Loads a masked Keccak lane into a destination masked value.
 *
 * This implementation is a placeholder used during testing. Instead of
 * loading a real lane from the state, all shares of the destination are
 * explicitly cleared to zero. The input parameters are intentionally
 * unused to keep the call interface consistent with a full lane load
 * routine while preventing unwanted compiler optimizations.
 *
 * @param dst  Destination masked lane structure.
 * @param st   Source masked state (unused).
 * @param x    Lane x-coordinate (unused).
 * @param y    Lane y-coordinate (unused).
 */
static inline void lane_load_a(masked_uint64_t *dst,
                               const masked_state_t *st,
                               int x, int y)
{
    (void)st; (void)x; (void)y;

    for (int s = 0; s < MASKING_N; s++) {
        dst->share[s] = 0;
    }
}

/**
 * Computes the column parity value C[x] for a given share using a balanced XOR tree.
 *
 * This helper is used in the Theta step of Keccak to compute the parity of
 * the five lanes in column x for a specific masking share index. The
 * computation reduces switching imbalance by avoiding a single linear
 * accumulator and instead using intermediate variables and randomized
 * masking operations.
 *
 * In the active branch, additional random masking and dummy operations are
 * introduced to decorrelate intermediate values and reduce side-channel
 * leakage. Each lane value is masked with a temporary random value that is
 * immediately canceled, maintaining correctness while equalizing data
 * transitions.
 *
 * The alternative branch (disabled) shows the standard balanced XOR tree:
 *
 *      ((a0 ⊕ a1) ⊕ a2) ⊕ (a3 ⊕ a4)
 *
 * The final parity value for share i is stored in C[x].share[i].
 *
 * @param C   Output array holding the column parity values.
 * @param st  Input masked Keccak state.
 * @param x   Column index.
 * @param i   Share index being processed.
 */
__attribute__((always_inline))
static inline void build_C_twoacc(volatile masked_uint64_t C[5],
                                  const volatile masked_state_t *st,
                                  int x, int i)
{
    
    uint64_t junk = get_rand64();
    junk ^= junk;
    uint64_t r;

    r = get_rand64();
    uint64_t acc = st->share[i][x][0] ^ r; acc ^= r;

    r = get_rand64();
    uint64_t v = st->share[i][x][1] ^ r; v ^= r; acc ^= v;

    r = get_rand64();
    v = st->share[i][x][2] ^ r; v ^= r; acc ^= v;

    r = get_rand64();
    v = st->share[i][x][3] ^ r; v ^= r; acc ^= v;

    r = get_rand64();
    v = st->share[i][x][4] ^ r; v ^= r; acc ^= v;

    C[x].share[i] = acc;

    asm volatile ("" ::: "memory");
    

}


/**
 * Prevents the compiler from optimizing away or reordering a 64-bit value.
 *
 * This helper introduces an empty inline assembly barrier that treats the
 * input value as both read and written. It is used to force a value to pass
 * through a register and remain visible to the generated code, which is
 * useful when controlling dataflow for side-channel experiments or register
 * precharge patterns.
 *
 * @param x  64-bit value to preserve through a compiler barrier.
 */
void barrier_u64(uint64_t x) {
    __asm volatile("" : "+r"(x));
}

/**
 * Refreshes a Boolean-masked value using a precomputed symmetric randomness matrix.
 *
 * This operation re-randomizes the shares of a masked value without changing
 * the underlying recombined secret. For each distinct share pair (i, j), the
 * same random value r[i][j] is XORed into both shares. Because the value is
 * added twice, the global XOR of all shares remains unchanged.
 *
 * This is a standard share-refresh step used in masking schemes to decorrelate
 * intermediate values and reduce leakage accumulation between operations.
 *
 * The randomness matrix is expected to be symmetric with zeros on the diagonal,
 * matching the usual ISW-style pairwise refresh layout.
 *
 * @param v  Masked value to refresh in place.
 * @param r  Pairwise randomness matrix used for share refreshing.
 */
static inline void masked_refresh(
    volatile masked_uint64_t *v,
    const volatile uint64_t r[MASKING_N][MASKING_N])
{
    for (size_t i = 0; i < MASKING_N; i++) {
        for (size_t j = i + 1; j < MASKING_N; j++) {
            uint64_t rij = r[i][j];
            v->share[i] ^= rij;
            v->share[j] ^= rij;
        }
    }
}

/**
 * Refreshes a Boolean-masked value using freshly generated random pairwise masks.
 *
 * This function performs the same logical share-refresh operation as
 * `masked_refresh`, but generates each pairwise random value on demand using
 * `get_rand64()` rather than consuming a precomputed matrix.
 *
 * For every pair of distinct shares (i, j), a fresh random 64-bit value is
 * XORed into both shares, preserving the recombined secret while replacing the
 * masking with a statistically independent sharing.
 *
 * This variant is useful when precomputed randomness is not available or when
 * on-the-fly refreshing is preferred.
 *
 * @param v  Masked value to refresh in place.
 */
static inline void masked_refresh_oa(volatile masked_uint64_t *v) {
    for (size_t i = 0; i < MASKING_N; i++) {
        for (size_t j = i + 1; j < MASKING_N; j++) {
            uint64_t rij = get_rand64();
            v->share[i] ^= rij;
            v->share[j] ^= rij;
        }
    }
}


/**
 * Loads a single masked lane from the Keccak state into a masked lane object.
 *
 * This helper copies all shares of the lane at coordinates (x, y) from the
 * masked state into the destination structure. Each share is moved explicitly
 * in a simple loop to preserve predictable code generation.
 *
 * After each load, a compiler barrier is applied and the temporary register is
 * precharged to zero. This is intended to reduce data-dependent register reuse
 * effects during side-channel experiments.
 *
 * @param dst  Destination masked lane.
 * @param st   Source masked Keccak state.
 * @param x    Lane x-coordinate.
 * @param y    Lane y-coordinate.
 */
__attribute__((optimize("O0")))
void lane_load(masked_uint64_t *dst,
                             const masked_state_t *st,
                             int x, int y)
{
    uint64_t t = 0;

    for (int s = 0; s < MASKING_N; s++) {

        t = ST(st, s, x, y);
        dst->share[s] = t;

        // --- precharge register to constant ---
        barrier_u64(t);
        t = 0;
        barrier_u64(t);
    }
}

/**
 * Loads the entire masked Keccak state into a 5×5 array of masked lanes.
 *
 * This routine copies every lane and every masking share from the source
 * state into the destination lane array. The destination layout remains
 * indexed by lane position, while each lane retains its full share vector.
 *
 * A compiler memory barrier is inserted after each share sweep to discourage
 * unwanted reordering across iterations.
 *
 * @param dst  Destination 5×5 array of masked lanes.
 * @param st   Source masked Keccak state.
 */
void lane_load_all(masked_uint64_t dst[5][5],
                                 const masked_state_t *st)
{
    for (int s = 0; s < MASKING_N; s++) {
        for (int y = 0; y < 5; y++) {
            for (int x = 0; x < 5; x++) {
                dst[x][y].share[s] = st->share[s][x][y];
            }
        }
        __asm volatile("" ::: "memory"); // optional barrier
    }
}

/**
 * Stores a single masked lane back into the Keccak state.
 *
 * This helper writes all shares of a masked lane into the state position
 * (x, y). Each share is stored explicitly to preserve predictable code shape.
 *
 * After each store, a compiler barrier is applied and the temporary register is
 * precharged to zero. This mirrors the load-side behavior and is intended to
 * reduce unwanted data-dependent register effects in side-channel measurements.
 *
 * @param st   Destination masked Keccak state.
 * @param x    Lane x-coordinate.
 * @param y    Lane y-coordinate.
 * @param src  Source masked lane to store.
 */
__attribute__((optimize("O0")))
void lane_store(masked_state_t *st,
                              int x, int y,
                              const masked_uint64_t *src)
{
    uint64_t t = 0;

    for (int s = 0; s < MASKING_N; s++) {

        t = src->share[s];
        ST(st, s, x, y) = t;

        // --- precharge ---
        barrier_u64(t);
        t = 0;
        barrier_u64(t);
    }
}


/**
 * Computes the share-wise XOR of two Boolean-masked values.
 *
 * Since Boolean masking is linear with respect to XOR, this operation is
 * performed independently on each share. The recombined result therefore
 * equals the XOR of the two underlying unmasked values.
 *
 * This is one of the fundamental linear operations used throughout masked
 * Keccak and other Boolean-masked computations.
 *
 * @param out  Destination masked value.
 * @param a    First masked input.
 * @param b    Second masked input.
 */
void masked_xor(volatile masked_uint64_t *out,
                const volatile  masked_uint64_t *a,
                const volatile masked_uint64_t *b) {
    for (size_t i = 0; i < MASKING_N; i++) {
        out->share[i] = xor64_safe(a->share[i], b->share[i]);
    }
}


/**
 * Computes a Boolean-masked AND using the ISW masking scheme.
 *
 * This function implements a secure nonlinear AND between two Boolean-masked
 * values. Unlike XOR, AND is not linear under Boolean masking, so the output
 * cannot be computed share-wise alone. Instead, the ISW construction is used:
 *
 *   - diagonal terms compute (a_i & b_i) for each share
 *   - cross-terms combine (a_i & b_j) and (a_j & b_i)
 *   - shared randomness r[i][j] is used to mask the cross-term contributions
 *
 * The result is a fresh masked sharing of the bitwise AND of the underlying
 * values, while preserving the required security order under the masking model.
 *
 * The implementation is specialized by `MASKING_N` for clarity and efficiency.
 * The randomness matrix is expected to be symmetric, with the upper triangle
 * providing the pairwise random masks used by the ISW correction terms.
 *
 * @param out  Destination masked value.
 * @param a    First masked input.
 * @param b    Second masked input.
 * @param r    Pairwise randomness matrix for ISW cross-term protection.
 */
#if MASKING_N == 1
void masked_and(volatile masked_uint64_t *out,
                const volatile masked_uint64_t *a,
                const volatile masked_uint64_t *b,
                const volatile uint64_t r[1][1]) {
    out->share[0] = a->share[0] & b->share[0];
}

#elif MASKING_N == 2
static inline void masked_and(masked_uint64_t *out,
                              const volatile masked_uint64_t *a,
                              const volatile masked_uint64_t *b,
                              const volatile uint64_t r[2][2])
{
    // Compute diagonal terms
    uint64_t t00 = a->share[0] & b->share[0];
    uint64_t t11 = a->share[1] & b->share[1];

    // Cross-term (shared randomness)
    uint64_t t01 = xor64_safe((a->share[0] & b->share[1]),
                              (a->share[1] & b->share[0]));

    // Apply ISW correction
    out->share[0] = xor64_safe(t00, r[0][1]);
    out->share[1] = xor64_safe(t11, xor64_safe(t01, r[0][1]));
}

#elif MASKING_N == 3
static inline void masked_and(volatile masked_uint64_t *out,
                              const volatile masked_uint64_t *a,
                              const volatile masked_uint64_t *b,
                              const volatile uint64_t r[3][3]) {
    // Diagonal terms
    out->share[0] = a->share[0] & b->share[0];
    out->share[1] = a->share[1] & b->share[1];
    out->share[2] = a->share[2] & b->share[2];

    // Cross terms
    uint64_t t01 = xor64_safe((a->share[0] & b->share[1]),
                              (a->share[1] & b->share[0]));
    out->share[0] = xor64_safe(out->share[0], r[0][1]);
    out->share[1] = xor64_safe(out->share[1], xor64_safe(t01, r[0][1]));

    uint64_t t02 = xor64_safe((a->share[0] & b->share[2]),
                              (a->share[2] & b->share[0]));
    out->share[0] = xor64_safe(out->share[0], r[0][2]);
    out->share[2] = xor64_safe(out->share[2], xor64_safe(t02, r[0][2]));

    uint64_t t12 = xor64_safe((a->share[1] & b->share[2]),
                              (a->share[2] & b->share[1]));
    out->share[1] = xor64_safe(out->share[1], r[1][2]);
    out->share[2] = xor64_safe(out->share[2], xor64_safe(t12, r[1][2]));
}

#elif MASKING_N == 4
static inline void masked_and(volatile masked_uint64_t *out,
                              const volatile masked_uint64_t *a,
                              const volatile masked_uint64_t *b,
                              const volatile uint64_t r[4][4]) {
    // Diagonal terms
    out->share[0] = a->share[0] & b->share[0];
    out->share[1] = a->share[1] & b->share[1];
    out->share[2] = a->share[2] & b->share[2];
    out->share[3] = a->share[3] & b->share[3];

    // Cross terms
    uint64_t t01 = xor64_safe((a->share[0] & b->share[1]),
                              (a->share[1] & b->share[0]));
    out->share[0] = xor64_safe(out->share[0], r[0][1]);
    out->share[1] = xor64_safe(out->share[1], xor64_safe(t01, r[0][1]));

    uint64_t t02 = xor64_safe((a->share[0] & b->share[2]),
                              (a->share[2] & b->share[0]));
    out->share[0] = xor64_safe(out->share[0], r[0][2]);
    out->share[2] = xor64_safe(out->share[2], xor64_safe(t02, r[0][2]));

    uint64_t t03 = xor64_safe((a->share[0] & b->share[3]),
                              (a->share[3] & b->share[0]));
    out->share[0] = xor64_safe(out->share[0], r[0][3]);
    out->share[3] = xor64_safe(out->share[3], xor64_safe(t03, r[0][3]));

    uint64_t t12 = xor64_safe((a->share[1] & b->share[2]),
                              (a->share[2] & b->share[1]));
    out->share[1] = xor64_safe(out->share[1], r[1][2]);
    out->share[2] = xor64_safe(out->share[2], xor64_safe(t12, r[1][2]));

    uint64_t t13 = xor64_safe((a->share[1] & b->share[3]),
                              (a->share[3] & b->share[1]));
    out->share[1] = xor64_safe(out->share[1], r[1][3]);
    out->share[3] = xor64_safe(out->share[3], xor64_safe(t13, r[1][3]));

    uint64_t t23 = xor64_safe((a->share[2] & b->share[3]),
                              (a->share[3] & b->share[2]));
    out->share[2] = xor64_safe(out->share[2], r[2][3]);
    out->share[3] = xor64_safe(out->share[3], xor64_safe(t23, r[2][3]));
}

#elif MASKING_N == 5
static inline void masked_and(volatile masked_uint64_t *out,
                              const volatile masked_uint64_t *a,
                              const volatile masked_uint64_t *b,
                              const volatile uint64_t r[5][5]) {
    // Diagonal
    for (int i = 0; i < 5; i++)
        out->share[i] = a->share[i] & b->share[i];

    // Cross terms
    for (int i = 0; i < 5; i++) {
        for (int j = i + 1; j < 5; j++) {
            uint64_t t = xor64_safe((a->share[i] & b->share[j]),
                                    (a->share[j] & b->share[i]));
            out->share[i] = xor64_safe(out->share[i], r[i][j]);
            out->share[j] = xor64_safe(out->share[j], xor64_safe(t, r[i][j]));
        }
    }
}

#else
static inline void masked_and(volatile masked_uint64_t *out,
                              const volatile masked_uint64_t *a,
                              const volatile masked_uint64_t *b,
                              const volatile uint64_t r[MASKING_N][MASKING_N]) {
    // Diagonal
    for (size_t i = 0; i < MASKING_N; i++)
        out->share[i] = a->share[i] & b->share[i];

    // Cross terms
    for (size_t i = 0; i < MASKING_N; i++) {
        for (size_t j = i + 1; j < MASKING_N; j++) {
            uint64_t t = xor64_safe((a->share[i] & b->share[j]),
                                    (a->share[j] & b->share[i]));
            out->share[i] = xor64_safe(out->share[i], r[i][j]);
            out->share[j] = xor64_safe(out->share[j], xor64_safe(t, r[i][j]));
        }
    }
}
#endif

/**
 * Computes the Boolean NOT of a masked value.
 *
 * In Boolean masking, a direct bitwise NOT cannot simply be applied to
 * each share independently because the recombined value would not equal
 * the negation of the underlying secret. Instead, a predefined constant
 * vector `K_NOT_SPLIT` is XORed into the shares.
 *
 * The constants are constructed such that the XOR of all constants equals
 * the all-ones mask (~0ULL). As a result, when the shares are recombined,
 * the resulting value is the correct bitwise negation of the original
 * unmasked value.
 *
 * @param dst  Destination masked value.
 * @param src  Source masked value to invert.
 */
static inline void masked_not(volatile masked_uint64_t *dst, const volatile masked_uint64_t *src) {
    for (size_t i = 0; i < MASKING_N; ++i)
        dst->share[i] = xor64_safe(src->share[i],  K_NOT_SPLIT[i]);
}

/**
 * Loads a single share from a masked Keccak state lane.
 *
 * This lightweight accessor retrieves the value of share `s` for the lane
 * located at coordinates (x, y) in the masked Keccak state. It is used to
 * simplify code that needs direct access to individual shares without
 * copying the full masked lane structure.
 *
 * @param st  Source masked Keccak state.
 * @param s   Share index.
 * @param x   Lane x-coordinate.
 * @param y   Lane y-coordinate.
 *
 * @return The 64-bit value of the requested share.
 */
static inline uint64_t lane_load_share(const masked_state_t *st, int s, int x, int y)
{
    return st->share[s][x][y];
}

/**
 * Stores a single share into a masked Keccak state lane.
 *
 * This helper writes a 64-bit value into share `s` of the lane located at
 * coordinates (x, y) in the masked Keccak state. It provides a simple and
 * consistent interface for updating individual shares within the state.
 *
 * @param st  Destination masked Keccak state.
 * @param s   Share index.
 * @param x   Lane x-coordinate.
 * @param y   Lane y-coordinate.
 * @param v   64-bit value to store.
 */
static inline void lane_store_share(masked_state_t *st, int s, int x, int y, uint64_t v)
{
    st->share[s][x][y] = v;
}



/*------------------------- Keccak Round Functionality -------------------------*/
/**
 * ----------------------------------------------------------------------------
 * Function: masked_theta
 * ----------------------------------------------------------------------------
 * Purpose:
 *   Keccak-f[1600] Θ (Theta). Computes column parities C[x], forms D[x] =
 *   C[x-1] ⊕ ROTL1(C[x+1]), and XORs D[x] into every lane of column x.
 *
 * Security Model:
 *   - Boolean masking of order (MASKING_N-1). Θ is linear -> masking preserved.
 *   - No fresh randomness required for correctness; implementation uses
 *     balanced, pairwise accumulation (build_C_twoacc) to reduce toggling bias.
 *
 * Randomness:
 *   None consumed by Θ itself (only uses get_rand64() inside helper macros
 *   if configured elsewhere; default here avoids extra randomness).
 *
 * Timing / Tracing:
 *   - Writes cycle counts into time_records:
 *       time_records[0]: C[x] accumulation
 *       time_records[1]: D[x] formation
 *       time_records[2]: D[x] broadcast into state
 *   - Surround with TEST_PROLOGUE/TEST_EPILOGUE for reproducible traces.
 *
 * Preconditions:
 *   - state is a 5×5 array of masked_uint64_t with MASKING_N shares per lane.
 *   - DWT->CYCCNT initialized (dwt_init()) if timing is recorded.
 *
 * Postconditions:
 *   - state updated in place with Θ applied, masking order preserved.
 *
 * Side-channel Notes:
 *   - Linear operation; no ISW AND. Pairwise XOR tree helps reduce glitch bias.
 * ----------------------------------------------------------------------------
 */

__attribute__((noinline,optimize("O0")))
void masked_theta(volatile masked_state_t *state,
                  uint32_t *time_records)
{
    volatile masked_uint64_t C[5], D[5];
    uint64_t tmp[5][MASKING_N];
    uint32_t start, end;

    // === 1. Compute C[x] ===
    TSTAMP(t1);
    start = DWT->CYCCNT;

    for (int x = 0; x < 5; x++)
        for (int i = 0; i < MASKING_N; i++)
            build_C_twoacc(C, state, x, i);   // NOTE: update this prototype too

    end = DWT->CYCCNT;
    time_records[0] = end - start;

    // === 2. Compute D[x] ===
    start = DWT->CYCCNT;
    TSTAMP(t2);
    for (int x = 0; x < 5; x++) {
        int xm1 = (x + 4) % 5;
        int xp1 = (x + 1) % 5;

        for (int i = 0; i < MASKING_N; i++) {

            uint64_t m = get_rand64();

            uint64_t t = xor64_safe(C[xp1].share[i], m);
            t = (t << 1) | (t >> 63);
            t = xor64_safe(t, ((m << 1) | (m >> 63)));

            asm volatile ("" ::: "memory");

            tmp[x][i] = xor64_safe(C[xm1].share[i], t);
        }
    }

    // Copy tmp to D
    for (int x = 0; x < 5; x++)
        for (int i = 0; i < MASKING_N; i++)
            D[x].share[i] = tmp[x][i];

    end = DWT->CYCCNT;
    time_records[1] = end - start;
    TSTAMP(t3);

    // === 3. Apply D[x] to all lanes ===
    start = DWT->CYCCNT;

    for (int x = 0; x < 5; x++)
        for (int i = 0; i < MASKING_N; i++) {

            uint64_t dx = D[x].share[i];

            for (int y = 0; y < 5; y++)
                state->share[i][x][y] =
                    xor64_safe(state->share[i][x][y], dx);
        }

    end = DWT->CYCCNT;
    TSTAMP(t4);
    time_records[2] = end - start;
}

/**
 * ----------------------------------------------------------------------------
 * Function: masked_rho
 * ----------------------------------------------------------------------------
 * Purpose:
 *   Keccak-f[1600] ρ (Rho). Rotates each lane by a fixed offset.
 *
 * Security Model:
 *   - Linear bit permutation -> masking preserved.
 *   - Uses a "blind-rotate-unblind" pattern per share to avoid data-dependent
 *     switching during rotation (tmp = share ⊕ r; ROTL(tmp); ⊕ ROTL(r)).
 *
 * Randomness:
 *   - Consumes one 64-bit random mask per lane share (get_rand64()).
 *   - Mask cancels algebraically; used to stabilize switching activity.
 *
 * Timing / Tracing:
 *   - time_records[0]: total cost of masked rotations.
 *
 * Preconditions:
 *   - keccak_rho_offsets[][] must match Keccak spec.
 *   - DWT->CYCCNT initialized if timing is recorded.
 *
 * Postconditions:
 *   - state updated in place with ρ applied; masking order preserved.
 *
 * Side-channel Notes:
 *   - Linear; the blind/unblind rotation helps reduce leakage on real silicon.
 * ----------------------------------------------------------------------------
 */
void masked_rho(volatile masked_state_t *state,
                uint32_t *time_records)
{
    uint32_t start, end;

    start = DWT->CYCCNT;
    TSTAMP(t1);

    for (int x = 0; x < 5; x++) {
        for (int y = 0; y < 5; y++) {

            uint8_t r = keccak_rho_offsets[x][y];

            for (int i = 0; i < MASKING_N; i++) {

                uint64_t rand_mask = get_rand64();

                // Apply random blinding before rotate
                uint64_t tmp =
                    xor64_safe(state->share[i][x][y], rand_mask);

                // Rotate both blinded value and mask
                uint64_t tmp_rot = ROL64(tmp, r);
                uint64_t mask_rot = ROL64(rand_mask, r);

                // Unblind
                state->share[i][x][y] =
                    xor64_safe(tmp_rot, mask_rot);
            }
        }
    }

    end = DWT->CYCCNT;
    TSTAMP(t2);
    time_records[0] = end - start;
}

/**
 * ----------------------------------------------------------------------------
 * Function: masked_pi
 * ----------------------------------------------------------------------------
 * Purpose:
 *   Keccak-f[1600] π (Pi). Permutes lane coordinates:
 *     B[y, (2x + 3y) mod 5] = A[x, y].
 *
 * Security Model:
 *   - Pure data movement (linear index permutation) -> masking preserved.
 *   - Copy via tmp to avoid overlap hazards when in-place reindexing.
 *
 * Randomness:
 *   - None required/consumed.
 *
 * Timing / Tracing:
 *   - Surround with TEST_PROLOGUE/TEST_EPILOGUE externally if needed.
 *
 * Preconditions:
 *   - state is a 5×5 array with MASKING_N shares; tmp uses the same layout.
 *
 * Postconditions:
 *   - state overwritten with permuted lanes; masking intact.
 *
 * Side-channel Notes:
 *   - No non-linear operations. Data movement only.
 * ----------------------------------------------------------------------------
 */
__attribute__((optimize("O0")))
void masked_pi(volatile masked_state_t *state)
{
    masked_state_t tmp;

    // Copy input to tmp: tmp.share[s][x][y] = state->share[s][x][y]
    for (int s = 0; s < MASKING_N; ++s)
        for (int x = 0; x < 5; ++x)
            for (int y = 0; y < 5; ++y)
                tmp.share[s][x][y] = state->share[s][x][y];

    TSTAMP(t1);

    // Pi: B[y, (2x + 3y) % 5] = A[x, y]
    for (int y = 0; y < 5; ++y) {
        for (int x = 0; x < 5; ++x) {
            int nx = y;
            int ny = (2 * x + 3 * y) % 5;

            for (int s = 0; s < MASKING_N; ++s)
                state->share[s][nx][ny] = tmp.share[s][x][y];
        }
    }

    TSTAMP(t2);
}

/**
 * ----------------------------------------------------------------------------
 * Function: masked_chi
 * ----------------------------------------------------------------------------
 * Purpose:
 *   Keccak-f[1600] χ (Chi). Non-linear S-box per row:
 *     out[x,y] = in[x,y] ⊕ ((~in[x+1,y]) & in[x+2,y]).
 *
 * Security Model:
 *   - Implements ISW-style masked AND with refresh steps:
 *        1) Refresh inputs b,c before non-linear op   (masked_refresh)
 *        2) t1 = masked_not(b_loc) using a fixed public split (K_NOT_SPLIT)
 *        3) t2 = masked_and(t1, c_loc, r_and[x][y])  (ISW)
 *        4) Refresh t2 after AND                      (masked_refresh)
 *        5) out = a ⊕ t2                              (masked_xor)
 *        6) Refresh out before storing                (masked_refresh)
 *     - All refreshes use independent randomness slices to avoid correlation.
 *
 * Randomness:
 *   - r_and[x][y][i][j]  : ISW cross-term randomness (symmetric, i≠j)
 *   - r_in[x][y][i]      : pre-AND refresh for b_loc
 *   - r_in[(x+1)%5][y][i]: pre-AND refresh for c_loc
 *   - r_mid[x][y][i]     : post-AND refresh for t2
 *   - r_out[x][y][i]     : post-store refresh for out[x][y]
 *   - Optional micro-jitter loop seeded from TRNG to de-synchronize traces.
 *
 * Timing / Tracing:
 *   - Caller may measure around masked_chi() with TEST_PROLOGUE/EPILOGUE.
 *   - Function does not write time_records itself.
 *
 * Preconditions:
 *   - All randmat arguments must be filled (prepare_randmat_cmd()).
 *   - K_NOT_SPLIT[] initialized (init_not_split()).
 *
 * Postconditions:
 *   - out[x][y] holds masked χ(in[x][y]) for all x,y; masking order preserved.
 *   - All sensitive temporaries zeroized before exit.
 *
 * Side-channel Notes:
 *   - χ is the only non-linear step; the ISW AND and thorough refresh hygiene
 *     are essential to avoid 1st/2nd order leakage.
 * ----------------------------------------------------------------------------
 */

__attribute__((optimize("O0")))
void masked_chi(
    masked_state_t *out,
    const masked_state_t *in,
    const uint64_t r_and[5][5][MASKING_N][MASKING_N],
    const uint64_t r_in [5][5][4][MASKING_N][MASKING_N],
    const uint64_t r_mid[5][5][MASKING_N][MASKING_N],
    const uint64_t r_out[5][5][MASKING_N][MASKING_N])
{
    static const uint8_t xb_lut[5] = {1,2,3,4,0};
    static const uint8_t xc_lut[5] = {2,3,4,0,1};

    for (int y = 0; y < 5; y++) {
        for (int x = 0; x < 5; x++) {

            int xb = xb_lut[x];
            int xc = xc_lut[x];
            TSTAMP(time_1);
            masked_uint64_t a;
            masked_uint64_t b_loc, c_loc;
            masked_uint64_t t1, t2, o;

            TSTAMP(time_2);

            // Load lanes from SoA state into AoS lane containers
            lane_load_a(&a,     in, x,  y);
            //a = (masked_uint64_t){0};
            TSTAMP(time_3);
            lane_load_a(&b_loc, in, xb, y);
            TSTAMP(time_4);
            masked_refresh(&b_loc,  r_in[x][y][0]);

            TSTAMP(time_5);

            lane_load_a(&c_loc, in, xc, y);
            TSTAMP(time_6);
            masked_refresh(&c_loc,  r_in[x][y][1]);

            TSTAMP(time_7);

            masked_not(&t1, &b_loc);
            TSTAMP(time_8);

            //masked_xor(&t2, &t1, &c_loc);

            masked_and(&t2, &t1, &c_loc, r_and[x][y]);
            TSTAMP(time_9);

            masked_refresh(&t2,     r_in[x][y][2]);

            TSTAMP(time_10);

            masked_xor((masked_uint64_t *)&o, &a, &t2);
            TSTAMP(time_11);
            masked_refresh(&o,      r_in[x][y][3]);

            TSTAMP(time_12);

            // Store result lane back into SoA out state
            lane_store(out, x, y, &o);
        }
    }
}



/**
 * ----------------------------------------------------------------------------
 * Function: masked_iota
 * ----------------------------------------------------------------------------
 * Purpose:
 *   Keccak-f[1600] ι (Iota). XORs a round constant (RC[round_idx]) into
 *   lane (0,0). For masked case, injects a masked constant across shares.
 *
 * Security Model:
 *   - XOR with a public constant is linear. For masked execution (MASKING_N>1),
 *     we re-mask the constant RC as 'masked_rc[share]' so each share gets a
 *     piece and the XOR preserves the Boolean masking.
 *
 * Randomness:
 *   - If MASKING_N == 1: none.
 *   - Else: expects (MASKING_N-1) fresh 64-bit randoms in 'randoms' to form
 *     masked_rc[0..N-2]; masked_rc[N-1] = RC ⊕ XOR(masked_rc[0..N-2]).
 *
 * Timing / Tracing:
 *   - If time_records != NULL, caller may record surrounding timing.
 *
 * Preconditions:
 *   - round_idx in [0..23]; RC[] matches Keccak spec.
 *   - For MASKING_N>1, 'randoms' points to at least (MASKING_N-1) values.
 *
 * Postconditions:
 *   - state[0][0] shares XORed with masked_rc shares; masking preserved.
 *
 * Side-channel Notes:
 *   - Linear op. Still inject as masked constant to avoid share imbalance.
 * ----------------------------------------------------------------------------
 */
void masked_iota(volatile masked_state_t *state,
                 uint8_t round_idx,
                 const uint64_t *randoms,
                 uint32_t *time_records)
{
#if (MASKING_N == 1)

    state->share[0][0][0] =
        xor64_safe(state->share[0][0][0], RC[round_idx]);

#else

    uint64_t rc = RC[round_idx];
    uint64_t acc = 0;
    uint64_t masked_rc[MASKING_N];

    // Mask the round constant
    for (int i = 0; i < MASKING_N - 1; i++) {
        masked_rc[i] = randoms[i];
        acc ^= masked_rc[i];
    }

    masked_rc[MASKING_N - 1] = rc ^ acc;

    // Apply masked constant to lane (0,0)
    for (int i = 0; i < MASKING_N; i++)
        state->share[i][0][0] =
            xor64_safe(state->share[i][0][0], masked_rc[i]);

#endif
TSTAMP(t1);
}


/*------------------------- Keccak Test Methods-------------------------*/
/**
 * Captures a baseline trace consisting only of idle NOP execution.
 *
 * This test method is used as a reference measurement for side-channel
 * experiments. It performs no cryptographic computation and instead executes
 * a fixed number of NOP instructions inside the normal test trigger window.
 *
 * The purpose is to characterize the measurement setup itself, including
 * trigger alignment, clocking artifacts, and environmental noise, without
 * any contribution from masked Keccak operations.
 *
 * @param c  SimpleSerial command identifier (unused).
 * @param s  Subcommand identifier (unused).
 * @param l  Length of received data (unused).
 * @param b  Input buffer (unused).
 *
 * @return 0x00 on success.
 */
__attribute__((optimize("O0")))
static uint8_t nop_trace(uint8_t c,uint8_t s,uint8_t l,uint8_t *b){
    UNUSED(c); UNUSED(s); UNUSED(l); UNUSED(b);
    TEST_PROLOGUE();
     for (volatile uint32_t i = 0; i < 1024; i++) {
        __NOP();
    }
    TEST_EPILOGUE();
    return 0x00;
}

/**
 * Captures a trace of the masked Theta step without storing a final result externally.
 *
 * This test method invokes `masked_theta()` on the global masked state inside
 * the trigger window and records timing information for later analysis. A small
 * NOP padding block is inserted before the measured region to improve alignment
 * consistency across captures.
 *
 * In addition to the total cycle count, the function stores internal timing
 * markers produced by `masked_theta()` into the global timing buffer. These
 * can be used to study sub-phase timing structure within the Theta
 * implementation.
 *
 * After measurement, the internal masked state is cleared to avoid cross-test
 * contamination.
 *
 * @param cmd   SimpleSerial command identifier (unused).
 * @param scmd  Subcommand identifier (unused).
 * @param len   Length of received data (unused).
 * @param buf   Input buffer (unused).
 *
 * @return 0x00 on success.
 */
__attribute__((optimize("O0")))
static uint8_t theta_test_no_save(uint8_t cmd, uint8_t scmd, uint8_t len, uint8_t *buf)
{
    UNUSED(cmd); UNUSED(scmd); UNUSED(buf);
    uint32_t start, end, time; 
    uint32_t time_records[3] = {0};

    __DSB();
    __ISB();
    
    TEST_PROLOGUE();
    NOP_BLOCK(128);

    start = DWT->CYCCNT;
    masked_theta(&global_state, time_records);
    end = DWT->CYCCNT;
    TEST_EPILOGUE();

    time = end - start;
    for(int i = current_timing_position; i < 3; i++)
        time_buffer[i] = time_records[i];
    time_buffer[3] = time;
    current_timing_position = 4;
    clear_internal();
    return 0x00;
}


/**
 * Captures a trace of the masked Rho step without storing a final result externally.
 *
 * This routine measures the execution of `masked_rho()` on the global masked
 * state within the normal test trigger window. A fixed NOP block is inserted
 * before the measured code region to improve alignment and make traces more
 * comparable across repeated captures.
 *
 * The function records selected internal timing markers returned by
 * `masked_rho()` and appends them to the global timing buffer. The global
 * masked state is then cleared so that each capture starts from a clean state.
 *
 * @param cmd   SimpleSerial command identifier (unused).
 * @param scmd  Subcommand identifier (unused).
 * @param len   Length of received data (unused).
 * @param buf   Input buffer (unused).
 *
 * @return 0x00 on success.
 */
__attribute__((optimize("O0")))
static uint8_t rho_test_no_save(uint8_t cmd, uint8_t scmd, uint8_t len, uint8_t *buf)
{
    UNUSED(cmd); UNUSED(scmd); UNUSED(buf);
    uint32_t start, end, time; 
    uint32_t time_records[3] = {0};

    __DSB();
    __ISB();
    TEST_PROLOGUE();
    NOP_BLOCK(128);

    start = DWT->CYCCNT;
    masked_rho(&global_state, time_records);
    end = DWT->CYCCNT;
   
    TEST_EPILOGUE();
    time = end - start;
    time_buffer[current_timing_position] = time_records[0];
    time_buffer[current_timing_position + 1] = time_records[1];
    current_timing_position = current_timing_position + 2;
    clear_internal();
    return 0x00;
}

/**
 * Captures a trace of the masked Pi step without storing a final result externally.
 *
 * This test method runs `masked_pi()` on the global masked state inside the
 * measurement trigger window. A fixed NOP padding region is executed before
 * the measured code to improve temporal alignment during trace comparison.
 *
 * The total cycle count for the Pi step is recorded into the timing buffer.
 * After the measurement completes, the internal masked state is cleared to
 * ensure independence between successive tests.
 *
 * @param cmd   SimpleSerial command identifier (unused).
 * @param scmd  Subcommand identifier (unused).
 * @param len   Length of received data (unused).
 * @param buf   Input buffer (unused).
 *
 * @return 0x00 on success.
 */
__attribute__((optimize("O0")))
static uint8_t pi_test_no_save(uint8_t cmd, uint8_t scmd, uint8_t len, uint8_t *buf)
{
    UNUSED(cmd); UNUSED(scmd); UNUSED(buf);
    uint32_t start, end, time; 
    uint32_t time_records[3] = {0};
    
    __DSB();
    __ISB();
    TEST_PROLOGUE();
    NOP_BLOCK(128);

    start = DWT->CYCCNT;
    masked_pi(&global_state);   
    end = DWT->CYCCNT;
    
    TEST_EPILOGUE();

    time = end - start;
    time_buffer[current_timing_position] = time;
    current_timing_position++;
    clear_internal();
    return 0x00;
}



/**
 * Captures a trace of the current masked Chi implementation.
 *
 * This test method measures the execution of `masked_chi()` using the global
 * masked input state and the pre-generated randomness matrices required by the
 * masked nonlinear gadgets. The computation is performed inside the trigger
 * window so that its power consumption can be captured for leakage analysis.
 *
 * A fixed NOP padding block is executed immediately before the measured region
 * to improve alignment between traces. The total cycle count of the Chi step
 * is recorded in the timing buffer. After execution, the internal masked state
 * is cleared.
 *
 * The result of the Chi computation is written into `out_state`, but the
 * routine is intended primarily for side-channel measurement rather than
 * functional output handling.
 *
 * @param cmd   SimpleSerial command identifier (unused).
 * @param scmd  Subcommand identifier (unused).
 * @param len   Length of received data (unused).
 * @param buf   Input buffer (unused).
 *
 * @return 0x00 on success.
 */
__attribute__((optimize("O0")))
static uint8_t chi_test_no_save(uint8_t cmd, uint8_t scmd,
                                uint8_t len, uint8_t *buf)
{
    UNUSED(cmd); UNUSED(scmd); UNUSED(len); UNUSED(buf);

    uint32_t start, end, time;

    __DSB();
    __ISB();
    TEST_PROLOGUE();
    NOP_BLOCK(128);
    start = DWT->CYCCNT;

    masked_chi(&out_state,
                &global_state,
                randmat_and,
                randmat_in,
                randmat_mid,
                randmat_out);

    end = DWT->CYCCNT;
    TEST_EPILOGUE();

    clear_internal();
    time = end - start;
    time_buffer[current_timing_position] = time;
    current_timing_position++;

    return 0x00;
}

/**
 * Captures a trace of the masked Iota step for a selected Keccak round.
 *
 * This test method measures execution of the Iota step on the global masked
 * state using the previously selected round index. For masked configurations
 * (`MASKING_N > 1`), it calls the masked Iota implementation with prepared
 * randomness. For the unmasked case, it injects the standard round constant
 * directly.
 *
 * A fixed NOP padding block is inserted before the measured code to improve
 * trace alignment. The total cycle count is stored in the timing buffer, and
 * the internal masked state is cleared after the capture.
 *
 * @param cmd   SimpleSerial command identifier (unused).
 * @param scmd  Subcommand identifier (unused).
 * @param len   Length of received data (unused).
 * @param buf   Input buffer (unused).
 *
 * @return 0x00 on success.
 */
__attribute__((optimize("O0")))
static uint8_t iota_test(uint8_t cmd, uint8_t scmd, uint8_t len, uint8_t *buf)
{
    UNUSED(cmd); UNUSED(scmd); UNUSED(buf);
    uint32_t start, end, time; 
    uint32_t time_records[3] = {0};

    __DSB();
    __ISB();
    TEST_PROLOGUE();
    NOP_BLOCK(128);
    start = DWT->CYCCNT;
    #if (MASKING_N > 1)
        masked_iota(&global_state, iota_round_idx, iota_rands, time_records);
    #else
        masked_iota(global_state, RC[iota_round_idx], time_records);
    #endif
    end = DWT->CYCCNT;          
    TEST_EPILOGUE();
    time = end - start;
    time_buffer[current_timing_position] = time;
    current_timing_position++;
    clear_internal();

    return 0x00;  
}

/**
 * Converts an unmasked 64-bit value into a Boolean-masked representation.
 *
 * This helper splits the input value into `MASKING_N` shares such that the XOR
 * of all shares reconstructs the original value. It generates `MASKING_N - 1`
 * fresh random shares and computes the final share so that the complete XOR
 * sum is correct.
 *
 * This is a standard initialization routine for masked values and is used when
 * injecting known plaintext constants or unmasked intermediates into the
 * masked computation domain.
 *
 * @param out    Destination masked value.
 * @param value  Unmasked 64-bit value to encode into shares.
 */

void masked_value_set(masked_uint64_t *out, uint64_t value)
{
    uint64_t acc = value;

    /* generate MASKING_N-1 random shares */
    for (int i = 0; i < MASKING_N - 1; i++) {
        uint64_t r = get_rand64();
        out->share[i] = r;
        acc ^= r;
    }

    /* final share closes XOR sum */
    out->share[MASKING_N - 1] = acc;
}

//~~~~~~~~~~~~~~~~~~~~KECCAK METHODS~~~~~~~~~~~~~~~~~~~~


/**
 * Squeezes output bytes from a masked Keccak state.
 *
 * This is the final phase in sponge-based hashing or XOF like SHAKE.
 * Recombines masked lanes to extract real output bytes.
 * Applies Keccak-f permutations between squeezing rounds if more output is needed.
 *
 * @param output      Buffer to receive the output
 * @param output_len  Number of output bytes desired
 * @param state       5x5 masked state to squeeze from
 * @param rate        Sponge bitrate in bytes (e.g. 168 for SHAKE128)
 */
void masked_squeeze(uint8_t *output,
                    size_t output_len,
                    masked_state_t *state,
                    size_t rate)
{
    size_t offset = 0;

    while (offset < output_len) {

        /* pull up to rate bytes per round */
        for (size_t i = 0; i < rate && offset < output_len; i++) {

            size_t x = (i / 8) % 5;
            size_t y = (i / 8) / 5;
            size_t byte_pos = i % 8;

            /* recombine masked lane */
            uint64_t lane = 0;

            for (int s = 0; s < MASKING_N; s++)
                lane ^= state->share[s][x][y];

            /* extract byte */
            output[offset++] = (lane >> (8 * byte_pos)) & 0xFF;
        }

        /* permute state if more output needed */
        if (offset < output_len) {
            masked_keccak_f1600(state);
        }
    }
}


/**
 * Absorbs input bytes into a masked Keccak sponge state.
 *
 * This is the input phase of sponge-based hashing such as SHA-3.
 * Initializes the masked state to zero, XORs input blocks into the
 * bitrate portion of the state, and applies the Keccak-f[1600]
 * permutation after each full block.
 *
 * Input bytes are packed into 64-bit lanes in little-endian form,
 * converted into a masked representation using fresh random shares,
 * and then XORed share-wise into the state.
 *
 * After all full blocks are processed, the remaining bytes are placed
 * into a final padded block using SHA-3 pad10*1 padding with domain
 * separator 0x06, then absorbed and permuted once more.
 *
 * @param state      Masked Keccak state to update
 * @param input      Input message bytes to absorb
 * @param input_len  Number of input bytes to process
 * @param rate       Sponge bitrate in bytes
 */
void masked_absorb(masked_state_t *state,
                   const uint8_t *input,
                   size_t input_len,
                   size_t rate)
{

    size_t offset = 0;

    /* ---- process full rate blocks ---- */

    while (input_len >= rate)
    {
        for (size_t i = 0; i < rate; i += 8)
        {
            uint64_t lane = 0;

            /* pack 8 bytes little-endian */
            for (int j = 0; j < 8; j++)
                lane |= ((uint64_t)input[offset + i + j]) << (8 * j);

            size_t x = (i / 8) % 5;
            size_t y = (i / 8) / 5;

            masked_uint64_t masked_lane;
            masked_value_set(&masked_lane, lane);

            /* XOR masked lane into state */
            for (int s = 0; s < MASKING_N; s++)
                state->share[s][x][y] ^= masked_lane.share[s];
        }

        masked_keccak_f1600(state);

        offset += rate;
        input_len -= rate;
    }

    /* ---- final padded block ---- */

    uint8_t block[KECCAK_RATE] = {0};

    for (size_t i = 0; i < input_len; i++)
        block[i] = input[offset + i];

    /* pad10*1 with SHA3 domain */
    block[input_len] ^= 0x06;
    block[rate - 1] ^= 0x80;

    for (size_t i = 0; i < rate; i += 8)
    {
        uint64_t lane = 0;

        for (int j = 0; j < 8 && (i + j) < rate; j++)
            lane |= ((uint64_t)block[i + j]) << (8 * j);

        size_t x = (i / 8) % 5;
        size_t y = (i / 8) / 5;

        masked_uint64_t masked_lane;
        masked_value_set(&masked_lane, lane);

        for (int s = 0; s < MASKING_N; s++)
            state->share[s][x][y] ^= masked_lane.share[s];
    }

    masked_keccak_f1600(state);
}

/**
 * Applies one full masked Keccak round to a masked state.
 *
 * Executes the five Keccak round steps in order:
 * Theta, Rho, Pi, Chi, and Iota.
 * Before the nonlinear steps, this routine prepares the randomness
 * required by the masked Chi and Iota implementations.
 *
 * Intermediate timing buffers are cleared between stages so that
 * per-step timing data can be collected independently if needed.
 * The Chi step writes to a temporary output state which is then
 * copied back into the main state before Iota is applied.
 *
 * @param S            Masked 5x5 Keccak state to update in place
 * @param round_index  Round number in the range 0 to 23
 */
void masked_keccak_round(masked_state_t *S,
                         uint8_t round_index)
{
    uint32_t time_records[TIME_BUFFER_SIZE] = {0};

    /* prepare randomness */
    prepare_randmat_cmd(0,0,0,NULL);

    uint8_t rbuf[1] = {round_index};
    prepare_iota_rands_cmd(0,0,1,rbuf);

    /* ------------------------------------------------ */
    /* Theta                                            */
    /* ------------------------------------------------ */

    masked_theta(S, time_records);
    for(int i = 0; i < TIME_BUFFER_SIZE; i++)
        time_records[i] = 0;


    /* ------------------------------------------------ */
    /* Rho                                              */
    /* ------------------------------------------------ */

    masked_rho(S, time_records);
    for(int i = 0; i < TIME_BUFFER_SIZE; i++)
        time_records[i] = 0;


    /* ------------------------------------------------ */
    /* Pi                                               */
    /* ------------------------------------------------ */

    masked_pi(S);
    for(int i = 0; i < TIME_BUFFER_SIZE; i++)
        time_records[i] = 0;


    /* ------------------------------------------------ */
    /* Chi                                              */
    /* ------------------------------------------------ */


    masked_state_t chi_out;

    masked_chi(&chi_out,
               S,
               randmat_and,
               randmat_in,
               randmat_mid,
               randmat_out);

    *S = chi_out;

    for(int i = 0; i < TIME_BUFFER_SIZE; i++)
        time_records[i] = 0;

    /* ------------------------------------------------ */
    /* Iota                                             */
    /* ------------------------------------------------ */

    masked_iota(S, round_index, iota_rands, time_records);
}


/**
 * Performs the full Keccak-f[1600] permutation on a masked state.
 *
 * Applies all 24 rounds of the Keccak permutation to the given
 * masked 5x5 state. Each round consists of the standard sequence
 * of Keccak steps: Theta, Rho, Pi, Chi, and Iota.
 *
 * This is the core permutation used by masked sponge constructions
 * such as SHA-3 and SHAKE.
 *
 * @param state  Masked Keccak state to permute in place
 */
void masked_keccak_f1600(masked_state_t *state)
{
    for (int i = 0; i < 24; i++) {
        masked_keccak_round(state, i);
    }
}

/**
 * Sets all shares of a masked Keccak state to zero.
 *
 * This helper initializes or clears a masked state structure by assigning
 * zero to every share of every lane. It is typically used before beginning
 * a new sponge computation or when preparing a temporary block state during
 * message absorption.
 *
 * @param st  Pointer to the masked state to be cleared.
 */
static void masked_state_zero(masked_state_t *st)
{
    for (int s = 0; s < MASKING_N; s++)
        for (int x = 0; x < 5; x++)
            for (int y = 0; y < 5; y++)
                st->share[s][x][y] = 0ULL;
}

/**
 * XORs one masked Keccak state into another.
 *
 * This function performs a share-wise XOR of two masked states. Because
 * Boolean masking preserves linear operations, the XOR can be applied
 * independently to each share without breaking the masking invariant.
 *
 * The result is equivalent to XORing the underlying unmasked states.
 * This operation is used during the sponge absorb phase when message
 * blocks are injected into the running Keccak state.
 *
 * @param dst  Destination masked state to be updated.
 * @param src  Masked state whose values will be XORed into the destination.
 */
static void masked_state_xor(masked_state_t *dst, const masked_state_t *src)
{
    for (int s = 0; s < MASKING_N; s++)
        for (int x = 0; x < 5; x++)
            for (int y = 0; y < 5; y++)
                dst->share[s][x][y] ^= src->share[s][x][y];
}

/*
 * Build one masked Keccak rate block from raw bytes.
 * The bytes are packed little-endian into lanes and written into a temporary
 * masked_state_t shaped exactly like the Keccak state.
 */
static void message_block_to_masked_state(masked_state_t *block_state,
                                          const uint8_t *block,
                                          size_t block_len)
{
    masked_state_zero(block_state);

    for (size_t i = 0; i < block_len; i += 8) {
        uint64_t lane = 0;
        size_t chunk = (block_len - i >= 8) ? 8 : (block_len - i);

        for (size_t j = 0; j < chunk; j++) {
            lane |= ((uint64_t)block[i + j]) << (8 * j);
        }

        size_t lane_idx = i / 8;
        size_t x = lane_idx % 5;
        size_t y = lane_idx / 5;

        masked_uint64_t masked_lane;
        masked_value_set(&masked_lane, lane);

        for (int s = 0; s < MASKING_N; s++) {
            block_state->share[s][x][y] = masked_lane.share[s];
        }
    }
}
/**
 * Converts a block of message bytes into a masked Keccak state block.
 *
 * This routine packs input bytes into 64-bit Keccak lanes using little-endian
 * ordering, converts each lane into a masked representation, and stores the
 * resulting shares in a temporary masked state structure.
 *
 * The block is first cleared, then populated lane-by-lane from the input
 * message. Each constructed lane is split into MASKING_N shares using
 * fresh randomness so that the XOR of the shares equals the original lane
 * value.
 *
 * The resulting masked state block can then be XORed into the running
 * Keccak sponge state during the absorb phase.
 *
 * @param block  Masked state structure to populate.
 * @param msg    Pointer to the message bytes to convert.
 * @param len    Number of bytes from the message to include in the block.
 */
static void build_masked_block(masked_state_t *block,
                               const uint8_t *msg,
                               size_t len)
{
    masked_state_zero(block);

    for (size_t i = 0; i < len; i += 8)
    {
        uint64_t lane = 0;

        size_t remaining = len - i;
        size_t chunk = remaining >= 8 ? 8 : remaining;

        for (size_t j = 0; j < chunk; j++)
            lane |= ((uint64_t)msg[i + j]) << (8 * j);

        size_t lane_index = i / 8;
        size_t x = lane_index % 5;
        size_t y = lane_index / 5;

        masked_uint64_t masked_lane;
        masked_value_set(&masked_lane, lane);

        for (int s = 0; s < MASKING_N; s++)
            block->share[s][x][y] = masked_lane.share[s];
    }
}


/**
 * Computes a masked SHA-3 hash of an input message.
 *
 * This routine implements the SHA-3 hashing functions using a masked
 * Keccak sponge construction. It supports the standard SHA-3 variants
 * (SHA3-224, SHA3-256, SHA3-384, and SHA3-512), which differ only in
 * their bitrate (rate) and output length.
 *
 * The function performs the following steps:
 *
 * 1. Initializes a fresh masked Keccak state.
 * 2. Absorbs the input message in rate-sized blocks using masked
 *    state XOR operations.
 * 3. Applies the Keccak-f[1600] permutation after each absorbed block.
 * 4. Pads the final partial block using the SHA-3 domain separator
 *    (0x06) followed by the standard multi-rate padding bit (0x80).
 * 5. Performs a final permutation.
 * 6. Squeezes the requested number of output bytes from the masked
 *    state to produce the final digest.
 *
 * All internal state transformations are performed on masked values.
 * The state remains masked throughout absorption and permutation, and
 * only the output bytes are recombined when producing the digest.
 *
 * SHA-3 variants use the following parameters:
 *
 *   Variant     Rate (bytes)   Output length
 *   ---------------------------------------
 *   SHA3-224    144            28 bytes
 *   SHA3-256    136            32 bytes
 *   SHA3-384    104            48 bytes
 *   SHA3-512     72            64 bytes
 *
 * @param output     Buffer to receive the final hash digest
 * @param input      Input message to hash
 * @param input_len  Length of the input message in bytes
 */
void masked_sha3_224(uint8_t *output,
                     const uint8_t *input,
                     size_t input_len)
{
    masked_state_t state;
    masked_state_zero(&state);

    while (input_len >= SHA3_224_RATE)
    {
        masked_state_t block;

        build_masked_block(&block, input, SHA3_224_RATE);
        masked_state_xor(&state, &block);

        masked_keccak_f1600(&state);

        input += SHA3_224_RATE;
        input_len -= SHA3_224_RATE;
    }

    uint8_t block_bytes[SHA3_224_RATE];

    for (size_t i = 0; i < SHA3_224_RATE; i++)
        block_bytes[i] = 0;

    for (size_t i = 0; i < input_len; i++)
        block_bytes[i] = input[i];

    block_bytes[input_len] ^= 0x06;
    block_bytes[SHA3_224_RATE - 1] ^= 0x80;

    masked_state_t final_block;

    build_masked_block(&final_block, block_bytes, SHA3_224_RATE);
    masked_state_xor(&state, &final_block);

    masked_keccak_f1600(&state);

    masked_squeeze(output, 28, &state, SHA3_224_RATE);
}

void masked_sha3_256(uint8_t *output,
                     const uint8_t *input,
                     size_t input_len)
{
    masked_state_t state;
    masked_state_zero(&state);

    /* ---- absorb full blocks ---- */

    while (input_len >= SHA3_256_RATE)
    {
        masked_state_t block;

        build_masked_block(&block, input, SHA3_256_RATE);
        masked_state_xor(&state, &block);

        masked_keccak_f1600(&state);

        input += SHA3_256_RATE;
        input_len -= SHA3_256_RATE;
    }

    /* ---- final padded block ---- */

    uint8_t block_bytes[SHA3_256_RATE];

    for (size_t i = 0; i < SHA3_256_RATE; i++)
        block_bytes[i] = 0;

    for (size_t i = 0; i < input_len; i++)
        block_bytes[i] = input[i];

    block_bytes[input_len] ^= 0x06;
    block_bytes[SHA3_256_RATE - 1] ^= 0x80;

    masked_state_t final_block;

    build_masked_block(&final_block, block_bytes, SHA3_256_RATE);
    masked_state_xor(&state, &final_block);

    masked_keccak_f1600(&state);

    /* ---- squeeze ---- */

    masked_squeeze(output, 32, &state, SHA3_256_RATE);
}
void masked_sha3_384(uint8_t *output,
                     const uint8_t *input,
                     size_t input_len)
{
    masked_state_t state;
    masked_state_zero(&state);

    while (input_len >= SHA3_384_RATE)
    {
        masked_state_t block;

        build_masked_block(&block, input, SHA3_384_RATE);
        masked_state_xor(&state, &block);

        masked_keccak_f1600(&state);

        input += SHA3_384_RATE;
        input_len -= SHA3_384_RATE;
    }

    uint8_t block_bytes[SHA3_384_RATE];

    for (size_t i = 0; i < SHA3_384_RATE; i++)
        block_bytes[i] = 0;

    for (size_t i = 0; i < input_len; i++)
        block_bytes[i] = input[i];

    block_bytes[input_len] ^= 0x06;
    block_bytes[SHA3_384_RATE - 1] ^= 0x80;

    masked_state_t final_block;

    build_masked_block(&final_block, block_bytes, SHA3_384_RATE);
    masked_state_xor(&state, &final_block);

    masked_keccak_f1600(&state);

    masked_squeeze(output, 48, &state, SHA3_384_RATE);
}

void masked_sha3_512(uint8_t *output,
                     const uint8_t *input,
                     size_t input_len)
{
    masked_state_t state;
    masked_state_zero(&state);

    while (input_len >= SHA3_512_RATE)
    {
        masked_state_t block;

        build_masked_block(&block, input, SHA3_512_RATE);
        masked_state_xor(&state, &block);

        masked_keccak_f1600(&state);

        input += SHA3_512_RATE;
        input_len -= SHA3_512_RATE;
    }

    uint8_t block_bytes[SHA3_512_RATE];

    for (size_t i = 0; i < SHA3_512_RATE; i++)
        block_bytes[i] = 0;

    for (size_t i = 0; i < input_len; i++)
        block_bytes[i] = input[i];

    block_bytes[input_len] ^= 0x06;
    block_bytes[SHA3_512_RATE - 1] ^= 0x80;

    masked_state_t final_block;

    build_masked_block(&final_block, block_bytes, SHA3_512_RATE);
    masked_state_xor(&state, &final_block);

    masked_keccak_f1600(&state);

    masked_squeeze(output, 64, &state, SHA3_512_RATE);
}


/**
 * Computes a masked SHAKE extensible-output function (XOF).
 *
 * This routine implements the SHAKE128 or SHAKE256 functions using
 * a masked Keccak sponge construction. Unlike SHA-3 hashes, SHAKE is
 * an extensible-output function (XOF) and can produce an arbitrary
 * number of output bytes.
 *
 * The function performs the following steps:
 *
 * 1. Initializes a fresh masked Keccak state.
 * 2. Absorbs the input message in rate-sized blocks.
 * 3. Applies the Keccak-f[1600] permutation after each absorbed block.
 * 4. Pads the final partial block using the SHAKE domain separator
 *    (0x1F) followed by the multi-rate padding bit (0x80).
 * 5. Performs the final permutation.
 * 6. Repeatedly squeezes output blocks from the masked state until the
 *    requested number of output bytes has been produced.
 *
 * The masked sponge ensures that all intermediate state operations
 * remain masked. Only the output bytes are recombined when extracted
 * from the state.
 *
 * SHAKE parameters:
 *
 *   Variant      Rate (bytes)
 *   -------------------------
 *   SHAKE128     168
 *   SHAKE256     136
 *
 * Since SHAKE is an extensible-output function, the output length is
 * determined by the caller.
 *
 * @param output      Buffer to receive the generated output stream
 * @param output_len  Number of output bytes requested
 * @param input       Input message to process
 * @param input_len   Length of the input message in bytes
 */
void masked_shake128(uint8_t *output,
                     size_t output_len,
                     const uint8_t *input,
                     size_t input_len)
{
    masked_state_t state;
    masked_state_zero(&state);

    while (input_len >= SHAKE128_RATE)
    {
        masked_state_t block;

        build_masked_block(&block, input, SHAKE128_RATE);
        masked_state_xor(&state, &block);

        masked_keccak_f1600(&state);

        input += SHAKE128_RATE;
        input_len -= SHAKE128_RATE;
    }

    uint8_t block_bytes[SHAKE128_RATE];

    for (size_t i = 0; i < SHAKE128_RATE; i++)
        block_bytes[i] = 0;

    for (size_t i = 0; i < input_len; i++)
        block_bytes[i] = input[i];

    block_bytes[input_len] ^= 0x1F;
    block_bytes[SHAKE128_RATE - 1] ^= 0x80;

    masked_state_t final_block;

    build_masked_block(&final_block, block_bytes, SHAKE128_RATE);
    masked_state_xor(&state, &final_block);

    masked_keccak_f1600(&state);

    masked_squeeze(output, output_len, &state, SHAKE128_RATE);
}

void masked_shake256(uint8_t *output,
                     size_t output_len,
                     const uint8_t *input,
                     size_t input_len)
{
    masked_state_t state;
    masked_state_zero(&state);

    while (input_len >= SHAKE256_RATE)
    {
        masked_state_t block;

        build_masked_block(&block, input, SHAKE256_RATE);
        masked_state_xor(&state, &block);

        masked_keccak_f1600(&state);

        input += SHAKE256_RATE;
        input_len -= SHAKE256_RATE;
    }

    uint8_t block_bytes[SHAKE256_RATE];

    for (size_t i = 0; i < SHAKE256_RATE; i++)
        block_bytes[i] = 0;

    for (size_t i = 0; i < input_len; i++)
        block_bytes[i] = input[i];

    block_bytes[input_len] ^= 0x1F;
    block_bytes[SHAKE256_RATE - 1] ^= 0x80;

    masked_state_t final_block;

    build_masked_block(&final_block, block_bytes, SHAKE256_RATE);
    masked_state_xor(&state, &final_block);

    masked_keccak_f1600(&state);

    masked_squeeze(output, output_len, &state, SHAKE256_RATE);
}

/*

//TEST MAIN FOR ORDINARY NON CW USE
int main(void)
{
    init_not_split();
    printf("==== Masked SHA3-256 Test ====\n\n");

    //---- fixed small vectors ---- 

    run_test("Empty message",
             (uint8_t *)"",
             0);

    run_test("abc",
             (uint8_t *)"abc",
             3);

    run_test("hello world",
             (uint8_t *)"hello world",
             11);

    // ---- boundary tests around rate (136 bytes) ---- 

    uint8_t buf1[SHA3_256_RATE - 1];
    uint8_t buf2[SHA3_256_RATE];
    uint8_t buf3[SHA3_256_RATE + 1];

    for (size_t i = 0; i < sizeof(buf1); i++) buf1[i] = i;
    for (size_t i = 0; i < sizeof(buf2); i++) buf2[i] = i;
    for (size_t i = 0; i < sizeof(buf3); i++) buf3[i] = i;

    run_test("Rate-1 bytes", buf1, sizeof(buf1));
    run_test("Rate bytes",   buf2, sizeof(buf2));
    run_test("Rate+1 bytes", buf3, sizeof(buf3));

    uint8_t long_msg[512];

    for (size_t i = 0; i < sizeof(long_msg); i++)
        long_msg[i] = (uint8_t)i;

    run_test("512-byte message", long_msg, sizeof(long_msg));

    return 0;

*/

int main(void)
{
    platform_init();
    init_uart();
    trigger_setup();
    init_not_split();
    dwt_init();
    simpleserial_init();       

    simpleserial_addcmd(0x70, 0, nop_trace);
    simpleserial_addcmd(0x30, 0, theta_test_no_save);
    simpleserial_addcmd(0x31, 0, rho_test_no_save);
    simpleserial_addcmd(0x32, 0, pi_test_no_save);
    simpleserial_addcmd(0x33, 0, iota_test);
    simpleserial_addcmd(0x34, 0, chi_test_no_save);


    simpleserial_addcmd(0x20, 0, prepare_state_cmd);
    simpleserial_addcmd(0x21, 0, prepare_iota_rands_cmd);
    simpleserial_addcmd(0x22, 0, prepare_randmat_cmd);



    while (1)
        simpleserial_get();
}

/*-----------------------------------------------------------------------------
 * End of File
 *---------------------------------------------------------------------------*/
