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
#include <string.h>
#include <stddef.h>
#include "stm32f4xx.h"      // <-- defines IRQn_Type, SysTick_IRQn, etc.
#include "core_cm4.h"       // <-- defines __get_PRIMASK(), __DSB(), etc.




#ifndef UNUSED
#define UNUSED(x) ((void)(x))
#endif

#ifndef MASKING_ORDER
#define MASKING_ORDER 1
#endif
#define SS_ERR_LEN 0xA1

#define MAX_ORDER 10 //arbitrary limit - this can be removed

#define NROUNDS 24 //Not used
#define MASKING_N_VAL (MASKING_ORDER + 1)
#define MASKING_N MASKING_N_VAL


// --- Exact-NOP macro using GAS .rept (no loop overhead) -----------------
#define NOP_BLOCK(N) asm volatile (".rept " #N "\n\tnop\n\t.endr\n" ::: "memory")


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
// --- set up macros to memory barriers and distinct trace
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

//  trigger_low();    //trigger_high();   
#define ROL64(x, n) (((x) << (n)) | ((x) >> (64 - (n))))

//Spaces in time record buffer
#define TIME_BUFFER_SIZE 32

//fundamental masked object
typedef struct {
    uint64_t share[MASKING_N];
} masked_uint64_t;


//global state for keccak 
static volatile masked_uint64_t global_state[5][5];

// Random matrices for masked_chi()
static volatile uint64_t randmat[5][5][MASKING_N][MASKING_N];
static volatile uint64_t randmat_and [5][5][MASKING_N][MASKING_N];
static volatile uint64_t randmat_in  [5][5][MASKING_N - 1];
static volatile uint64_t randmat_mid [5][5][MASKING_N - 1];
static volatile uint64_t randmat_out [5][5][MASKING_N - 1];

//for determining pause length in chi
volatile uint64_t chi_jitter_seed = 0;


// ----------------------------------------------------------------------------
// K_NOT_SPLIT[] holds a fixed public Boolean share decomposition of 0xFFFF...FFFF.
// Used by masked_not() to implement bitwise NOT securely under Boolean masking.
// Constructed so that XOR of all shares == 0xFFFFFFFFFFFFFFFFULL.
// ----------------------------------------------------------------------------
uint64_t K_NOT_SPLIT[MASKING_N];



//Stores timing records from the test and implementation functions
static volatile uint32_t time_buffer[TIME_BUFFER_SIZE] = {0};
//number of time captures recorded during run
static int time_captures = 0;

//Specific randomness for Iota
static volatile uint64_t iota_rands[MASKING_N-1];
static volatile uint8_t iota_round_idx = 0; 



static const uint8_t keccak_rho_offsets[5][5] = {
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


// ----------------------------------------------------------------------------
// The chip used for capture was not returning accurate 64 bit xor sums 
// This method splits each 64 bit xor into two 32 bit ones
// ----------------------------------------------------------------------------
static inline uint64_t xor64_safe(uint64_t a, uint64_t b) {
    uint32_t a_lo = (uint32_t)(a & 0xFFFFFFFFULL);
    uint32_t a_hi = (uint32_t)(a >> 32);
    uint32_t b_lo = (uint32_t)(b & 0xFFFFFFFFULL);
    uint32_t b_hi = (uint32_t)(b >> 32);

    uint64_t res_lo = (uint64_t)(a_lo ^ b_lo);
    uint64_t res_hi = (uint64_t)(a_hi ^ b_hi);

    return (res_hi << 32) | res_lo;
}

/**
 * ----------------------------------------------------------------------------
 * Function: get_rand64
 * ----------------------------------------------------------------------------
 * Purpose:
 *   Compose a 64-bit random value from two 32-bit HAL RNG draws.
 *
 * Security / Side-channel Notes:
 *   - Relies on the platform HAL `get_rand()` primitive.
 *   - Stateless combiner (hi||lo); no whitening beyond concatenation.
 *
 * Returns:
 *   64-bit random value.
 * ----------------------------------------------------------------------------
 */
//Calling RNG function defined in CW HAL
uint64_t get_rand64(void)
{
    uint64_t hi = get_rand();
    uint64_t lo = get_rand();
    return (hi << 32) | lo;
}

/**
 * ----------------------------------------------------------------------------
 * Function: init_not_split
 * ----------------------------------------------------------------------------
 * Purpose:
 *   Initialize K_NOT_SPLIT[] — a fixed, public Boolean decomposition of
 *   0xFFFFFFFFFFFFFFFF used by masked_not(). The XOR of all shares equals
 *   all-ones, enabling NOT via XOR with public shares.
 *
 * Construction:
 *   - First N-1 shares are derived from distinct, public constants.
 *   - Final share is chosen so that XOR of all shares == ~0ULL.
 *
 * Security / Side-channel Notes:
 *   - Public constants only (no secret/random inputs).
 *   - Used to avoid data-dependent negation; XOR-only implementation.
 * ----------------------------------------------------------------------------
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

/**
 * ----------------------------------------------------------------------------
 * Function: dwt_init
 * ----------------------------------------------------------------------------
 * Purpose:
 *   Enable and reset the DWT cycle counter (DWT->CYCCNT) for cycle-accurate
 *   timing of masked primitives.
 *
 * Preconditions:
 *   - Core debug block accessible on this device.
 *
 * Postconditions:
 *   - DWT->CYCCNT = 0 and counting.
 *
 * Side-channel Notes:
 *   - Timing is for instrumentation only; does not affect data path.
 * ----------------------------------------------------------------------------
 */

void dwt_init(void) {
    // Enable TRC (trace)
    CoreDebug->DEMCR |= CoreDebug_DEMCR_TRCENA_Msk;
    
    // Reset the cycle counter
    DWT->CYCCNT = 0;
    
    // Enable the cycle counter
    DWT->CTRL |= DWT_CTRL_CYCCNTENA_Msk;
}

/**
 * ----------------------------------------------------------------------------
 * Function: time_buffer_extract
 * ----------------------------------------------------------------------------
 * Purpose:
 *   Return captured timing records over SimpleSerial as a raw byte stream.
 *
 * Protocol:
 *   - Sends `time_captures * sizeof(uint32_t)` bytes as payload ('r' channel).
 *
 * Parameters:
 *   (cmd, scmd, len, buf) — SimpleSerial framing (unused).
 *
 * Returns:
 *   0x00 (SSv2 ack).
 * ----------------------------------------------------------------------------
 */

static uint8_t time_buffer_extract(uint8_t cmd, uint8_t scmd, uint8_t len, uint8_t *buf) {
    //returns up to 255 bytes from big_buffer, starting at big_buf_rd

    UNUSED(cmd); UNUSED(scmd); (void)len; UNUSED(buf);

    simpleserial_put('r', time_captures * sizeof(uint32_t), (uint8_t *)time_buffer);  
    return 0x00; 
}


/*------------------------- Prepare random Matrices and Fill Values -------------------------*/

/**
 * ----------------------------------------------------------------------------
 * Function: prepare_randmat_cmd
 * ----------------------------------------------------------------------------
 * Purpose:
 *   Populate all randomness tables used by masked χ:
 *     - randmat_and[x][y][i][j] : ISW pairwise masks for AND (i<j symmetric).
 *     - randmat_in / randmat_mid / randmat_out : refresh masks for operands
 *       before AND, after AND, and after final store, respectively.
 *
 * Randomness:
 *   - Uses true randomness via fill_pattern(..., pattern=0).
 *   - Diagonals in randmat_and are zero by ISW construction.
 *
 * Side-channel Notes:
 *   - All randomness is fresh per invocation to prevent reuse.
 *   - No secret input is processed here.
 *
 * Returns:
 *   0x00 (SSv2 ack).
 * ----------------------------------------------------------------------------
 */

static uint8_t prepare_randmat_cmd(uint8_t cmd, uint8_t scmd,
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
                        uint64_t temp[1];
                        fill_pattern(temp, 1, 0);  // pattern 0 = true random
                        randmat_and[x][y][i][j] = temp[0];
                        randmat_and[x][y][j][i] = temp[0];
                    }
                }
            }
        }
    }

    // --- Fill pre-AND refresh randomness (r_in) ---
    for (int y = 0; y < 5; y++) {
        for (int x = 0; x < 5; x++) {
            for (int i = 0; i < MASKING_N - 1; i++) {
                uint64_t temp[1];
                fill_pattern(temp, 1, 0);
                randmat_in[x][y][i] = temp[0];
            }
        }
    }

    // --- Fill post-AND refresh randomness (r_mid) ---
    for (int y = 0; y < 5; y++) {
        for (int x = 0; x < 5; x++) {
            for (int i = 0; i < MASKING_N - 1; i++) {
                uint64_t temp[1];
                fill_pattern(temp, 1, 0);
                randmat_mid[x][y][i] = temp[0];
            }
        }
    }

    // --- Fill post-store refresh randomness (r_out) ---
    for (int y = 0; y < 5; y++) {
        for (int x = 0; x < 5; x++) {
            for (int i = 0; i < MASKING_N - 1; i++) {
                uint64_t temp[1];
                fill_pattern(temp, 1, 0);
                randmat_out[x][y][i] = temp[0];
            }
        }
    }

    return 0x00;
}

/**
 * ----------------------------------------------------------------------------
 * Function: prepare_state_cmd
 * ----------------------------------------------------------------------------
 * Purpose:
 *   Initialize the global masked Keccak state with a chosen public pattern.
 *   Each 64-bit lane is split into MASKING_N Boolean shares with fresh
 *   randomness so that XOR of shares == lane_val.
 *
 * Inputs (buf[0]):
 *   0: random lane value (fresh per lane)
 *   1: 0x0000...0000
 *   2: 0xFFFF...FFFF
 *   3: 0xAAAAAAAAAAAAAAAA
 *   4: 0xA5A5A5A5A5A5A5A5
 *   5: 0x123456789ABCDEF0
 *  else: 0x0000...0000
 *
 * Side-channel Notes:
 *   - True randomness consumed for MASKING_N-1 shares per lane.
 *   - Resulting masking is first-order secure for linear ops.
 *
 * Returns:
 *   0x00 (SSv2 ack).
 * ----------------------------------------------------------------------------
 */

static uint8_t prepare_state_cmd(uint8_t cmd, uint8_t scmd,
                                 uint8_t len, uint8_t *buf)
{
    UNUSED(cmd); UNUSED(scmd);

    int pattern = buf[0];
    uint64_t lane_val = 0;
    chi_jitter_seed = get_rand64(); 

    switch (pattern) {
        case 0: { // random lane value
            lane_val = get_rand64();
        } break;
        case 1:  lane_val = 0x0ULL; break;
        case 2:  lane_val = ~0ULL; break;
        case 3:  lane_val = 0xAAAAAAAAAAAAAAAAULL; break;
        case 4:  lane_val = 0xA5A5A5A5A5A5A5A5ULL; break;
        case 5:  lane_val = 0x123456789ABCDEF0ULL; break;
        default: lane_val = 0ULL; break;
    }

    for (int y = 0; y < 5; y++) {
        for (int x = 0; x < 5; x++) {

            uint64_t acc = 0;

            // First MASKING_N-1 shares = random
            for (int s = 0; s < MASKING_N - 1; s++) {
                uint64_t r = get_rand64();
                global_state[x][y].share[s] = r;
                acc = xor64_safe(acc, r);

                for (int j = 0; j < 8; j++)
                    //random_data_buffer[random_data_buffer_length++] =
                        (uint8_t)(r >> (8 * j));
            }

            // Last share ensures XOR = lane_val
            uint64_t last = xor64_safe(lane_val, acc);
            global_state[x][y].share[MASKING_N - 1] = last;

            for (int j = 0; j < 8; j++)
                //random_data_buffer[random_data_buffer_length++] =
                    (uint8_t)(last >> (8 * j));
        }
    }

    return 0x00;
}


/**
 * ----------------------------------------------------------------------------
 * Function: prepare_state_real
 * ----------------------------------------------------------------------------
 * Purpose:
 *   Load up to 200 bytes of *real* input into the global masked Keccak state.
 *   Lanes are filled sequentially (8 bytes each), then Boolean-masked with
 *   fresh randomness (XOR shares).
 *
 * Parameters:
 *   buf : pointer to input data (SimpleSerial payload).
 *   len : number of bytes provided in buf (<= 200).
 *
 * Security / Side-channel Notes:
 *   - Fresh randomness for MASKING_N-1 shares per lane.
 *   - No permutation executed here; only state initialization.
 *
 * Returns:
 *   0x00 (SSv2 ack).
 * ----------------------------------------------------------------------------
 */

__attribute__((optimize("O0")))
__attribute__((noinline))
__attribute__((used))
static uint8_t prepare_state_real(uint8_t cmd, uint8_t scmd,
                                  uint8_t len, uint8_t *buf)
{
    UNUSED(cmd); UNUSED(scmd);

    // zero entire state
    for (int y = 0; y < 5; y++)
        for (int x = 0; x < 5; x++)
            for (int s = 0; s < MASKING_N; s++)
                global_state[x][y].share[s] = 0;

    // copy up to 200 bytes of real data (Keccak input)
    int total = len; // len from simpleserial frame
    uint8_t *p = buf;
    volatile uint8_t t1 = buf[0];
   volatile  uint8_t t2 = buf[1];
    volatile uint8_t t3 = buf[2];
    volatile uint8_t t4 = buf[3];
    int idx = 0;

    // fill lanes sequentially, 8 bytes per lane
    for (int y = 0; y < 5 && idx < total; y++) {
        for (int x = 0; x < 5 && idx < total; x++) {
            uint64_t lane_val = 0;
            int bytes_this_lane = (total - idx >= 8) ? 8 : (total - idx);
            memcpy(&lane_val, &p[idx], bytes_this_lane);
            idx += bytes_this_lane;

            // mask the lane_val as before
            uint64_t acc = 0;
            for (int s = 0; s < MASKING_N - 1; s++) {
                uint64_t r = get_rand64();
                global_state[x][y].share[s] = r;
                acc = xor64_safe(acc, r); 
            }
            global_state[x][y].share[MASKING_N - 1] = xor64_safe(lane_val, acc);
        }
    }

    return 0x00; // SSv2 ack
}

/**
 * ----------------------------------------------------------------------------
 * Function: prepare_iota_rands_cmd
 * ----------------------------------------------------------------------------
 * Purpose:
 *   Select the upcoming Iota round index and prefill the (N-1) random shares
 *   used to mask the round constant before XORing into lane (0,0).
 *
 * Inputs:
 *   buf[0] : round index in [0..23].
 *
 * Returns:
 *   0x00 on success; SS_ERR_LEN on malformed input.
 *
 * Side-channel Notes:
 *   - For MASKING_N > 1, RC is never applied as a single value; it is split
 *     into shares using the provided randomness to avoid leakage.
 * ----------------------------------------------------------------------------
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
 * ----------------------------------------------------------------------------
 * Function: fill_pattern
 * ----------------------------------------------------------------------------
 * Purpose:
 *   Fill an array of 64-bit words with either true randomness (pattern=0) or
 *   deterministic public patterns for test/calibration.
 *
 * Patterns:
 *   0: true random (from get_rand64)
 *   1: all zeros
 *   2: all ones
 *   3: 0xAAAAAAAAAAAAAAAA
 *   4: 0xA5A5A5A5A5A5A5A5
 *   5: 0x123456789ABCDEF0
 *
 * Side-channel Notes:
 *   - Always exercises RNG even when overwriting (balances power draw).
 * ----------------------------------------------------------------------------
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
 * ----------------------------------------------------------------------------
 * Function: clear_global_state
 * ----------------------------------------------------------------------------
 * Purpose:
 *   Zeroize the entire masked Keccak state (all shares of all 25 lanes).
 *
 * Side-channel Notes:
 *   - Clearing is constant-structure; no secret-dependent control flow.
 *
 * Returns:
 *   0x00 (SSv2 ack).
 * ----------------------------------------------------------------------------
 */

static uint8_t clear_global_state(uint8_t cmd ,uint8_t scmd,uint8_t len,uint8_t *buf)
{
    UNUSED(cmd); UNUSED(scmd); UNUSED(len); UNUSED(buf);

    for (int y = 0; y < 5; y++)
        for (int x = 0; x < 5; x++)
            for (int j = 0; j < MASKING_N; j++)
                global_state[x][y].share[j] = 0ULL;
    __DSB(); __ISB();
     return 0x00;
}

/**
 * ----------------------------------------------------------------------------
 * Function: clear_internal
 * ----------------------------------------------------------------------------
 * Purpose:
 *   Convenience wrapper to zeroize global state without SS framing.
 * ----------------------------------------------------------------------------
 */

static void clear_internal()
{
    clear_global_state(0, 0, 0, NULL);
}





// -----------------------------------------------------------------------
__attribute__((optimize("O0")))
static uint8_t nop_trace(uint8_t c,uint8_t s,uint8_t l,uint8_t *b){
    UNUSED(c); UNUSED(s); UNUSED(l); UNUSED(b);
    TEST_PROLOGUE();
     for (volatile uint32_t i = 0; i < 2400; i++) {
        __NOP();
    }
    TEST_EPILOGUE();
    return 0x00;
}


/**
 * ----------------------------------------------------------------------------
 * Function: rotl1_u64
 * ----------------------------------------------------------------------------
 * Purpose:
 *   Rotate a 64-bit value left by 1 bit (ROTL1).
 *   Used throughout Keccak (Θ and ρ steps) where a single-bit rotation
 *   implements the inter-column diffusion.
 *
 * Definition:
 *   ROTL1(x) = ((x << 1) | (x >> 63))  (mod 64)
 *
 * Security / Side-channel Notes:
 *   - Pure bitwise operation, no branching or data-dependent timing.
 *   - Safe under masking since rotation is a linear permutation of bits.
 *
 * Parameters:
 *   x  —  64-bit input word.
 *
 * Returns:
 *   64-bit value with bits rotated left by one position.
 * ----------------------------------------------------------------------------
 */

static inline uint64_t rotl1_u64(uint64_t x) {
    return (x << 1) | (x >> 63);
}


/**
 * build_C_twoacc:
 *   Balanced XOR tree for C[x] parity: ((a0 ⊕ a1) ⊕ a2) ⊕ (a3 ⊕ a4).
 *   Reduces single-accumulator toggling skew; ephemeral temps only.
 */

__attribute__((always_inline))
static inline void build_C_twoacc(volatile masked_uint64_t C[5],
                                  volatile masked_uint64_t st[5][5],
                                  int x, int i)
{
    // all temps die within this call
    uint64_t a0 = st[x][0].share[i];
    uint64_t a1 = st[x][1].share[i];
    uint64_t a2 = st[x][2].share[i];
    uint64_t a3 = st[x][3].share[i];
    uint64_t a4 = st[x][4].share[i];

    uint64_t accL = xor64_safe(xor64_safe(a0, a1), a2);   // lives briefly
    uint64_t accR = xor64_safe(a3, a4);        // lives briefly
    C[x].share[i] = xor64_safe(accL,  accR);

    // kill any temptation to carry regs across calls
    asm volatile ("" ::: "memory");
}

/**
 * masked_refresh:
 *   ISW refresh for a masked value v with (MASKING_N-1) random words r[i]:
 *     v[i]   ^= r[i] for i=0..N-2
 *     v[N-1] ^= r[i] for i=0..N-2
 *   Preserves XOR of shares; re-randomizes intermediate distribution.
 */

static inline void masked_refresh(volatile masked_uint64_t *v,
                                  const volatile uint64_t r_vec[MASKING_N - 1]) {
    // ISW Refresh: v[i] ^= r[i], v[N-1] ^= r[i] for i = 0..N-2
    for (size_t i = 0; i < MASKING_N - 1; i++) {
        uint64_t ri = r_vec[i];
        uint64_t vi = v->share[i];
        uint64_t vn1 = v->share[MASKING_N - 1];

        v->share[i]           = xor64_safe(vi, ri);
        v->share[MASKING_N-1] = xor64_safe(vn1, ri);
    }
}



void masked_xor(volatile masked_uint64_t *out,
                const volatile  masked_uint64_t *a,
                const volatile masked_uint64_t *b) {
    for (size_t i = 0; i < MASKING_N; i++) {
        out->share[i] = xor64_safe(a->share[i], b->share[i]);
    }
}

/**
 * masked_and:
 *   ISW-secure AND for Boolean-masked values. Diagonal terms (a_i & b_i)
 *   plus cross-terms corrected by shared randomness r[i][j], i<j.
 *   Implemented per MASKING_N specialization for clarity and efficiency.
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
 * masked_not:
 *   Masked bitwise NOT using a fixed public split K_NOT_SPLIT[share]:
 *     ~x  ≡  x ⊕ (all-ones). Here (all-ones) is represented as XOR of
 *     share-wise constants that sum to 0xFFFF...FFFF. Preserves masking.
 */
static inline void masked_not(volatile masked_uint64_t *dst, const volatile masked_uint64_t *src) {
    for (size_t i = 0; i < MASKING_N; ++i)
        dst->share[i] = xor64_safe(src->share[i],  K_NOT_SPLIT[i]);
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
void masked_theta(volatile masked_uint64_t state[5][5], uint32_t *time_records)
{
    volatile masked_uint64_t C[5], D[5];
    uint64_t tmp[5][MASKING_N];
    uint32_t start, end;

    // === 1. Compute C[x] ===
    start = DWT->CYCCNT;
    for (int x = 0; x < 5; x++)
        for (int i = 0; i < MASKING_N; i++)
            build_C_twoacc(C, state, x, i);
    end = DWT->CYCCNT;
    time_records[0] = end - start;   // e.g., ~1935 cycles - please note these may not be correct at release.

    // === 2. Compute D[x] ===
    start = DWT->CYCCNT;
    for (int x = 0; x < 5; x++) {
        int xm1 = (x + 4) % 5;
        int xp1 = (x + 1) % 5;

        for (int i = 0; i < MASKING_N; i++) {
            // temporary random mask for ISW safety (cancels algebraically)
            uint64_t m = get_rand64();

            // emulate ROTL1(C[xp1]) with masking linearity preserved
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
    time_records[1] = end - start;   // e.g., ~5105 cycles

    // === 3. Apply D[x] to all lanes ===
    start = DWT->CYCCNT;
    for (int x = 0; x < 5; x++)
        for (int i = 0; i < MASKING_N; i++) {
            uint64_t dx = D[x].share[i];
            for (int y = 0; y < 5; y++)
                state[x][y].share[i] = xor64_safe(state[x][y].share[i], dx);
        }
    end = DWT->CYCCNT;
    time_records[2] = end - start;   // e.g., ~3579 cycles
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
void masked_rho(volatile masked_uint64_t state[5][5], uint32_t *time_records) {
    uint32_t start, end;

    start = DWT->CYCCNT;
    for (int x = 0; x < 5; x++) {
        for (int y = 0; y < 5; y++) {
            uint8_t r = keccak_rho_offsets[x][y];
            for (int i = 0; i < MASKING_N; i++) {
                uint64_t rand_mask = get_rand64();

                // Apply random blinding before rotate
                uint64_t tmp = xor64_safe(state[x][y].share[i], rand_mask);

                // Rotate both the blinded value and the mask
                uint64_t tmp_rot = ROL64(tmp, r);
                uint64_t mask_rot = ROL64(rand_mask, r);

                // Unblind
                state[x][y].share[i] = xor64_safe(tmp_rot, mask_rot);
            }
        }
    }
    end = DWT->CYCCNT;
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

__attribute__((optimize("O0"))) //RUN CAPTURE 
void masked_pi(masked_uint64_t state[5][5]) {
    masked_uint64_t tmp[5][5];

    // Copy input to tmp
    for (int y = 0; y < 5; ++y)
        for (int x = 0; x < 5; ++x)
            for (int s = 0; s < MASKING_N; ++s)
                tmp[x][y].share[s] = state[x][y].share[s];

    // Pi: B[y, (2x + 3y) % 5] = A[x, y]
    for (int y = 0; y < 5; ++y) {
        for (int x = 0; x < 5; ++x) {
            int nx = y;
            int ny = (2 * x + 3 * y) % 5;
            for (int s = 0; s < MASKING_N; ++s)
                state[nx][ny].share[s] = tmp[x][y].share[s];
        }
    }
}

/*
void masked_pi(volatile masked_uint64_t state[5][5], uint32_t *time_records) {
    masked_uint64_t tmp[5][5];
    uint32_t start,end;

    start = DWT->CYCCNT;        
    for (int x = 0; x < 5; ++x)
        for (int y = 0; y < 5; ++y)
            for (int i = 0; i < MASKING_N; ++i) {
                uint64_t r = get_rand64();                    // new random mask
                tmp[x][y].share[i] = xor64_safe(state[x][y].share[i], r); // XOR once
                tmp[x][y].share[i] = xor64_safe(tmp[x][y].share[i], r); // XOR again to restore
            }

    end = DWT->CYCCNT;
    time_records[0] = end - start;
    start = DWT->CYCCNT;  
    for (int x = 0; x < 5; ++x)
        for (int y = 0; y < 5; ++y) {
            int new_x = y;
            int new_y = (2 * x + 3 * y) % 5;
            state[new_x][new_y] = tmp[x][y];
        }

    end = DWT->CYCCNT;
    time_records[1] = end - start;
}
*/


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

void masked_iota(masked_uint64_t state[5][5],
                 uint8_t round_idx,
                 const uint64_t *randoms,
                 uint32_t *time_records)
{

    #if (MASKING_N == 1)
        state[0][0].share[0] = xor64_safe(state[0][0].share[0], RC[round_idx]);
    #else
        // XOR only one share with a masked constant.
        uint64_t rc = RC[round_idx];
        uint64_t acc = 0;
        uint64_t masked_rc[MASKING_N];

        // Use provided randomness to mask the RC
        for (int i = 0; i < MASKING_N - 1; i++) {
            masked_rc[i] = randoms[i];
            acc = xor64_safe(acc,  masked_rc[i]);
        }
        masked_rc[MASKING_N - 1] = xor64_safe(rc, acc);

        // XOR each share of lane (0,0) with its masked constant
        for (int i = 0; i < MASKING_N; i++)
            state[0][0].share[i] = xor64_safe(state[0][0].share[i],  masked_rc[i]);
    #endif

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
    masked_uint64_t out[5][5],
    const masked_uint64_t in[5][5],
    const volatile uint64_t r_and[5][5][MASKING_N][MASKING_N],
    const volatile uint64_t r_in [5][5][MASKING_N - 1],
    const volatile uint64_t r_mid[5][5][MASKING_N - 1],
    const volatile uint64_t r_out[5][5][MASKING_N - 1])
{
    uint32_t t0 = DWT->CYCCNT;

    for (int y = 0; y < 5; y++) {
        for (int x = 0; x < 5; x++) {

            // Aliases to the three input lanes for this (x,y)
            const masked_uint64_t *a_ptr = &in[x][y];

            // Make local copies for operands we will mutate/refresh
            masked_uint64_t b_loc = in[(x + 1) % 5][y];
            masked_uint64_t c_loc = in[(x + 2) % 5][y];

            // === Refresh BOTH operands BEFORE any non-linear op ===
            // Use different randomness slices to avoid correlation
            masked_refresh(&b_loc, r_in[x][y]);
            masked_refresh(&c_loc, r_in[(x + 1) % 5][y]);
            __asm__ __volatile__("" ::: "memory");

            // Optional micro-jitter (mix in a fresh seed so it's not constant)
            uint64_t delay_seed;
            fill_pattern(&delay_seed, 1, 0);    // true random
            chi_jitter_seed = xor64_safe(chi_jitter_seed, delay_seed);
            chi_jitter_seed = chi_jitter_seed * 6364136223846793005ULL + 1ULL;
            int delay = (int)((chi_jitter_seed >> 8) & 0x7ULL);
            for (volatile int k = 0; k < delay; k++) {
                __asm__ __volatile__("nop");
            }

            // Temps must be real objects (not pointers into const input)
            __attribute__((aligned(16))) volatile masked_uint64_t t1, t2;

            // ~b (masked NOT using public split) on the REFRESHED b_loc
            masked_not(&t1, &b_loc);
            __asm__ __volatile__("" ::: "memory");

            // ISW AND with refreshed operands
            masked_and(&t2, &t1, &c_loc, r_and[x][y]);

            // Post-AND refresh + write out
            masked_refresh(&t2, r_mid[x][y]);
            __asm__ __volatile__("" ::: "memory");

            masked_xor(&out[x][y], a_ptr, &t2);
            masked_refresh(&out[x][y], r_out[x][y]);
            __asm__ __volatile__("" ::: "memory");

            // Zeroize all sensitive temporaries (without HD leakage)
            for (size_t i = 0; i < MASKING_N; i++) {
                t1.share[i] = 0;
                t2.share[i] = 0;
                b_loc.share[i] = 0;
                c_loc.share[i] = 0;
            }
            __asm__ __volatile__("" ::: "memory");
        }
    }

}

/*------------------------- Keccak Test Methods-------------------------*/

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
    masked_theta(global_state, time_records);
   end = DWT->CYCCNT;
    TEST_EPILOGUE();

    time = end - start;
    clear_internal();
    return 0x00;
}



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
    masked_rho(global_state, time_records);
   end = DWT->CYCCNT;
   
    TEST_EPILOGUE();
    time = end - start;

    clear_internal();
    return 0x00;
}


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
    masked_pi(global_state);   
    end = DWT->CYCCNT;
    
    TEST_EPILOGUE();

    time = end - start;
    clear_internal();
    return 0x00;
}


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
        masked_iota(global_state, iota_round_idx, iota_rands, time_records);
    #else
        masked_iota(global_state, RC[iota_round_idx], NULL);
    #endif
    end = DWT->CYCCNT;          
    TEST_EPILOGUE();
    time = end - start;
    clear_internal();

    return 0x00;  
}


__attribute__((optimize("O0")))
static uint8_t chi_test_no_save(uint8_t cmd, uint8_t scmd,
                                uint8_t len, uint8_t *buf)
{
    UNUSED(cmd); UNUSED(scmd); UNUSED(buf);

    masked_uint64_t out[5][5];
    uint32_t time_records[2] = {0};
    uint32_t start, end, time;

    __DSB();
    __ISB();
    TEST_PROLOGUE();
    NOP_BLOCK(128);
    start = DWT->CYCCNT;

    masked_chi(out, global_state,
               randmat_and,
               randmat_in,
               randmat_mid,
               randmat_out);

    end = DWT->CYCCNT;
    TEST_EPILOGUE();

    clear_internal();
    time = end - start;
    return 0x00;
}




// -----------------------------------------------------------------------
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
