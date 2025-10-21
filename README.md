Each of these methods has been designed for side-channel security and has been evaluated at first order masking with TVLA order 1-4 for the linnear operations 
theta, rho, pi and iota and at masking orders 1-3 with TVLA orders 1-4 for the only non-linnear method chi. 

This is not a claim of security, simply that it is expected not to leak information with a T value greater than 4.5.

1st order evaluation 50k vs 50k fixed vs random

2nd order evaluation 200k vs 200k fixed vs random

3rd order evaluation 500k vs 500k fixed vs random


This is designed for use with the ChipWhisperer and assumes the ChipWhisperer HAL. Please download the version 6.0.0 of the ChipWhisperer software and place this folder inside /chipwhisperer/firmware/mcu and the contained headers and makefile will work.

# Masked Keccak-f[1600] Evaluation Firmware

**Author:** Adam Beattie  
**Platform:** STM32F4 (ChipWhisperer Husky target)  
**Date:** 2025  

---

## Overview
This firmware implements an **arbitrary-order Boolean-masked Keccak-f[1600] permutation** for cycle-accurate **side-channel evaluation**.  
It is designed as a minimal **capture target**, providing direct access to each round function component (Θ, ρ, π, χ, ι) via SimpleSerial commands.

Only **one capture file and Makefile** are included in this stage — the purpose is to verify the **core masked Keccak permutation** before integrating the full **SHA-3 sponge construction**.

---

## Current Capabilities
- ✔ Arbitrary-order Boolean masking (`MASKING_ORDER + 1` shares)
- ✔ Full Keccak step coverage: Θ, ρ, π, χ, ι
- ✔ Cycle-accurate timing via DWT->CYCCNT
- ✔ Controlled trigger boundaries for reproducible captures
- ✔ Randomness preparation commands for χ and ι
- ✔ Minimal dependencies (HAL + SimpleSerial)

---

## Commands (SimpleSerial)
| Command | Function               |
|----------|------------------------|
| `0x30`   | Run masked Θ (Theta)   |
| `0x31`   | Run masked ρ (Rho)     |
| `0x32`   | Run masked π (Pi)      |
| `0x33`   | Run masked ι (Iota)    |
| `0x34`   | Run masked χ (Chi)     |
| `0x20`   | Prepare test state     |
| `0x21`   | Prepare Iota randomness |
| `0x22`   | Prepare Chi randomness  |
| `0x70`   | Fixed-length NOP trace (alignment test) |

---

## Build and Flash
```bash
make PLATFORM=CW308_STM32F4 CRYPTO_TARGET=MASKED_KECCAK
