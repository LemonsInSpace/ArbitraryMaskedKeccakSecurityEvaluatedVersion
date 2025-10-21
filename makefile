# Hey Emacs, this is a -*- makefile -*-
#----------------------------------------------------------------------------
# Makefile for ChipWhisperer SimpleSerial-Keccak Program
#----------------------------------------------------------------------------

#This repo does not contain the associated chipwhisperer HAL. Download and install ChipWhisperer 6.0.0 and place this folder inside 
#/chipwhisperer/firmware/mcu

# Target file name (without extension)
TARGET = simpleserial-base

SS_VER = SS_VER_2_1

# List C source files here
SRC += keccak_tvla.c \

CFLAGS  += -Og -mlong-calls -g
LDFLAGS += -mlong-calls


# List additional headers if required
# (Not strictly necessary unless custom build rules)
# HDRS += masked_keccak.h masked_gadgets.h masked_types.h sha_shake.h params.h

# -----------------------------------------------------------------------------
# Add simpleserial project to build
include ../simpleserial/Makefile.simpleserial

FIRMWAREPATH = ../.
include $(FIRMWAREPATH)/Makefile.inc