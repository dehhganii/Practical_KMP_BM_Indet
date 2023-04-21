#!/bin/bash


./bm_indet $1 demofile_text demofile_pattern
./kmp_indet $1 demofile_text demofile_pattern
./brute_force $1 demofile_text demofile_pattern

