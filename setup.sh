#!/bin/bash

# loosely following the instructions from
# http://releases.llvm.org/7.0.1/docs/Benchmarking.html
# requires cpset

# Instructions:
# 0. Modify the following 2 variables (RESERVED_THREAD and N_REPS) 
# 1. Call "make" on the terminal
# 2. Execute this script (with root)
# 3. call python plot.py which generates plot.pdf

#Modify These 2 variables to suit your machine/needs

# indicate which thread(s) to reserve
#   - Example: cores number 3 and 7 they are an SMT pair on my machine, so if I
#              reserve both and only run 1 process, I don't think SMT
#              (hyperthreading) will be a problem
RESERVED_THREAD=3,7

# Number of repetitions to do per trial, per array
N_REPS=25


# disable address space randomization
default_random="$(cat /proc/sys/kernel/randomize_va_space)"
echo 0 > /proc/sys/kernel/randomize_va_space

# set the scaling governor to performance
index=0
declare -a default_scaling

for i in /sys/devices/system/cpu/cpu*/cpufreq/scaling_governor
do
    default_scaling[$index]="$(cat ${i})"
    echo performance > ${i}
    index=$((index+1))
done





# reserve threads 3 and 7

cset shield -c ${RESERVED_THREAD} -k on


cset shield --exec -- ./test 32 ${N_REPS} result_32.txt
cset shield --exec -- ./test 64 ${N_REPS} result_64.txt
cset shield --exec -- ./test 128 ${N_REPS} result_128.txt
cset shield --exec -- ./test 256 ${N_REPS} result_256.txt
cset shield --exec -- ./test 512 ${N_REPS} result_512.txt

# cleanup
cset shield --reset

# reset the scaling_governors back to the defaults (I don't think this works
# right)
index=0
for i in /sys/devices/system/cpu/cpu*/cpufreq/scaling_governor
do
    echo ${default_scaling[$index]} > ${i}
    cat ${i}
    index=$((index+1))
done

# re-enable address space randomization (I don't think this works quite right)
echo ${default_random} > /proc/sys/kernel/randomize_va_space
