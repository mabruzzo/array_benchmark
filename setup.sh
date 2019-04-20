#!/bin/bash

# following the instructions from
# http://releases.llvm.org/7.0.1/docs/Benchmarking.html

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
# they are an SMT pain on my machine, so as long as I use 1 thread, SMT won't
# be a problem
cset shield -c 3,7 -k on

cset shield --exec -- ./test 64 15 result_64.txt
cset shield --exec -- ./test 128 15 result_128.txt
cset shield --exec -- ./test 256 15 result_256.txt

# cleanup
cset shield --reset

# reset the scaling_governors
index=0
for i in /sys/devices/system/cpu/cpu*/cpufreq/scaling_governor
do
    echo ${default_scaling[$index]} > ${i}
    cat ${i}
    index=$((index+1))
done

# re-enable address space randomization
echo ${default_random} > /proc/sys/kernel/randomize_va_space
