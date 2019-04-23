This repository includes a very simple script to profile EnzoArray.

EnzoArray is the current implementation.
EnzoArrayN has a preliminary modification to allow for negative indexing.

The benchmarking compares the amount of time to access all elements in an
(N,N,N) for various array types. The benchmarked cases are:

- traditional pointer
- boost::multi_array using the traditional bracket notation [][][]
- boost::multi_array using the operator() method (where a container is passed)
- EnzoArray
- EnzoArrayN

This benchmark requires the cpuset utility and root access. The instructions
for running the benchmark are as follows:

1. modify N_REPS (number trials for each benchmark case and value of N) and RESERVED_THREAD (the cpu id to use for benchmarking)
2. Call "make" at the terminal.
3. Execute setup.sh as root (Note: the script turns off address space randomization and sets the scaling governor to performance - I don't think it successfully restores these to their initial values)
4. Run plot.py to plot the timings


The current plot saved to plot.pdf uses the minimum times to run each benchmark case take from 100 repetitions for each N. The top panel shows how the exectution of all cases scales with N, while the bottom panel shows the ratio of each case to the timings for the traditional pointer case
