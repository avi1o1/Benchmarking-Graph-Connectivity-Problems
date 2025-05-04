https://snap.stanford.edu/data/twitch_gamers.html
https://snap.stanford.edu/data/ego-Facebook.html

# Graph Biconnectivity Components Benchmarking

This project implements and benchmarks several algorithms for finding biconnected components in graphs.

## Implemented Algorithms

1. **Slota-Madduri**
   - Sequential implementation: `slota_madduri/slota_madduri_seq.c`
   - Parallel implementation: `slota_madduri/slota_madduri_par.c`

2. **Tarjan-Vishkin**
   - Sequential implementation: `tarjan-viskin/tarjan_vishkin_seq.c` 
   - Parallel implementation: `tarjan-viskin/tarjan_vishkin_par.c`

3. **Jen-Schmidt**
   - Sequential implementation: `jen-schmidt/jen_schmidt_seq.c`
   - Parallel implementation: `jen-schmidt/jen_schmidt_par.c`

## Performance Results

### Sequential Performance
**Slota Madduri**
- Small: 0.00554 seconds ; 0.00562 seconds
- Medium: 1.37 seconds ; 1.34 seconds
- Large: 44444.38883 seconds ; 44463.96789 seconds

Performance for the other algorithms will be measured and documented after running the benchmarks.

## Datasets
The algorithms are tested on three datasets of different sizes:
- Small: `datasets/small.txt`
- Medium: `datasets/medium.txt` 
- Large: `datasets/large.txt`

Dataset sources:
- https://snap.stanford.edu/data/twitch_gamers.html
- https://snap.stanford.edu/data/ego-Facebook.html

## Running the Benchmarks

The project includes two makefiles to automate the process of compiling and running the algorithms.

### Sequential Algorithms

```bash
# Compile and run all sequential implementations on all datasets
make -f Makefile.seq

# Compile only
make -f Makefile.seq slota_madduri_seq jen_schmidt_seq tarjan_vishkin_seq

# Run on specific dataset
make -f Makefile.seq run_small_jen_schmidt_seq  # For small dataset
make -f Makefile.seq run_medium_tarjan_vishkin_seq  # For medium dataset

# Generate a summary report of all results
make -f Makefile.seq report

# Clean up binaries and results
make -f Makefile.seq clean
```

### Parallel Algorithms

```bash
# Compile and run all parallel implementations on all datasets with various thread counts
make -f Makefile.par

# Compile only
make -f Makefile.par slota_madduri_par jen_schmidt_par tarjan_vishkin_par

# Run with specific thread count on specific dataset
make -f Makefile.par run_4_small_jen_schmidt_par  # 4 threads, small dataset
make -f Makefile.par run_8_large_tarjan_vishkin_par  # 8 threads, large dataset

# Generate summary and scaling analysis reports
make -f Makefile.par report

# Clean up binaries and results
make -f Makefile.par clean
```

## Results Analysis

The sequential results are saved in the `results_seq` directory, while parallel results are saved in the `results_par` directory.

For parallel implementations, a scaling analysis report is generated that compares performance across different thread counts (2, 4, 8, and 16 threads).

## Compiling Individual Programs

If you want to compile and run specific implementations directly:

```bash
# Sequential implementations
gcc -O3 -Wall -o slota_seq slota_madduri/slota_madduri_seq.c
gcc -O3 -Wall -o jen_seq jen-schmidt/jen_schmidt_seq.c
gcc -O3 -Wall -o tarjan_seq tarjan-viskin/tarjan_vishkin_seq.c

# Parallel implementations
gcc -O3 -Wall -fopenmp -o slota_par slota_madduri/slota_madduri_par.c
gcc -O3 -Wall -fopenmp -o jen_par jen-schmidt/jen_schmidt_par.c
gcc -O3 -Wall -fopenmp -o tarjan_par tarjan-viskin/tarjan_vishkin_par.c
```

And run them with:

```bash
./slota_seq datasets/small.txt
OMP_NUM_THREADS=4 ./slota_par datasets/medium.txt
```