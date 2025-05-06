# Graph Biconnectivity Components Benchmarking

This project is part of the Introduction to CS1.305 Algorithm Engineering Course at IIIT-H, as a course project. Here, we implement and benchmarks several algorithms across various dataset sizes and implementation for finding biconnected components in graphs.

## Implemented Algorithms

1. **Slota-Madduri**
   - Sequential implementation: [slota_madduri/slota_madduri_seq.c](./slota_madduri/slota_madduri_seq.c)
   - Parallel implementation: [slota_madduri/slota_madduri_par.c](./slota_madduri/slota_madduri_par.c)

2. **Tarjan-Vishkin**
   - Sequential implementation: [tarjan-viskin/tarjan_vishkin_seq.c](./tarjan-viskin/tarjan_vishkin_seq.c)
   - Parallel implementation: [tarjan-viskin/tarjan_vishkin_par.c](./tarjan-viskin/tarjan_vishkin_par.c)

3. **Jen-Schmidt**
   - Sequential implementation: [jen-schmidt/jen_schmidt_seq.c](./jen-schmidt/jen_schmidt_seq.c)
   - Parallel implementation: [jen-schmidt/jen_schmidt_par.c](./jen-schmidt/jen_schmidt_par.c)

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

## Performance Results

### Sequential Performance

[TODO] Performance for the other algorithms will be measured and documented after running the benchmarks.

## Results Analysis

[TODO] Analysis of the performance results will be added here, including comparisons between sequential and parallel implementations, and insights into the scalability of each algorithm.

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
./slota_seq datasets/ custom.txt
OMP_NUM_THREADS=4 ./slota_par datasets/ custom.txt
```

## Datasets
| Dataset | Web URL | Nodes | Edges | Description |
|-------------|----------|-------|-------|-------------|
| [custom-built](./datasets/custom.txt)       | - | 50 | 128 | Custom built small set for testing |
| [congress-Twitter](./datasets/congress.txt)       | https://snap.stanford.edu/data/congress-twitter.html | 475 | 13,289 | Twitter interaction network for the US Congress |
| [email-Eu-Core](./datasets/email.txt)      | https://snap.stanford.edu/data/email-Eu-core.html | 1,005 | 25,571 | E-mail network |
| [ego-Facebook](./datasets/facebook.txt)       | https://snap.stanford.edu/data/ego-Facebook.html | 4,039 | 88,234 | Facebook social circles network |
| [ca-HepPh](./datasets/hep.txt)       | https://snap.stanford.edu/data/ca-HepPh.html | 12,008 | 118,521 | Collaboration network of Arxiv High Energy Physics |
| [musae-Facebook](./datasets/musae.txt)       | https://snap.stanford.edu/data/facebook-large-page-page-network.html | 22,469 | 171,002 | Facebook page-page network with page names |
| [musae-github](./datasets/github.txt)       | https://snap.stanford.edu/data/github-social.html | 37,700 | 289,003 | Social network of Github developers |
| [soc-Slashdot0922](./datasets/slashdot.txt)       | https://snap.stanford.edu/data/soc-Slashdot0902.html | 82,167 | 948,464 | Slashdot social network from February 2009 |