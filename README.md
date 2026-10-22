# DNA-GARS
**DNA Genome Analysis and Reconstruction System**

A C implementation and empirical comparison of classic algorithms applied to DNA sequences, built for the Design and Analysis of Algorithms (DAA) course.

## Modules

| Module | Problem | Algorithms |
|---|---|---|
| 1. Pattern Search | Find all occurrences of a pattern in a genome | Naive, KMP, Rabin-Karp, Boyer-Moore, Z-Algorithm |
| 2. Sequence Alignment | Align two sequences and detect mutations | Needleman-Wunsch (global), Smith-Waterman (local) |
| 3. Genome Reconstruction | Shortest Common Superstring from fragments | Brute Force, Bitmask DP, Greedy (approximation) |

Every algorithm reports its result, number of comparisons/operations, and measured execution time. All counts and timings are measured by the program, not hard-coded.

## Features
- Input from keyboard, plain text, or FASTA files (validated to `A/C/G/T`)
- Random genome generator and fragment generator (for testing reconstruction)
- Cross-verification of all pattern-search algorithms against the naive baseline
- Alignment display with DP matrix (small inputs) and mutation report (substitutions, insertions, deletions, % identity)
- Benchmark mode: sweeps input sizes and exports CSV to `results/`
- Plot script for time/comparisons vs. input size

## Project Structure
```
dna-gars/
├── Makefile
├── include/          # headers
├── src/
│   ├── main.c        # menu-driven CLI
│   ├── search/       # Module 1
│   ├── align/        # Module 2
│   └── reconstruct/  # Module 3
├── tests/            # correctness tests
├── data/             # sample genomes and fragment sets
├── results/          # generated CSVs (not committed)
└── scripts/plot.py   # graphs from CSV
```

## Build & Run
Requires `gcc` and `make` (Linux/macOS, or MinGW/WSL on Windows).

```bash
make            # build
./dna-gars      # run (dna-gars.exe on Windows)
make test       # run tests
make clean      # remove build files
```

Plot benchmark results (optional, needs Python 3 + matplotlib):
```bash
python3 scripts/plot.py results/
```

## Conventions
- Positions are reported **0-based**.
- Time is the average over multiple runs, in milliseconds.

## Complexity Summary

| Algorithm | Time | Space |
|---|---|---|
| Naive | O(n·m) | O(1) |
| KMP | O(n + m) | O(m) |
| Rabin-Karp | O(n + m) avg, O(n·m) worst | O(1) |
| Boyer-Moore | O(n/m) best, O(n·m) worst | O(σ) |
| Z-Algorithm | O(n + m) | O(n + m) |
| Needleman-Wunsch | O(n·m) | O(n·m) |
| Smith-Waterman | O(n·m) | O(n·m) |
| SCS Brute Force | O(k!·L) | O(k·L) |
| SCS Bitmask DP | O(k²·2ᵏ) | O(k·2ᵏ) |
| SCS Greedy | O(k³·L) | O(k²) |

*n = genome length, m = pattern length, σ = alphabet size, k = fragments, L = fragment length*

## Author
Aditi — DAA Course Project