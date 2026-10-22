/*
 * benchmark.h — algorithm registry, timed runners and CSV benchmarks.
 */
#ifndef DNA_GARS_BENCHMARK_H
#define DNA_GARS_BENCHMARK_H

#include "search.h"
#include "align.h"
#include "reconstruct.h"

/* ---- Pattern search registry (Naive first: it is the reference) ---- */

typedef struct {
    const char *name;
    SearchFn fn;
} SearchAlgo;

extern const SearchAlgo SEARCH_ALGOS[];
extern const int SEARCH_ALGO_COUNT;

/* 1 if both results report exactly the same positions. */
int search_results_equal(const SearchResult *a, const SearchResult *b);

/* ---- Timed runners ----
 * Each runs the algorithm once for its result, then times repeated batches of
 * runs (each batch >= MIN_TIMING_MS) and stores the fastest per-run time in
 * time_ms. This keeps tiny inputs from showing 0.000 ms and filters out noise. */

SearchResult search_run(SearchFn fn, const char *text, int n, const char *pat, int m);
AlignResult  align_run(AlignFn fn, const char *a, const char *b, ScoreScheme s);
ReconResult  recon_run(ReconFn fn, char **frags, int k, int **ov);

/* ---- Benchmarks: sweep input sizes, print a table, write CSV ----
 * Each returns 0 on success, -1 if the CSV could not be written. */

int bench_search(const char *csv_path);
int bench_align(const char *csv_path);
int bench_reconstruct(const char *csv_path);

/* Runs all three, writing search.csv, align.csv and reconstruct.csv into dir. */
int bench_all(const char *dir);

#endif /* DNA_GARS_BENCHMARK_H */
