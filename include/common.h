/*
 * common.h — shared constants, macros and small helpers for DNA-GARS.
 */
#ifndef DNA_GARS_COMMON_H
#define DNA_GARS_COMMON_H

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

/* ---- Alphabet ---- */
#define DNA_ALPHABET "ACGT"

/* ---- Rabin-Karp ---- */
#define RK_BASE 4ULL            /* A=0, C=1, G=2, T=3 */
#define RK_MOD  1000000007ULL   /* large prime; set to e.g. 101 to see spurious hits */

/* ---- Timing ---- */
#define MIN_TIMING_MS   2.0     /* a timed batch of runs lasts at least this long */
#define MAX_TIMING_RUNS 100000  /* ...or contains this many runs, whichever comes first */
#define TIMING_BATCHES  5       /* batches timed; the fastest per-run time is reported */

/* ---- Alignment ---- */
#define DEFAULT_MATCH     1
#define DEFAULT_MISMATCH -1
#define DEFAULT_GAP      -2
#define MATRIX_PRINT_MAX 15     /* print the DP matrix only when both inputs are this short */
#define MAX_ALIGN_CELLS  100000000LL

/* ---- Reconstruction (Shortest Common Superstring) ---- */
#define BRUTE_FORCE_MAX_K 9     /* k! permutations */
#define DP_BITMASK_MAX_K  18    /* k^2 * 2^k states */

#define MIN(a, b) ((a) < (b) ? (a) : (b))
#define MAX(a, b) ((a) > (b) ? (a) : (b))

/* ---- Allocation wrappers: abort with a message instead of returning NULL ---- */
static inline void *xmalloc(size_t size) {
    void *p = malloc(size ? size : 1);
    if (!p) {
        fprintf(stderr, "Out of memory (requested %zu bytes)\n", size);
        exit(EXIT_FAILURE);
    }
    return p;
}

static inline void *xcalloc(size_t count, size_t size) {
    void *p = calloc(count ? count : 1, size ? size : 1);
    if (!p) {
        fprintf(stderr, "Out of memory (requested %zu x %zu bytes)\n", count, size);
        exit(EXIT_FAILURE);
    }
    return p;
}

static inline void *xrealloc(void *ptr, size_t size) {
    void *p = realloc(ptr, size ? size : 1);
    if (!p) {
        fprintf(stderr, "Out of memory (requested %zu bytes)\n", size);
        exit(EXIT_FAILURE);
    }
    return p;
}

static inline char *xstrdup(const char *s) {
    size_t n = strlen(s) + 1;
    char *d = xmalloc(n);
    memcpy(d, s, n);
    return d;
}

/* Map a base to 0..3 (used by Rabin-Karp hashing). */
static inline int dna_code(char c) {
    switch (c) {
    case 'A': return 0;
    case 'C': return 1;
    case 'G': return 2;
    case 'T': return 3;
    default:  return 0;
    }
}

#endif /* DNA_GARS_COMMON_H */
