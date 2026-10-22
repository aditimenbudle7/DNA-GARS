/*
 * search.h — Module 1: exact pattern search in a genome.
 *
 * Every algorithm has the same signature and returns a SearchResult, so the
 * menu, tests and benchmarks can treat them uniformly (see SEARCH_ALGOS in
 * benchmark.h). Positions are 0-based.
 */
#ifndef DNA_GARS_SEARCH_H
#define DNA_GARS_SEARCH_H

#include "common.h"

typedef struct {
    int *positions;          /* match start positions, ascending, 0-based */
    int count;
    int capacity;
    long long comparisons;   /* character comparisons, including preprocessing */
    long long hash_checks;   /* Rabin-Karp only: window hashes compared */
    long long spurious_hits; /* Rabin-Karp only: hash matched but text did not */
    double time_ms;          /* time per run, filled by search_run() */
} SearchResult;

typedef SearchResult (*SearchFn)(const char *text, int n, const char *pat, int m);

static inline SearchResult search_result_new(void) {
    SearchResult r;
    memset(&r, 0, sizeof r);
    return r;
}

static inline void search_result_add(SearchResult *r, int pos) {
    if (r->count == r->capacity) {
        r->capacity = r->capacity ? r->capacity * 2 : 16;
        r->positions = xrealloc(r->positions, (size_t)r->capacity * sizeof(int));
    }
    r->positions[r->count++] = pos;
}

static inline void search_result_free(SearchResult *r) {
    free(r->positions);
    r->positions = NULL;
    r->count = r->capacity = 0;
}

/* Brute force: try every shift.                        O(n*m) time, O(1) space */
SearchResult naive_search(const char *text, int n, const char *pat, int m);

/* Knuth-Morris-Pratt with prefix (failure) function.   O(n+m) time, O(m) space */
SearchResult kmp_search(const char *text, int n, const char *pat, int m);

/* Rabin-Karp rolling hash, base 4.                     O(n+m) avg, O(n*m) worst */
SearchResult rabin_karp_search(const char *text, int n, const char *pat, int m);

/* Boyer-Moore, bad character + strong good suffix.     O(n/m) best, O(n*m) worst */
SearchResult boyer_moore_search(const char *text, int n, const char *pat, int m);

/* Z-algorithm on pat + '$' + text.                     O(n+m) time and space */
SearchResult z_search(const char *text, int n, const char *pat, int m);

/* KMP prefix function: lps[i] = length of the longest proper prefix of p[0..i]
 * that is also a suffix of it. Adds character comparisons to *comparisons
 * (may be NULL). Reused by Module 3 to compute fragment overlaps. */
void kmp_prefix_function(const char *p, int m, int *lps, long long *comparisons);

#endif /* DNA_GARS_SEARCH_H */
