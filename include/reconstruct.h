/*
 * reconstruct.h — Module 3: genome reconstruction as the
 * Shortest Common Superstring (SCS) problem (NP-hard).
 *
 * Pipeline:
 *   1. fragments_remove_contained()  drop duplicates / substrings of other fragments
 *   2. overlap_matrix()              ov[i][j] = longest suffix of i that is a prefix of j
 *   3. scs_*()                       choose an order of fragments maximising total overlap
 *   4. superstring_from_order()      merge fragments in that order
 *
 * Maximising total overlap along a path == minimising superstring length.
 */
#ifndef DNA_GARS_RECONSTRUCT_H
#define DNA_GARS_RECONSTRUCT_H

#include "common.h"

typedef struct {
    char *superstring;
    int length;
    int *order;             /* fragment indices in merge order */
    int k;                  /* number of fragments */
    long long operations;   /* overlap lookups / edges examined */
    double time_ms;         /* filled by recon_run() */
} ReconResult;

typedef ReconResult (*ReconFn)(char **frags, int k, int **ov);

/* ---- overlap.c ---- */

/* Longest suffix of a that equals a prefix of b, via the KMP prefix function
 * on b + '#' + a. Adds comparisons to *ops (may be NULL). */
int overlap_length(const char *a, const char *b, long long *ops);

/* k x k overlap matrix (diagonal = 0). Free with overlap_matrix_free(). */
int **overlap_matrix(char **frags, int k, long long *ops);
void overlap_matrix_free(int **ov, int k);

/* Remove fragments that are duplicates of, or contained in, another fragment
 * (frees them and compacts the array). Returns the new count. */
int fragments_remove_contained(char **frags, int k);

/* Merge fragments in the given order using the overlap matrix. Caller frees. */
char *superstring_from_order(char **frags, int k, int **ov, const int *order);

/* Fill a ReconResult (superstring, length, copy of order) from an order. */
ReconResult recon_result_from_order(char **frags, int k, int **ov, const int *order,
                                    long long operations);

void recon_result_free(ReconResult *r);

/* ---- Algorithms ---- */

/* Exact: try all k! orders.                 O(k! * k) time, O(k) space */
ReconResult scs_brute_force(char **frags, int k, int **ov);

/* Exact: Held-Karp style bitmask DP.         O(k^2 * 2^k) time, O(k * 2^k) space */
ReconResult scs_dp_bitmask(char **frags, int k, int **ov);

/* Approximate: take overlap edges in decreasing order, keeping a single path.
 *                                            O(k^2 log k) time after overlaps */
ReconResult scs_greedy(char **frags, int k, int **ov);

#endif /* DNA_GARS_RECONSTRUCT_H */
