/*
 * brute_force.c — exact SCS by trying every order of the fragments.
 *
 * For each of the k! permutations, the superstring length is
 *     sum of fragment lengths - sum of overlaps between consecutive fragments,
 * so the shortest superstring is the permutation with the largest total overlap.
 * Permutations are generated in lexicographic order (next_permutation).
 *
 * Time:  O(k! * k) after the overlap matrix
 * Space: O(k)
 * Practical limit: k <= BRUTE_FORCE_MAX_K (9! = 362,880 orders)
 */
#include "reconstruct.h"

/* Rearrange a[] into the next lexicographic permutation; 0 when done. */
static int next_permutation(int *a, int n) {
    int i = n - 2;
    while (i >= 0 && a[i] >= a[i + 1])
        i--;
    if (i < 0)
        return 0;
    int j = n - 1;
    while (a[j] <= a[i])
        j--;
    int t = a[i]; a[i] = a[j]; a[j] = t;
    for (int l = i + 1, r = n - 1; l < r; l++, r--) {
        t = a[l]; a[l] = a[r]; a[r] = t;
    }
    return 1;
}

ReconResult scs_brute_force(char **frags, int k, int **ov) {
    int *perm = xmalloc((size_t)(k ? k : 1) * sizeof(int));
    int *best = xmalloc((size_t)(k ? k : 1) * sizeof(int));
    long long ops = 0;
    int best_overlap = -1;

    for (int i = 0; i < k; i++)
        perm[i] = best[i] = i;

    if (k > 0) {
        do {
            int total = 0;
            for (int i = 0; i + 1 < k; i++)
                total += ov[perm[i]][perm[i + 1]];
            ops += k - 1;
            if (total > best_overlap) {
                best_overlap = total;
                memcpy(best, perm, (size_t)k * sizeof(int));
            }
        } while (next_permutation(perm, k));
    }

    ReconResult r = recon_result_from_order(frags, k, ov, best, ops);
    free(perm);
    free(best);
    return r;
}
