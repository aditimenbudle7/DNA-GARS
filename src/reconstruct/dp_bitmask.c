/*
 * dp_bitmask.c — exact SCS with bitmask dynamic programming (Held-Karp style).
 *
 * dp[mask][last] = maximum total overlap of a path that visits exactly the
 *                  fragments in `mask` and ends at fragment `last`.
 *
 *   dp[{i}][i] = 0
 *   dp[mask | {next}][next] = max over last in mask of dp[mask][last] + ov[last][next]
 *
 * Answer: max over last of dp[all][last]; the order is rebuilt from parent[].
 * Same idea as the DP for the Travelling Salesman Problem, without returning
 * to the start.
 *
 * Time:  O(k^2 * 2^k)
 * Space: O(k * 2^k)
 * Practical limit: k <= DP_BITMASK_MAX_K
 */
#include "reconstruct.h"

#define IDX(mask, i) ((size_t)(mask) * (size_t)k + (size_t)(i))

ReconResult scs_dp_bitmask(char **frags, int k, int **ov) {
    if (k == 0)
        return recon_result_from_order(frags, 0, ov, NULL, 0);

    size_t states = (size_t)1 << k;
    int full = (int)(states - 1);
    int *dp = xmalloc(states * (size_t)k * sizeof(int));
    int *parent = xmalloc(states * (size_t)k * sizeof(int));
    long long ops = 0;

    for (size_t t = 0; t < states * (size_t)k; t++) {
        dp[t] = -1;          /* unreachable */
        parent[t] = -1;
    }
    for (int i = 0; i < k; i++)
        dp[IDX(1 << i, i)] = 0;

    for (int mask = 1; mask <= full; mask++) {
        for (int last = 0; last < k; last++) {
            if (!(mask & (1 << last)))
                continue;
            int cur = dp[IDX(mask, last)];
            if (cur < 0)
                continue;
            for (int next = 0; next < k; next++) {
                if (mask & (1 << next))
                    continue;
                ops++;
                int nmask = mask | (1 << next);
                int val = cur + ov[last][next];
                if (val > dp[IDX(nmask, next)]) {
                    dp[IDX(nmask, next)] = val;
                    parent[IDX(nmask, next)] = last;
                }
            }
        }
    }

    /* Best end fragment, then walk parents backwards */
    int last = 0;
    for (int i = 1; i < k; i++)
        if (dp[IDX(full, i)] > dp[IDX(full, last)])
            last = i;

    int *order = xmalloc((size_t)k * sizeof(int));
    int mask = full;
    for (int pos = k - 1; pos >= 0; pos--) {
        order[pos] = last;
        int prev = parent[IDX(mask, last)];
        mask &= ~(1 << last);
        last = prev;
    }

    ReconResult r = recon_result_from_order(frags, k, ov, order, ops);
    free(order);
    free(dp);
    free(parent);
    return r;
}
