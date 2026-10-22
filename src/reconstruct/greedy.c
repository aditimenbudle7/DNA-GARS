/*
 * greedy.c — approximate SCS by greedy merging.
 *
 * Treat fragments as nodes and overlaps as weighted edges i -> j. Repeatedly
 * take the heaviest remaining edge (largest overlap) as long as
 *   - i has no outgoing edge yet,
 *   - j has no incoming edge yet,
 *   - it does not close a cycle (checked with union-find).
 * After k-1 edges the chosen edges form a single path: the merge order.
 * This is the classic "merge the pair with maximum overlap" greedy.
 *
 * Not always optimal; it is conjectured to be within 2x of the optimal length
 * (a factor of 3.5 is proven). The benchmark compares it to the exact DP.
 *
 * Time:  O(k^2 log k) after the overlap matrix (sorting the edges)
 * Space: O(k^2)
 */
#include "reconstruct.h"

typedef struct {
    int w, from, to;
} Edge;

static int edge_cmp(const void *pa, const void *pb) {
    const Edge *a = pa, *b = pb;
    if (a->w != b->w)
        return b->w - a->w;                 /* heavier first */
    if (a->from != b->from)
        return a->from - b->from;           /* deterministic tie-break */
    return a->to - b->to;
}

static int uf_find(int *p, int x) {
    while (p[x] != x) {
        p[x] = p[p[x]];                     /* path halving */
        x = p[x];
    }
    return x;
}

ReconResult scs_greedy(char **frags, int k, int **ov) {
    if (k == 0)
        return recon_result_from_order(frags, 0, ov, NULL, 0);

    int ne = k * (k - 1);
    Edge *edges = xmalloc((size_t)(ne ? ne : 1) * sizeof(Edge));
    int e = 0;
    for (int i = 0; i < k; i++)
        for (int j = 0; j < k; j++)
            if (i != j) {
                Edge ed = { ov[i][j], i, j };
                edges[e++] = ed;
            }
    qsort(edges, (size_t)ne, sizeof(Edge), edge_cmp);

    int *next = xmalloc((size_t)k * sizeof(int));
    int *has_in = xcalloc((size_t)k, sizeof(int));
    int *uf = xmalloc((size_t)k * sizeof(int));
    for (int i = 0; i < k; i++) {
        next[i] = -1;
        uf[i] = i;
    }

    long long ops = 0;
    int taken = 0;
    for (int t = 0; t < ne && taken < k - 1; t++) {
        ops++;
        int a = edges[t].from, b = edges[t].to;
        if (next[a] != -1 || has_in[b])
            continue;
        int ra = uf_find(uf, a), rb = uf_find(uf, b);
        if (ra == rb)
            continue;                       /* would create a cycle */
        next[a] = b;
        has_in[b] = 1;
        uf[ra] = rb;
        taken++;
    }

    /* Walk the single path from the fragment with no incoming edge */
    int start = 0;
    while (start < k && has_in[start])
        start++;
    int *order = xmalloc((size_t)k * sizeof(int));
    for (int pos = 0, cur = start; pos < k; pos++, cur = next[cur])
        order[pos] = cur;

    ReconResult r = recon_result_from_order(frags, k, ov, order, ops);
    free(order);
    free(edges);
    free(next);
    free(has_in);
    free(uf);
    return r;
}
