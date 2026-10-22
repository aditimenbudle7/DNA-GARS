/*
 * overlap.c — overlap computation and shared helpers for genome reconstruction.
 *
 * overlap(a, b) = length of the longest suffix of a that equals a prefix of b.
 * Computed with the KMP prefix function (Module 1) on  b + '#' + a :
 * the final lps value is the longest prefix of b that ends the string, i.e.
 * the longest prefix of b that is a suffix of a. '#' stops matches crossing over.
 *
 * Time:  O(|a| + |b|) per pair, O(k^2 * L) for the whole matrix
 */
#include "reconstruct.h"
#include "search.h"

int overlap_length(const char *a, const char *b, long long *ops) {
    int la = (int)strlen(a), lb = (int)strlen(b);
    int len = lb + 1 + la;
    char *s = xmalloc((size_t)len + 1);
    memcpy(s, b, (size_t)lb);
    s[lb] = '#';
    memcpy(s + lb + 1, a, (size_t)la);
    s[len] = '\0';

    int *lps = xmalloc((size_t)len * sizeof(int));
    kmp_prefix_function(s, len, lps, ops);
    int ov = lps[len - 1];

    free(lps);
    free(s);
    return ov;
}

int **overlap_matrix(char **frags, int k, long long *ops) {
    int **ov = xmalloc((size_t)k * sizeof(int *));
    for (int i = 0; i < k; i++) {
        ov[i] = xcalloc((size_t)k, sizeof(int));
        for (int j = 0; j < k; j++)
            if (i != j)
                ov[i][j] = overlap_length(frags[i], frags[j], ops);
    }
    return ov;
}

void overlap_matrix_free(int **ov, int k) {
    if (!ov)
        return;
    for (int i = 0; i < k; i++)
        free(ov[i]);
    free(ov);
}

/* A fragment inside another adds no information to the superstring and breaks
 * the "merge by overlap" model, so it is removed first. Uses KMP for the test. */
int fragments_remove_contained(char **frags, int k) {
    int *drop = xcalloc((size_t)k, sizeof(int));
    for (int i = 0; i < k; i++) {
        int li = (int)strlen(frags[i]);
        for (int j = 0; j < k && !drop[i]; j++) {
            if (i == j || drop[j])
                continue;
            int lj = (int)strlen(frags[j]);
            if (li > lj)
                continue;
            if (li == lj) {
                if (strcmp(frags[i], frags[j]) == 0)
                    drop[i] = 1;            /* duplicate: keep the other copy */
                continue;
            }
            SearchResult r = kmp_search(frags[j], lj, frags[i], li);
            if (r.count > 0)
                drop[i] = 1;
            search_result_free(&r);
        }
    }

    int w = 0;
    for (int i = 0; i < k; i++) {
        if (drop[i])
            free(frags[i]);
        else
            frags[w++] = frags[i];
    }
    free(drop);
    return w;
}

char *superstring_from_order(char **frags, int k, int **ov, const int *order) {
    if (k == 0)
        return xstrdup("");

    size_t total = strlen(frags[order[0]]);
    for (int t = 1; t < k; t++)
        total += strlen(frags[order[t]]) - (size_t)ov[order[t - 1]][order[t]];

    char *s = xmalloc(total + 1);
    size_t pos = strlen(frags[order[0]]);
    memcpy(s, frags[order[0]], pos);
    for (int t = 1; t < k; t++) {
        const char *f = frags[order[t]];
        int skip = ov[order[t - 1]][order[t]];
        size_t add = strlen(f) - (size_t)skip;
        memcpy(s + pos, f + skip, add);
        pos += add;
    }
    s[pos] = '\0';
    return s;
}

ReconResult recon_result_from_order(char **frags, int k, int **ov, const int *order,
                                    long long operations) {
    ReconResult r;
    memset(&r, 0, sizeof r);
    r.k = k;
    r.order = xmalloc((size_t)(k ? k : 1) * sizeof(int));
    if (k)
        memcpy(r.order, order, (size_t)k * sizeof(int));
    r.superstring = superstring_from_order(frags, k, ov, order);
    r.length = (int)strlen(r.superstring);
    r.operations = operations;
    return r;
}

void recon_result_free(ReconResult *r) {
    free(r->superstring);
    free(r->order);
    r->superstring = NULL;
    r->order = NULL;
}
