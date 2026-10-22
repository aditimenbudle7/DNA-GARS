/*
 * kmp.c — Knuth-Morris-Pratt pattern matching.
 *
 * Preprocess the pattern into the prefix function lps[], where lps[i] is the
 * length of the longest proper prefix of pat[0..i] that is also a suffix of it.
 * On a mismatch after j matched characters, the pattern jumps to lps[j-1]
 * instead of restarting, so the text pointer never moves backwards.
 *
 * Time:  O(m) preprocessing + O(n) search; at most 2n comparisons in the search
 * Space: O(m)
 */
#include "search.h"

void kmp_prefix_function(const char *p, int m, int *lps, long long *comparisons) {
    long long comps = 0;
    if (m <= 0)
        return;

    lps[0] = 0;
    int len = 0;    /* length of the current matched prefix */
    int i = 1;
    while (i < m) {
        comps++;
        if (p[i] == p[len]) {
            lps[i++] = ++len;
        } else if (len > 0) {
            len = lps[len - 1];     /* fall back, do not advance i */
        } else {
            lps[i++] = 0;
        }
    }
    if (comparisons)
        *comparisons += comps;
}

SearchResult kmp_search(const char *text, int n, const char *pat, int m) {
    SearchResult r = search_result_new();
    if (m <= 0 || m > n)
        return r;

    int *lps = xmalloc((size_t)m * sizeof(int));
    kmp_prefix_function(pat, m, lps, &r.comparisons);

    int i = 0, j = 0;   /* i: text index, j: number of pattern chars matched */
    while (i < n) {
        r.comparisons++;
        if (text[i] == pat[j]) {
            i++;
            j++;
            if (j == m) {
                search_result_add(&r, i - m);
                j = lps[j - 1];     /* allow overlapping matches */
            }
        } else if (j > 0) {
            j = lps[j - 1];
        } else {
            i++;
        }
    }
    free(lps);
    return r;
}
