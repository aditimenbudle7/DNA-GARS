/*
 * z_algo.c — pattern matching with the Z-algorithm.
 *
 * Build S = pattern + '$' + text. Z[i] is the length of the longest substring
 * starting at i that matches a prefix of S. The separator '$' never appears in
 * DNA, so Z[i] <= m, and Z[i] == m exactly where the pattern occurs.
 *
 * The algorithm keeps a window [l, r) = rightmost segment known to match a
 * prefix. For i inside it, Z[i - l] is reused, so each character is extended
 * past r at most once.
 *
 * Time:  O(n + m)
 * Space: O(n + m)
 */
#include "search.h"

SearchResult z_search(const char *text, int n, const char *pat, int m) {
    SearchResult r = search_result_new();
    if (m <= 0 || m > n)
        return r;

    int len = m + 1 + n;
    char *s = xmalloc((size_t)len + 1);
    memcpy(s, pat, (size_t)m);
    s[m] = '$';
    memcpy(s + m + 1, text, (size_t)n);
    s[len] = '\0';

    int *z = xcalloc((size_t)len, sizeof(int));
    int l = 0, rr = 0;  /* [l, rr) is the current Z-box */

    for (int i = 1; i < len; i++) {
        if (i < rr && z[i - l] < rr - i) {
            z[i] = z[i - l];            /* fully inside the box: copy, no comparisons */
        } else {
            z[i] = (i < rr) ? rr - i : 0;
            while (i + z[i] < len) {    /* extend past the box */
                r.comparisons++;
                if (s[z[i]] != s[i + z[i]])
                    break;
                z[i]++;
            }
            l = i;
            rr = i + z[i];
        }
        if (i > m && z[i] >= m)
            search_result_add(&r, i - m - 1);
    }

    free(z);
    free(s);
    return r;
}
