/*
 * naive.c — brute-force pattern matching.
 *
 * Slide the pattern over every shift s = 0..n-m and compare left to right
 * until a mismatch. This is the reference implementation every other
 * algorithm is checked against.
 *
 * Time:  O((n - m + 1) * m) worst case, e.g. text AAAA..., pattern AAA...C
 * Space: O(1) besides the output
 */
#include "search.h"

SearchResult naive_search(const char *text, int n, const char *pat, int m) {
    SearchResult r = search_result_new();
    if (m <= 0 || m > n)
        return r;

    for (int s = 0; s <= n - m; s++) {
        int j = 0;
        while (j < m) {
            r.comparisons++;
            if (text[s + j] != pat[j])
                break;
            j++;
        }
        if (j == m)
            search_result_add(&r, s);
    }
    return r;
}
