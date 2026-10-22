/*
 * boyer_moore.c — Boyer-Moore pattern matching with both heuristics.
 *
 * The pattern is compared right to left. On a mismatch at pattern index j,
 * the pattern shifts by the larger of:
 *
 *   Bad character rule:  align the mismatched text character with its last
 *                        occurrence in the pattern  ->  j - last[c]
 *
 *   Good suffix rule:    the matched suffix pat[j+1..m-1] must line up with
 *                        another occurrence of itself in the pattern, or with
 *                        a prefix of the pattern that is a suffix of it
 *                        ->  shift[j+1]  (strong good suffix table)
 *
 * Both shifts are safe, so taking the maximum never skips a match.
 *
 * Time:  O(m + sigma) preprocessing; O(n / m) best, O(n * m) worst search
 * Space: O(m + sigma), sigma = 256 here
 */
#include "search.h"

/* Bad character table: last index of each byte in the pattern, or -1. */
static void build_bad_char(const char *pat, int m, int last[256]) {
    for (int c = 0; c < 256; c++)
        last[c] = -1;
    for (int j = 0; j < m; j++)
        last[(unsigned char)pat[j]] = j;
}

/*
 * Strong good suffix table (shift has m+1 entries).
 * bpos[i] = start of the widest border of the suffix pat[i..m-1].
 *
 * Case 1: the matched suffix occurs elsewhere in the pattern, preceded by a
 *         different character.
 * Case 2: only a part of the matched suffix occurs as a prefix of the pattern.
 */
static void build_good_suffix(const char *pat, int m, int *shift, long long *comps) {
    int *bpos = xmalloc((size_t)(m + 1) * sizeof(int));
    for (int i = 0; i <= m; i++)
        shift[i] = 0;

    /* Case 1 */
    int i = m, j = m + 1;
    bpos[i] = j;
    while (i > 0) {
        while (j <= m) {
            (*comps)++;
            if (pat[i - 1] == pat[j - 1])
                break;
            if (shift[j] == 0)
                shift[j] = j - i;
            j = bpos[j];
        }
        i--;
        j--;
        bpos[i] = j;
    }

    /* Case 2 */
    j = bpos[0];
    for (i = 0; i <= m; i++) {
        if (shift[i] == 0)
            shift[i] = j;
        if (i == j)
            j = bpos[j];
    }
    free(bpos);
}

SearchResult boyer_moore_search(const char *text, int n, const char *pat, int m) {
    SearchResult r = search_result_new();
    if (m <= 0 || m > n)
        return r;

    int last[256];
    int *shift = xmalloc((size_t)(m + 1) * sizeof(int));
    build_bad_char(pat, m, last);
    build_good_suffix(pat, m, shift, &r.comparisons);

    int s = 0;  /* current shift of the pattern over the text */
    while (s <= n - m) {
        int j = m - 1;
        while (j >= 0) {
            r.comparisons++;
            if (pat[j] != text[s + j])
                break;
            j--;
        }
        if (j < 0) {
            search_result_add(&r, s);
            s += shift[0];                  /* shift by the period of the pattern */
        } else {
            int bad = j - last[(unsigned char)text[s + j]];
            int good = shift[j + 1];
            s += MAX(bad, good);
        }
    }
    free(shift);
    return r;
}
