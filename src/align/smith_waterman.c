/*
 * smith_waterman.c — local alignment (Smith-Waterman, 1981).
 *
 * Same recurrence as Needleman-Wunsch with two changes:
 *   1. Every cell is floored at 0, so a bad prefix never drags a good region down:
 *        H[i][j] = max(0, H[i-1][j-1] + sub, H[i-1][j] + gap, H[i][j-1] + gap)
 *   2. First row/column are 0, and the answer is the maximum cell anywhere.
 * Traceback starts at that maximum and stops at the first 0 cell, giving the
 * best-scoring pair of substrings.
 *
 * Time:  O(n * m)
 * Space: O(n * m)
 */
#include "align.h"

#define AT(i, j) dp[(size_t)(i) * (size_t)cols + (size_t)(j)]

static inline int sub_score(char x, char y, ScoreScheme s) {
    return x == y ? s.match : s.mismatch;
}

static void reverse(char *s, int len) {
    for (int i = 0, j = len - 1; i < j; i++, j--) {
        char t = s[i];
        s[i] = s[j];
        s[j] = t;
    }
}

AlignResult smith_waterman(const char *a, const char *b, ScoreScheme s) {
    AlignResult r;
    memset(&r, 0, sizeof r);

    int n = (int)strlen(a), m = (int)strlen(b);
    int rows = n + 1, cols = m + 1;
    int *dp = xmalloc((size_t)rows * (size_t)cols * sizeof(int));

    /* Fill, remembering the first maximum cell */
    int best = 0, bi = 0, bj = 0;
    for (int i = 0; i <= n; i++)
        AT(i, 0) = 0;
    for (int j = 0; j <= m; j++)
        AT(0, j) = 0;
    for (int i = 1; i <= n; i++) {
        for (int j = 1; j <= m; j++) {
            int diag = AT(i - 1, j - 1) + sub_score(a[i - 1], b[j - 1], s);
            int up = AT(i - 1, j) + s.gap;
            int left = AT(i, j - 1) + s.gap;
            int v = MAX(0, MAX(diag, MAX(up, left)));
            AT(i, j) = v;
            if (v > best) {
                best = v;
                bi = i;
                bj = j;
            }
        }
    }

    /* Traceback from the best cell until a 0 cell */
    char *ra = xmalloc((size_t)(n + m) + 1);
    char *rb = xmalloc((size_t)(n + m) + 1);
    int len = 0, i = bi, j = bj;
    while (i > 0 && j > 0 && AT(i, j) > 0) {
        if (AT(i, j) == AT(i - 1, j - 1) + sub_score(a[i - 1], b[j - 1], s)) {
            ra[len] = a[i - 1];
            rb[len] = b[j - 1];
            i--;
            j--;
        } else if (AT(i, j) == AT(i - 1, j) + s.gap) {
            ra[len] = a[i - 1];
            rb[len] = '-';
            i--;
        } else {
            ra[len] = '-';
            rb[len] = b[j - 1];
            j--;
        }
        len++;
    }
    reverse(ra, len);
    reverse(rb, len);
    ra[len] = rb[len] = '\0';

    r.aligned_a = ra;
    r.aligned_b = rb;
    r.length = len;
    r.score = best;
    r.a_start = i;
    r.a_end = bi;
    r.b_start = j;
    r.b_end = bj;
    r.matrix = dp;
    r.rows = rows;
    r.cols = cols;
    r.cells = (long long)n * m;
    return r;
}
