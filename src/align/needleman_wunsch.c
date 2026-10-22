/*
 * needleman_wunsch.c — global alignment (Needleman-Wunsch, 1970).
 *
 * DP[i][j] = best score aligning a[0..i) with b[0..j)
 *
 *   DP[0][j] = j * gap        DP[i][0] = i * gap
 *   DP[i][j] = max( DP[i-1][j-1] + (a[i-1] == b[j-1] ? match : mismatch),
 *                   DP[i-1][j]   + gap,      (a[i-1] aligned to a gap: deletion)
 *                   DP[i][j-1]   + gap )     (b[j-1] aligned to a gap: insertion)
 *
 * The answer is DP[n][m]; the alignment is recovered by tracing back from
 * (n, m) to (0, 0), preferring diagonal, then up, then left on ties.
 *
 * Time:  O(n * m)
 * Space: O(n * m) (full matrix kept for traceback and display)
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

AlignResult needleman_wunsch(const char *a, const char *b, ScoreScheme s) {
    AlignResult r;
    memset(&r, 0, sizeof r);

    int n = (int)strlen(a), m = (int)strlen(b);
    int rows = n + 1, cols = m + 1;
    int *dp = xmalloc((size_t)rows * (size_t)cols * sizeof(int));

    /* Fill */
    for (int i = 0; i <= n; i++)
        AT(i, 0) = i * s.gap;
    for (int j = 0; j <= m; j++)
        AT(0, j) = j * s.gap;
    for (int i = 1; i <= n; i++) {
        for (int j = 1; j <= m; j++) {
            int diag = AT(i - 1, j - 1) + sub_score(a[i - 1], b[j - 1], s);
            int up = AT(i - 1, j) + s.gap;
            int left = AT(i, j - 1) + s.gap;
            AT(i, j) = MAX(diag, MAX(up, left));
        }
    }

    /* Traceback */
    char *ra = xmalloc((size_t)(n + m) + 1);
    char *rb = xmalloc((size_t)(n + m) + 1);
    int len = 0, i = n, j = m;
    while (i > 0 || j > 0) {
        if (i > 0 && j > 0 && AT(i, j) == AT(i - 1, j - 1) + sub_score(a[i - 1], b[j - 1], s)) {
            ra[len] = a[i - 1];
            rb[len] = b[j - 1];
            i--;
            j--;
        } else if (i > 0 && AT(i, j) == AT(i - 1, j) + s.gap) {
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
    r.score = AT(n, m);
    r.a_start = 0;
    r.a_end = n;
    r.b_start = 0;
    r.b_end = m;
    r.matrix = dp;
    r.rows = rows;
    r.cols = cols;
    r.cells = (long long)n * m;
    return r;
}
