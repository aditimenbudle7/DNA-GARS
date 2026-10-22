/*
 * test_align.c — correctness tests for Module 2.
 * Hand-checked cases plus randomized property checks.
 */
#include "align.h"
#include "seq_io.h"

static int checks = 0, failures = 0;

#define CHECK(cond, ...)                                        \
    do {                                                        \
        checks++;                                               \
        if (!(cond)) {                                          \
            failures++;                                         \
            printf("FAIL %s:%d: ", __FILE__, __LINE__);         \
            printf(__VA_ARGS__);                                \
            printf("\n");                                       \
        }                                                       \
    } while (0)

/* Remove '-' from an aligned string. Caller frees. */
static char *strip_gaps(const char *s) {
    char *d = xmalloc(strlen(s) + 1);
    int k = 0;
    for (; *s; s++)
        if (*s != '-')
            d[k++] = *s;
    d[k] = '\0';
    return d;
}

static char *slice(const char *s, int start, int end) {
    char *d = xmalloc((size_t)(end - start) + 1);
    memcpy(d, s + start, (size_t)(end - start));
    d[end - start] = '\0';
    return d;
}

int main(void) {
    ScoreScheme s = score_default();   /* +1 / -1 / -2 */

    /* Identical sequences */
    {
        AlignResult r = needleman_wunsch("ACGT", "ACGT", s);
        MutationReport m = mutations_detect(&r);
        CHECK(r.score == 4, "NW identical: score %d, expected 4", r.score);
        CHECK(m.count == 0 && m.identity == 100.0, "NW identical: no mutations, 100%% identity");
        mutations_free(&m);
        align_result_free(&r);
    }

    /* One deletion: ACGT vs AGT -> A C G T / A - G T = 1 - 2 + 1 + 1 = 1 */
    {
        AlignResult r = needleman_wunsch("ACGT", "AGT", s);
        MutationReport m = mutations_detect(&r);
        CHECK(r.score == 1, "NW deletion: score %d, expected 1", r.score);
        CHECK(m.deletions == 1 && m.count == 1, "NW deletion: one deletion");
        CHECK(m.count == 1 && m.items[0].ref_pos == 1 && m.items[0].ref_base == 'C',
              "NW deletion: C at reference position 1");
        mutations_free(&m);
        align_result_free(&r);
    }

    /* One substitution: T -> A at position 3. Score 7 - 1 = 6 */
    {
        AlignResult r = needleman_wunsch("ACGTACGT", "ACGAACGT", s);
        MutationReport m = mutations_detect(&r);
        CHECK(r.score == 6, "NW substitution: score %d, expected 6", r.score);
        CHECK(m.substitutions == 1 && m.count == 1, "NW substitution: one substitution");
        CHECK(m.count == 1 && m.items[0].ref_pos == 3 && m.items[0].ref_base == 'T' &&
              m.items[0].alt_base == 'A', "NW substitution: T -> A at 3");
        mutations_free(&m);
        align_result_free(&r);
    }

    /* One insertion: ACGT vs ACCGT */
    {
        AlignResult r = needleman_wunsch("ACGT", "ACCGT", s);
        MutationReport m = mutations_detect(&r);
        CHECK(r.score == 2, "NW insertion: score %d, expected 2", r.score);
        CHECK(m.insertions == 1 && m.count == 1, "NW insertion: one insertion");
        mutations_free(&m);
        align_result_free(&r);
    }

    /* Local alignment finds the shared core ACG */
    {
        AlignResult r = smith_waterman("TTACGTT", "GGACGGG", s);
        CHECK(r.score == 3, "SW: score %d, expected 3", r.score);
        CHECK(r.a_start == 2 && r.a_end == 5, "SW: reference region [%d,%d), expected [2,5)",
              r.a_start, r.a_end);
        CHECK(r.b_start == 2 && r.b_end == 5, "SW: query region [%d,%d), expected [2,5)",
              r.b_start, r.b_end);
        align_result_free(&r);
    }

    /* No similarity at all: local score 0, empty alignment */
    {
        AlignResult r = smith_waterman("AAAA", "CCCC", s);
        CHECK(r.score == 0 && r.length == 0, "SW: unrelated sequences give an empty alignment");
        align_result_free(&r);
    }

    /* Randomized properties */
    unsigned long long st = 99;
    for (int iter = 0; iter < 300; iter++) {
        int n = 1 + (int)(rng_next(&st) % 60);
        char *a = seq_random(n, &st);
        char *b = seq_mutate(a, rng_uniform(&st) * 0.4, &st, NULL, NULL, NULL);

        AlignResult g = needleman_wunsch(a, b, s);
        char *ga = strip_gaps(g.aligned_a), *gb = strip_gaps(g.aligned_b);
        CHECK(strcmp(ga, a) == 0 && strcmp(gb, b) == 0, "NW: alignment must contain both full inputs");
        CHECK(align_rescore(&g, s) == g.score, "NW: rescored alignment must equal DP score");
        CHECK(g.matrix[(size_t)g.rows * g.cols - 1] == g.score, "NW: score is the last DP cell");

        AlignResult l = smith_waterman(a, b, s);
        char *la = strip_gaps(l.aligned_a), *lb = strip_gaps(l.aligned_b);
        char *ra = slice(a, l.a_start, l.a_end), *rb = slice(b, l.b_start, l.b_end);
        CHECK(strcmp(la, ra) == 0 && strcmp(lb, rb) == 0, "SW: alignment must match its region");
        CHECK(align_rescore(&l, s) == l.score, "SW: rescored alignment must equal DP score");
        CHECK(l.score >= 0 && l.score >= g.score, "SW: local score >= max(0, global score)");

        free(ga); free(gb); free(la); free(lb); free(ra); free(rb);
        align_result_free(&g);
        align_result_free(&l);
        free(a);
        free(b);
    }

    printf("test_align: %d checks, %d failure(s)\n", checks, failures);
    return failures ? 1 : 0;
}
