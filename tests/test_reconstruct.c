/*
 * test_reconstruct.c — correctness tests for Module 3.
 * Exact methods must agree; every result must contain every fragment.
 */
#include "reconstruct.h"
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

static int contains_all(const char *super, char **frags, int k) {
    for (int i = 0; i < k; i++)
        if (!strstr(super, frags[i]))
            return 0;
    return 1;
}

static char **copy_list(const char *const *src, int k) {
    char **d = xmalloc((size_t)k * sizeof(char *));
    for (int i = 0; i < k; i++)
        d[i] = xstrdup(src[i]);
    return d;
}

/* Run all three algorithms on a fragment set and check the invariants. */
static void check_set(char **frags, int k, const char *label) {
    k = fragments_remove_contained(frags, k);
    int **ov = overlap_matrix(frags, k, NULL);

    ReconResult bf = scs_brute_force(frags, k, ov);
    ReconResult dp = scs_dp_bitmask(frags, k, ov);
    ReconResult gr = scs_greedy(frags, k, ov);

    CHECK(contains_all(bf.superstring, frags, k), "%s: brute force misses a fragment", label);
    CHECK(contains_all(dp.superstring, frags, k), "%s: DP misses a fragment", label);
    CHECK(contains_all(gr.superstring, frags, k), "%s: greedy misses a fragment", label);
    CHECK(bf.length == dp.length, "%s: brute force (%d) and DP (%d) must agree", label,
          bf.length, dp.length);
    CHECK(gr.length >= dp.length, "%s: greedy (%d) cannot beat the optimum (%d)", label,
          gr.length, dp.length);

    recon_result_free(&bf);
    recon_result_free(&dp);
    recon_result_free(&gr);
    overlap_matrix_free(ov, k);
    seq_free_list(frags, k);
}

int main(void) {
    /* Overlap via KMP prefix function */
    CHECK(overlap_length("ATGCGTAC", "GTACGATC", NULL) == 4, "overlap ATGCGTAC -> GTACGATC = 4");
    CHECK(overlap_length("AAAA", "CCCC", NULL) == 0, "overlap with nothing in common = 0");
    CHECK(overlap_length("ACGTAC", "TACGGG", NULL) == 3, "overlap ACGTAC -> TACGGG = 3 (TAC)");

    /* Containment and duplicate removal */
    {
        const char *src[] = { "ACGTACGT", "GTAC", "ACGTACGT", "TTTT" };
        char **f = copy_list(src, 4);
        int k = fragments_remove_contained(f, 4);
        CHECK(k == 2, "remove_contained: expected 2 fragments left, got %d", k);
        seq_free_list(f, k);
    }

    /* Fragments from data/fragments.txt (genome from the project spec, 29 bases) */
    {
        const char *src[] = { "ATGCGTAC", "GTACGATC", "GATCGATC", "GATCGTAC",
                              "GTACGTAG", "GTAGCTAG", "GCTAGC" };
        char **f = copy_list(src, 7);
        int k = fragments_remove_contained(f, 7);
        int **ov = overlap_matrix(f, k, NULL);
        ReconResult dp = scs_dp_bitmask(f, k, ov);
        CHECK(dp.length <= 29, "spec fragments: optimum (%d) <= original length 29", dp.length);
        recon_result_free(&dp);
        overlap_matrix_free(ov, k);
        seq_free_list(f, k);

        check_set(copy_list(src, 7), 7, "spec fragments");
    }

    /* Clean tiling of a random genome is reconstructed exactly */
    {
        unsigned long long st = 7;
        char *g = seq_random(80, &st);
        int k = 0;
        char **f = seq_fragment(g, 20, 6, &k);
        k = fragments_remove_contained(f, k);
        int **ov = overlap_matrix(f, k, NULL);
        ReconResult dp = scs_dp_bitmask(f, k, ov);
        CHECK(dp.length <= 80 && contains_all(dp.superstring, f, k),
              "tiled genome: optimum within original length");
        recon_result_free(&dp);
        overlap_matrix_free(ov, k);
        seq_free_list(f, k);
        free(g);
    }

    /* Randomized: brute force == DP, greedy >= DP */
    unsigned long long st = 31337;
    for (int iter = 0; iter < 200; iter++) {
        int k = 1 + (int)(rng_next(&st) % 7);
        int glen = 12 + (int)(rng_next(&st) % 30);
        char *g = seq_random(glen, &st);
        char **f = xmalloc((size_t)k * sizeof(char *));
        for (int i = 0; i < k; i++) {
            int len = 3 + (int)(rng_next(&st) % 8);
            if (len > glen)
                len = glen;
            int start = (int)(rng_next(&st) % (unsigned long long)(glen - len + 1));
            f[i] = xmalloc((size_t)len + 1);
            memcpy(f[i], g + start, (size_t)len);
            f[i][len] = '\0';
        }
        check_set(f, k, "random set");
        free(g);
    }

    printf("test_reconstruct: %d checks, %d failure(s)\n", checks, failures);
    return failures ? 1 : 0;
}
