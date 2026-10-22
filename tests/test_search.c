/*
 * test_search.c — correctness tests for Module 1.
 * Known cases plus randomized cross-checks of every algorithm against Naive.
 */
#include "benchmark.h"
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

/* Every algorithm must return exactly `expected` (count positions). */
static void check_expected(const char *text, const char *pat, const int *expected, int count) {
    int n = (int)strlen(text), m = (int)strlen(pat);
    for (int a = 0; a < SEARCH_ALGO_COUNT; a++) {
        SearchResult r = SEARCH_ALGOS[a].fn(text, n, pat, m);
        int ok = (r.count == count);
        for (int i = 0; ok && i < count; i++)
            ok = (r.positions[i] == expected[i]);
        CHECK(ok, "%s: '%s' in '%s' expected %d match(es), got %d",
              SEARCH_ALGOS[a].name, pat, text, count, r.count);
        search_result_free(&r);
    }
}

/* Every algorithm must agree with Naive. */
static void check_agree(const char *text, const char *pat) {
    int n = (int)strlen(text), m = (int)strlen(pat);
    SearchResult base = naive_search(text, n, pat, m);
    for (int a = 1; a < SEARCH_ALGO_COUNT; a++) {
        SearchResult r = SEARCH_ALGOS[a].fn(text, n, pat, m);
        CHECK(search_results_equal(&base, &r), "%s disagrees with Naive: pattern '%s' (n=%d)",
              SEARCH_ALGOS[a].name, pat, n);
        search_result_free(&r);
    }
    search_result_free(&base);
}

/* Random string over the first `sigma` letters of ACGT (small sigma = many matches). */
static char *random_over(int n, int sigma, unsigned long long *st) {
    char *s = xmalloc((size_t)n + 1);
    for (int i = 0; i < n; i++)
        s[i] = DNA_ALPHABET[rng_next(st) % (unsigned long long)sigma];
    s[n] = '\0';
    return s;
}

int main(void) {
    /* Example from the project spec */
    const int spec[] = { 3, 15, 19 };
    check_expected("ATGCGTACGATCGATCGTACGTAGCTAGC", "CGT", spec, 3);

    /* Overlapping matches */
    const int overlap[] = { 0, 1, 2, 3, 4 };
    check_expected("AAAAAAA", "AAA", overlap, 5);

    /* Pattern equal to text, longer than text, absent, single base */
    const int zero[] = { 0 };
    check_expected("ACGT", "ACGT", zero, 1);
    check_expected("ACG", "ACGT", NULL, 0);
    check_expected("ACGTACGT", "TTT", NULL, 0);
    const int single[] = { 1, 5 };
    check_expected("ACGTACGT", "C", single, 2);

    /* Periodic patterns stress the KMP and good-suffix tables */
    check_agree("ACACACACACACAC", "ACAC");
    check_agree("AACAACAACAACAAC", "AACAAC");
    check_agree("GGGGGGGGGGGGGGGGGG", "GGGG");

    /* Randomized cross-check */
    unsigned long long st = 2026;
    for (int iter = 0; iter < 2000; iter++) {
        int sigma = 1 + (int)(rng_next(&st) % 4);
        int n = 1 + (int)(rng_next(&st) % 300);
        int m = 1 + (int)(rng_next(&st) % 10);
        char *text = random_over(n, sigma, &st);
        char *pat = random_over(m, sigma, &st);
        check_agree(text, pat);
        free(text);
        free(pat);
    }

    /* Rabin-Karp counters are consistent */
    {
        const char *t = "ATGCGTACGATCGATCGTACGTAGCTAGC";
        SearchResult r = rabin_karp_search(t, (int)strlen(t), "CGT", 3);
        CHECK(r.hash_checks == (long long)strlen(t) - 3 + 1, "RK hash checks = n - m + 1");
        CHECK(r.spurious_hits == 0, "RK no spurious hits for a short pattern");
        search_result_free(&r);
    }

    printf("test_search: %d checks, %d failure(s)\n", checks, failures);
    return failures ? 1 : 0;
}
