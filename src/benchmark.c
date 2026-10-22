/*
 * benchmark.c — algorithm registry, timed runners and CSV benchmarks.
 *
 * All numbers in the CSVs are measured by the program. Inputs are generated
 * with fixed seeds, so every run benchmarks the same data.
 */
#include "benchmark.h"
#include "seq_io.h"
#include "timer.h"

/* ------------------------------------------------------------------ */
/* Registry                                                            */
/* ------------------------------------------------------------------ */

const SearchAlgo SEARCH_ALGOS[] = {
    { "Naive",       naive_search },
    { "KMP",         kmp_search },
    { "Rabin-Karp",  rabin_karp_search },
    { "Boyer-Moore", boyer_moore_search },
    { "Z-Algorithm", z_search },
};
const int SEARCH_ALGO_COUNT = (int)(sizeof SEARCH_ALGOS / sizeof SEARCH_ALGOS[0]);

int search_results_equal(const SearchResult *a, const SearchResult *b) {
    if (a->count != b->count)
        return 0;
    for (int i = 0; i < a->count; i++)
        if (a->positions[i] != b->positions[i])
            return 0;
    return 1;
}

/* ------------------------------------------------------------------ */
/* Timed runners                                                       */
/* ------------------------------------------------------------------ */

/*
 * TIME_IT runs `stmt` in batches. The batch size doubles until one batch takes
 * at least MIN_TIMING_MS (so tiny inputs are measurable), then TIMING_BATCHES
 * batches are timed and the fastest per-run time is kept. The minimum is the
 * least noisy estimate: other processes can only make a run slower.
 * `stmt` must free whatever it allocates.
 */
#define TIME_IT(out_ms, stmt)                                               \
    do {                                                                    \
        int batch_ = 1;                                                     \
        double t_;                                                          \
        for (;;) {                                                          \
            double s_ = timer_now_ms();                                     \
            for (int q_ = 0; q_ < batch_; q_++) { stmt; }                   \
            t_ = timer_now_ms() - s_;                                       \
            if (t_ >= MIN_TIMING_MS || batch_ >= MAX_TIMING_RUNS)           \
                break;                                                      \
            batch_ *= 2;                                                    \
        }                                                                   \
        double best_ = t_ / batch_;                                         \
        for (int b_ = 1; b_ < TIMING_BATCHES; b_++) {                       \
            double s_ = timer_now_ms();                                     \
            for (int q_ = 0; q_ < batch_; q_++) { stmt; }                   \
            t_ = (timer_now_ms() - s_) / batch_;                            \
            if (t_ < best_)                                                 \
                best_ = t_;                                                 \
        }                                                                   \
        (out_ms) = best_;                                                   \
    } while (0)

SearchResult search_run(SearchFn fn, const char *text, int n, const char *pat, int m) {
    SearchResult r = fn(text, n, pat, m);
    TIME_IT(r.time_ms, { SearchResult tmp = fn(text, n, pat, m); search_result_free(&tmp); });
    return r;
}

AlignResult align_run(AlignFn fn, const char *a, const char *b, ScoreScheme s) {
    AlignResult r = fn(a, b, s);
    TIME_IT(r.time_ms, { AlignResult tmp = fn(a, b, s); align_result_free(&tmp); });
    return r;
}

ReconResult recon_run(ReconFn fn, char **frags, int k, int **ov) {
    ReconResult r = fn(frags, k, ov);
    TIME_IT(r.time_ms, { ReconResult tmp = fn(frags, k, ov); recon_result_free(&tmp); });
    return r;
}

/* ------------------------------------------------------------------ */
/* Module 1: pattern search                                            */
/* ------------------------------------------------------------------ */

/*
 * Two scenarios per genome length n (pattern length 16):
 *   random      random genome, pattern copied from its middle (>= 1 match)
 *   repetitive  genome AAAA...A, pattern AAA...AC: worst case for Naive,
 *               zero matches, every shift nearly matches
 */
int bench_search(const char *csv_path) {
    static const int sizes[] = { 1000, 5000, 10000, 50000, 100000, 500000, 1000000 };
    const int nsizes = (int)(sizeof sizes / sizeof sizes[0]);
    const int m = 16;
    static const char *scenarios[] = { "random", "repetitive" };

    FILE *f = fopen(csv_path, "w");
    if (!f) {
        fprintf(stderr, "  Cannot write '%s' (does the folder exist?)\n", csv_path);
        return -1;
    }
    fprintf(f, "scenario,algorithm,n,m,matches,comparisons,hash_checks,time_ms\n");

    printf("\n[Pattern search]  m = %d\n", m);
    printf("%-11s %-12s %9s %8s %14s %12s\n",
           "Scenario", "Algorithm", "n", "Matches", "Comparisons", "Time (ms)");

    for (int sc = 0; sc < 2; sc++) {
        for (int si = 0; si < nsizes; si++) {
            int n = sizes[si];
            char *text, *pat;
            if (sc == 0) {
                unsigned long long st = 1000003ULL * (unsigned long long)n;
                text = seq_random(n, &st);
                pat = xmalloc((size_t)m + 1);
                memcpy(pat, text + n / 2, (size_t)m);
                pat[m] = '\0';
            } else {
                text = xmalloc((size_t)n + 1);
                memset(text, 'A', (size_t)n);
                text[n] = '\0';
                pat = xmalloc((size_t)m + 1);
                memset(pat, 'A', (size_t)m - 1);
                pat[m - 1] = 'C';
                pat[m] = '\0';
            }

            for (int a = 0; a < SEARCH_ALGO_COUNT; a++) {
                SearchResult r = search_run(SEARCH_ALGOS[a].fn, text, n, pat, m);
                fprintf(f, "%s,%s,%d,%d,%d,%lld,%lld,%.6f\n", scenarios[sc],
                        SEARCH_ALGOS[a].name, n, m, r.count, r.comparisons, r.hash_checks,
                        r.time_ms);
                printf("%-11s %-12s %9d %8d %14lld %12.4f\n", scenarios[sc],
                       SEARCH_ALGOS[a].name, n, r.count, r.comparisons, r.time_ms);
                search_result_free(&r);
            }
            free(text);
            free(pat);
        }
    }
    fclose(f);
    printf("Saved %s\n", csv_path);
    return 0;
}

/* ------------------------------------------------------------------ */
/* Module 2: alignment                                                 */
/* ------------------------------------------------------------------ */

/* Reference: random genome of length n. Query: copy with 10% mutations. */
int bench_align(const char *csv_path) {
    static const int sizes[] = { 100, 200, 400, 800, 1600, 3200 };
    const int nsizes = (int)(sizeof sizes / sizeof sizes[0]);
    static const struct { const char *name; AlignFn fn; } algos[] = {
        { "Needleman-Wunsch", needleman_wunsch },
        { "Smith-Waterman",   smith_waterman },
    };
    ScoreScheme s = score_default();

    FILE *f = fopen(csv_path, "w");
    if (!f) {
        fprintf(stderr, "  Cannot write '%s' (does the folder exist?)\n", csv_path);
        return -1;
    }
    fprintf(f, "algorithm,n,m,cells,score,identity,time_ms\n");

    printf("\n[Alignment]  query = reference with 10%% mutations\n");
    printf("%-17s %6s %6s %12s %7s %9s %12s\n",
           "Algorithm", "n", "m", "Cells", "Score", "Identity", "Time (ms)");

    for (int si = 0; si < nsizes; si++) {
        int n = sizes[si];
        unsigned long long st = 7919ULL * (unsigned long long)n;
        char *a = seq_random(n, &st);
        char *b = seq_mutate(a, 0.10, &st, NULL, NULL, NULL);
        int m = (int)strlen(b);

        for (int k = 0; k < 2; k++) {
            AlignResult r = align_run(algos[k].fn, a, b, s);
            MutationReport mr = mutations_detect(&r);
            fprintf(f, "%s,%d,%d,%lld,%d,%.4f,%.6f\n", algos[k].name, n, m,
                    r.cells, r.score, mr.identity, r.time_ms);
            printf("%-17s %6d %6d %12lld %7d %8.2f%% %12.4f\n", algos[k].name, n, m,
                   r.cells, r.score, mr.identity, r.time_ms);
            mutations_free(&mr);
            align_result_free(&r);
        }
        free(a);
        free(b);
    }
    fclose(f);
    printf("Saved %s\n", csv_path);
    return 0;
}

/* ------------------------------------------------------------------ */
/* Module 3: reconstruction                                            */
/* ------------------------------------------------------------------ */

/*
 * For k = 3..DP_BITMASK_MAX_K: k reads of length 10 at random positions of a
 * random genome (about 2x coverage), so overlaps vary and greedy can miss the
 * optimum. Duplicate reads are redrawn, so each k uses exactly k distinct reads.
 * Brute force only runs while k <= BRUTE_FORCE_MAX_K.
 */
int bench_reconstruct(const char *csv_path) {
    const int read_len = 10;

    FILE *f = fopen(csv_path, "w");
    if (!f) {
        fprintf(stderr, "  Cannot write '%s' (does the folder exist?)\n", csv_path);
        return -1;
    }
    fprintf(f, "algorithm,k,superstring_len,operations,time_ms\n");

    printf("\n[Reconstruction]  reads of length %d at random positions\n", read_len);
    printf("%-12s %4s %8s %14s %12s\n", "Algorithm", "k", "Length", "Operations", "Time (ms)");

    for (int target = 3; target <= DP_BITMASK_MAX_K; target++) {
        unsigned long long st = 104729ULL * (unsigned long long)target;
        int glen = 5 * target + read_len;
        char *genome = seq_random(glen, &st);

        /* Draw reads until there are exactly `target` distinct ones */
        char **frags = xmalloc((size_t)target * sizeof(char *));
        int k = 0;
        while (k < target) {
            int start = (int)(rng_next(&st) % (unsigned long long)(glen - read_len + 1));
            frags[k] = xmalloc((size_t)read_len + 1);
            memcpy(frags[k], genome + start, (size_t)read_len);
            frags[k][read_len] = '\0';
            k = fragments_remove_contained(frags, k + 1);
        }
        int **ov = overlap_matrix(frags, k, NULL);

        static const struct { const char *name; ReconFn fn; int max_k; } algos[] = {
            { "BruteForce", scs_brute_force, BRUTE_FORCE_MAX_K },
            { "DP-Bitmask", scs_dp_bitmask,  DP_BITMASK_MAX_K },
            { "Greedy",     scs_greedy,      1000000 },
        };
        for (int a = 0; a < 3; a++) {
            if (k > algos[a].max_k)
                continue;
            ReconResult r = recon_run(algos[a].fn, frags, k, ov);
            fprintf(f, "%s,%d,%d,%lld,%.6f\n", algos[a].name, k, r.length, r.operations, r.time_ms);
            printf("%-12s %4d %8d %14lld %12.4f\n", algos[a].name, k, r.length, r.operations, r.time_ms);
            recon_result_free(&r);
        }

        overlap_matrix_free(ov, k);
        seq_free_list(frags, k);
        free(genome);
    }
    fclose(f);
    printf("Saved %s\n", csv_path);
    return 0;
}

int bench_all(const char *dir) {
    char path[1024];
    int rc = 0;

    snprintf(path, sizeof path, "%s/search.csv", dir);
    rc |= bench_search(path);
    snprintf(path, sizeof path, "%s/align.csv", dir);
    rc |= bench_align(path);
    snprintf(path, sizeof path, "%s/reconstruct.csv", dir);
    rc |= bench_reconstruct(path);
    return rc ? -1 : 0;
}
