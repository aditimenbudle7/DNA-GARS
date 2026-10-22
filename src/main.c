/*
 * main.c — DNA-GARS: DNA Genome Analysis and Reconstruction System.
 *
 * Menu-driven CLI for the three modules:
 *   1. Pattern search        Naive, KMP, Rabin-Karp, Boyer-Moore, Z-Algorithm
 *   2. Sequence alignment    Needleman-Wunsch, Smith-Waterman, mutation report
 *   3. Genome reconstruction Brute force, bitmask DP, greedy (Shortest Common Superstring)
 *   4. Benchmarks            CSV output for scripts/plot.py
 *
 * All positions are 0-based.
 */
#include <ctype.h>

#include "common.h"
#include "timer.h"
#include "seq_io.h"
#include "search.h"
#include "align.h"
#include "reconstruct.h"
#include "benchmark.h"

#define MAX_POSITIONS_SHOWN 50
#define MAX_MUTATIONS_SHOWN 30
#define MAX_ALIGN_BLOCKS    10
#define MAX_SEQ_SHOWN       200
#define MAX_LIST_SHOWN      12

/* ------------------------------------------------------------------ */
/* Input helpers                                                       */
/* ------------------------------------------------------------------ */

static char *prompt_line(const char *msg) {
    printf("%s", msg);
    fflush(stdout);
    char *s = seq_read_line(stdin);
    if (!s) {               /* EOF (Ctrl+D / end of piped input) */
        printf("\n");
        exit(0);
    }
    return s;
}

static char *trim(char *s) {
    while (isspace((unsigned char)*s))
        s++;
    char *e = s + strlen(s);
    while (e > s && isspace((unsigned char)e[-1]))
        *--e = '\0';
    return s;
}

static int prompt_int(const char *msg, int def) {
    char *line = prompt_line(msg);
    char *s = trim(line);
    int v = def;
    if (*s) {
        char *end;
        long x = strtol(s, &end, 10);
        if (*end == '\0')
            v = (int)x;
        else
            printf("  Not a number, using %d\n", def);
    }
    free(line);
    return v;
}

static double prompt_double(const char *msg, double def) {
    char *line = prompt_line(msg);
    char *s = trim(line);
    double v = def;
    if (*s) {
        char *end;
        double x = strtod(s, &end);
        if (*end == '\0')
            v = x;
        else
            printf("  Not a number, using %g\n", def);
    }
    free(line);
    return v;
}

static int prompt_yes(const char *msg, int def_yes) {
    char *line = prompt_line(msg);
    char c = *trim(line);
    free(line);
    if (!c)
        return def_yes;
    return c == 'y' || c == 'Y';
}

static char *prompt_path(const char *msg) {
    char *line = prompt_line(msg);
    char *s = xstrdup(trim(line));
    free(line);
    return s;
}

static unsigned long long prompt_seed(void) {
    return (unsigned long long)prompt_int("  Random seed [42]: ", 42);
}

static void print_truncated(const char *label, const char *s) {
    int n = (int)strlen(s);
    if (n <= MAX_SEQ_SHOWN)
        printf("%s%s\n", label, s);
    else
        printf("%s%.*s... (%d bases)\n", label, MAX_SEQ_SHOWN, s, n);
}

/*
 * Ask for a sequence: typed, from a file, random, or (if reference != NULL)
 * a mutated copy of the reference. Returns a validated sequence or NULL.
 */
static char *input_sequence(const char *what, const char *reference) {
    printf("\n%s input:\n", what);
    printf("  1. Type it\n");
    printf("  2. Load from file (plain text or FASTA)\n");
    printf("  3. Generate random\n");
    if (reference)
        printf("  4. Mutated copy of the reference\n");
    int c = prompt_int("Choice: ", 1);

    char *seq = NULL;
    switch (c) {
    case 1:
        seq = prompt_line("  Sequence: ");
        seq_normalize(seq);
        if (!seq_is_valid(seq)) {
            printf("  Invalid sequence: use only A, C, G, T.\n");
            free(seq);
            seq = NULL;
        }
        break;
    case 2: {
        char *path = prompt_path("  File path: ");
        seq = seq_read_file(path);
        free(path);
        break;
    }
    case 3: {
        int n = prompt_int("  Length [1000]: ", 1000);
        if (n < 1) {
            printf("  Length must be at least 1.\n");
            break;
        }
        unsigned long long st = prompt_seed();
        seq = seq_random(n, &st);
        break;
    }
    case 4: {
        if (!reference) {
            printf("  Invalid choice.\n");
            break;
        }
        double rate = prompt_double("  Mutation rate per base, 0 to 1 [0.05]: ", 0.05);
        if (rate < 0 || rate > 1) {
            printf("  Rate must be between 0 and 1.\n");
            break;
        }
        unsigned long long st = prompt_seed();
        int s, i, d;
        seq = seq_mutate(reference, rate, &st, &s, &i, &d);
        printf("  Introduced %d substitution(s), %d insertion(s), %d deletion(s)\n", s, i, d);
        break;
    }
    default:
        printf("  Invalid choice.\n");
    }

    if (seq) {
        printf("  %s length: %d\n", what, (int)strlen(seq));
        print_truncated("  ", seq);
    }
    return seq;
}

/* ------------------------------------------------------------------ */
/* 1. Pattern search                                                   */
/* ------------------------------------------------------------------ */

static void menu_search(void) {
    char *text = input_sequence("Genome", NULL);
    if (!text)
        return;
    char *pat = prompt_line("\nPattern: ");
    seq_normalize(pat);

    int n = (int)strlen(text), m = (int)strlen(pat);
    if (!seq_is_valid(pat)) {
        printf("Invalid pattern: use only A, C, G, T.\n");
        free(text);
        free(pat);
        return;
    }
    if (m > n) {
        printf("The pattern (%d) is longer than the genome (%d).\n", m, n);
        free(text);
        free(pat);
        return;
    }

    SearchResult *res = xmalloc((size_t)SEARCH_ALGO_COUNT * sizeof(SearchResult));
    for (int a = 0; a < SEARCH_ALGO_COUNT; a++)
        res[a] = search_run(SEARCH_ALGOS[a].fn, text, n, pat, m);
    const SearchResult *base = &res[0];   /* Naive is the reference */

    printf("\nPattern: %s\n", pat);
    printf("Occurrences (0-based positions):\n");
    if (base->count == 0)
        printf("  none\n");
    for (int i = 0; i < base->count && i < MAX_POSITIONS_SHOWN; i++)
        printf("%s%d%s", (i % 10 == 0) ? "  " : "", base->positions[i],
               (i % 10 == 9 || i == base->count - 1 || i == MAX_POSITIONS_SHOWN - 1) ? "\n" : ", ");
    if (base->count > MAX_POSITIONS_SHOWN)
        printf("  ... %d more\n", base->count - MAX_POSITIONS_SHOWN);
    printf("Total occurrences: %d\n", base->count);

    printf("\n%-12s %8s %14s %12s  %s\n", "Algorithm", "Matches", "Comparisons", "Time (ms)", "Same as Naive");
    printf("------------------------------------------------------------------\n");
    for (int a = 0; a < SEARCH_ALGO_COUNT; a++) {
        printf("%-12s %8d %14lld %12.5f  %s\n", SEARCH_ALGOS[a].name, res[a].count,
               res[a].comparisons, res[a].time_ms,
               search_results_equal(base, &res[a]) ? "yes" : "NO");
        if (SEARCH_ALGOS[a].fn == rabin_karp_search)
            printf("%-12s %8s   hash checks: %lld, spurious hits: %lld\n", "", "",
                   res[a].hash_checks, res[a].spurious_hits);
    }
    printf("Comparisons include preprocessing (KMP prefix function, BM good-suffix table, Z array).\n");

    for (int a = 0; a < SEARCH_ALGO_COUNT; a++)
        search_result_free(&res[a]);
    free(res);
    free(text);
    free(pat);
}

/* ------------------------------------------------------------------ */
/* 2. Sequence alignment                                               */
/* ------------------------------------------------------------------ */

static void run_alignment(const char *title, AlignFn fn, const char *a, const char *b,
                          ScoreScheme s) {
    AlignResult r = align_run(fn, a, b, s);

    printf("\n=== %s ===\n", title);
    printf("Score: %d | DP cells: %lld | Time: %.5f ms\n", r.score, r.cells, r.time_ms);

    if (r.rows - 1 <= MATRIX_PRINT_MAX && r.cols - 1 <= MATRIX_PRINT_MAX) {
        printf("\nDP matrix (rows: reference, columns: query):\n");
        align_print_matrix(&r, a, b, stdout);
    }

    if (r.length == 0) {
        printf("\nNo alignment with a positive score.\n");
    } else {
        printf("\nAlignment: reference %d..%d, query %d..%d\n", r.a_start, r.a_end - 1,
               r.b_start, r.b_end - 1);
        printf("('|' match, '*' substitution, '-' gap)\n\n");
        align_print(&r, stdout, MAX_ALIGN_BLOCKS);
        MutationReport mr = mutations_detect(&r);
        mutations_print(&mr, stdout, MAX_MUTATIONS_SHOWN);
        mutations_free(&mr);
    }
    align_result_free(&r);
}

static void menu_align(void) {
    char *a = input_sequence("Reference", NULL);
    if (!a)
        return;
    char *b = input_sequence("Query", a);
    if (!b) {
        free(a);
        return;
    }

    long long cells = (long long)(strlen(a) + 1) * (long long)(strlen(b) + 1);
    if (cells > MAX_ALIGN_CELLS) {
        printf("Inputs too large: the DP matrix would need %lld cells (limit %lld).\n",
               cells, MAX_ALIGN_CELLS);
        free(a);
        free(b);
        return;
    }

    printf("\nAlgorithm:\n  1. Needleman-Wunsch (global)\n  2. Smith-Waterman (local)\n  3. Both\n");
    int c = prompt_int("Choice [3]: ", 3);

    ScoreScheme s = score_default();
    if (!prompt_yes("Use default scores (match +1, mismatch -1, gap -2)? [Y/n]: ", 1)) {
        s.match = prompt_int("  Match [1]: ", DEFAULT_MATCH);
        s.mismatch = prompt_int("  Mismatch [-1]: ", DEFAULT_MISMATCH);
        s.gap = prompt_int("  Gap [-2]: ", DEFAULT_GAP);
    }

    if (c == 1 || c == 3)
        run_alignment("Needleman-Wunsch (global)", needleman_wunsch, a, b, s);
    if (c == 2 || c == 3)
        run_alignment("Smith-Waterman (local)", smith_waterman, a, b, s);
    if (c < 1 || c > 3)
        printf("Invalid choice.\n");

    free(a);
    free(b);
}

/* ------------------------------------------------------------------ */
/* 3. Genome reconstruction                                            */
/* ------------------------------------------------------------------ */

static void shuffle(char **list, int k, unsigned long long *st) {
    for (int i = k - 1; i > 0; i--) {
        int j = (int)(rng_next(st) % (unsigned long long)(i + 1));
        char *t = list[i];
        list[i] = list[j];
        list[j] = t;
    }
}

static void menu_reconstruct(void) {
    printf("\nFragment input:\n");
    printf("  1. Load from file (one fragment per line, or FASTA)\n");
    printf("  2. Type fragments\n");
    printf("  3. Simulate: cut a random genome into overlapping reads\n");
    int c = prompt_int("Choice: ", 3);

    char **frags = NULL;
    int k = 0, cap = 0;
    char *original = NULL;

    if (c == 1) {
        char *path = prompt_path("  File path: ");
        frags = seq_read_fragments(path, &k);
        free(path);
    } else if (c == 2) {
        printf("  Enter fragments, one per line. Empty line to finish.\n");
        for (;;) {
            char *line = prompt_line("  > ");
            seq_normalize(line);
            if (!*line) {
                free(line);
                break;
            }
            if (!seq_is_valid(line)) {
                printf("  Skipped: use only A, C, G, T.\n");
                free(line);
                continue;
            }
            seq_list_push(&frags, &k, &cap, line);
        }
    } else if (c == 3) {
        int glen = prompt_int("  Genome length [60]: ", 60);
        int rlen = prompt_int("  Read length [12]: ", 12);
        int step = prompt_int("  Step between reads [5]: ", 5);
        if (glen < 1 || rlen < 1 || step < 1) {
            printf("  All values must be at least 1.\n");
            return;
        }
        if (step >= rlen)
            printf("  Note: step >= read length, so reads will not overlap.\n");
        unsigned long long st = prompt_seed();
        original = seq_random(glen, &st);
        frags = seq_fragment(original, rlen, step, &k);
        shuffle(frags, k, &st);     /* the reads' order must not give the answer away */
        print_truncated("  Original genome: ", original);
    } else {
        printf("  Invalid choice.\n");
        return;
    }

    if (!frags || k == 0) {
        printf("No fragments to reconstruct.\n");
        seq_free_list(frags, k);
        free(original);
        return;
    }

    printf("\nFragments: %d\n", k);
    int k2 = fragments_remove_contained(frags, k);
    if (k2 < k)
        printf("Removed %d duplicate/contained fragment(s); %d remain.\n", k - k2, k2);
    k = k2;
    if (k <= MAX_LIST_SHOWN)
        for (int i = 0; i < k; i++)
            printf("  [%2d] %s\n", i, frags[i]);

    long long ov_ops = 0;
    double t0 = timer_now_ms();
    int **ov = overlap_matrix(frags, k, &ov_ops);
    double t_ov = timer_now_ms() - t0;
    printf("\nOverlap matrix (via KMP prefix function): %lld comparisons, %.4f ms\n", ov_ops, t_ov);
    if (k <= MAX_LIST_SHOWN) {
        printf("      ");
        for (int j = 0; j < k; j++)
            printf("%4d", j);
        printf("\n");
        for (int i = 0; i < k; i++) {
            printf("  %2d  ", i);
            for (int j = 0; j < k; j++)
                printf("%4d", ov[i][j]);
            printf("\n");
        }
    }

    static const struct { const char *name; ReconFn fn; int max_k; } algos[] = {
        { "Brute force", scs_brute_force, BRUTE_FORCE_MAX_K },
        { "DP bitmask",  scs_dp_bitmask,  DP_BITMASK_MAX_K },
        { "Greedy",      scs_greedy,      1000000 },
    };
    ReconResult res[3];
    int have[3] = { 0, 0, 0 };

    for (int a = 0; a < 3; a++) {
        if (k <= algos[a].max_k) {
            res[a] = recon_run(algos[a].fn, frags, k, ov);
            have[a] = 1;
        }
    }

    printf("\n%-12s %8s %14s %12s\n", "Algorithm", "Length", "Operations", "Time (ms)");
    printf("--------------------------------------------------\n");
    for (int a = 0; a < 3; a++) {
        if (have[a])
            printf("%-12s %8d %14lld %12.5f\n", algos[a].name, res[a].length,
                   res[a].operations, res[a].time_ms);
        else
            printf("%-12s %8s   skipped (k = %d > %d)\n", algos[a].name, "-", k, algos[a].max_k);
    }

    printf("\n");
    for (int a = 0; a < 3; a++) {
        if (!have[a])
            continue;
        char label[32];
        snprintf(label, sizeof label, "%-12s ", algos[a].name);
        print_truncated(label, res[a].superstring);
    }

    int optimal = have[1] ? res[1].length : (have[0] ? res[0].length : -1);
    if (optimal > 0)
        printf("\nGreedy / optimal length: %d / %d = %.3f\n", res[2].length, optimal,
               (double)res[2].length / optimal);

    if (original) {
        const char *best = have[1] ? res[1].superstring : res[2].superstring;
        printf("\nOriginal genome length: %d\n", (int)strlen(original));
        if (strcmp(best, original) == 0) {
            printf("Reconstruction exactly matches the original genome.\n");
        } else {
            printf("Reconstruction differs from the original genome.\n");
            printf("(Repeats, or reads that barely overlap, can make the shortest\n"
                   " superstring differ from the true genome.)\n");
        }
    }

    for (int a = 0; a < 3; a++)
        if (have[a])
            recon_result_free(&res[a]);
    overlap_matrix_free(ov, k);
    seq_free_list(frags, k);
    free(original);
}

/* ------------------------------------------------------------------ */
/* 4. Benchmarks                                                       */
/* ------------------------------------------------------------------ */

static void menu_benchmark(void) {
    printf("\nBenchmarks:\n  1. Pattern search\n  2. Alignment\n  3. Reconstruction\n  4. All\n");
    int c = prompt_int("Choice [4]: ", 4);
    char *dir = prompt_path("Output folder [results]: ");
    if (!*dir) {
        free(dir);
        dir = xstrdup("results");
    }

    char path[1024];
    switch (c) {
    case 1:
        snprintf(path, sizeof path, "%s/search.csv", dir);
        bench_search(path);
        break;
    case 2:
        snprintf(path, sizeof path, "%s/align.csv", dir);
        bench_align(path);
        break;
    case 3:
        snprintf(path, sizeof path, "%s/reconstruct.csv", dir);
        bench_reconstruct(path);
        break;
    case 4:
        bench_all(dir);
        break;
    default:
        printf("Invalid choice.\n");
        free(dir);
        return;
    }
    printf("\nPlot the results with:  python3 scripts/plot.py %s\n", dir);
    free(dir);
}

/* ------------------------------------------------------------------ */
/* 5. Utilities                                                        */
/* ------------------------------------------------------------------ */

static void menu_utilities(void) {
    printf("\nUtilities:\n");
    printf("  1. Generate a random genome and save it as FASTA\n");
    printf("  2. Cut a genome file into overlapping reads and save them\n");
    int c = prompt_int("Choice: ", 1);

    if (c == 1) {
        int n = prompt_int("  Length [10000]: ", 10000);
        if (n < 1) {
            printf("  Length must be at least 1.\n");
            return;
        }
        unsigned long long st = prompt_seed();
        char *path = prompt_path("  Save to [data/random_genome.fa]: ");
        if (!*path) {
            free(path);
            path = xstrdup("data/random_genome.fa");
        }
        char *g = seq_random(n, &st);
        if (seq_write_fasta(path, "random_genome", g) == 0)
            printf("  Saved %d bases to %s\n", n, path);
        free(g);
        free(path);
    } else if (c == 2) {
        char *in = prompt_path("  Genome file: ");
        char *g = seq_read_file(in);
        free(in);
        if (!g)
            return;
        int rlen = prompt_int("  Read length [12]: ", 12);
        int step = prompt_int("  Step between reads [5]: ", 5);
        int k = 0;
        char **frags = seq_fragment(g, rlen, step, &k);
        if (!frags) {
            printf("  Read length and step must be at least 1.\n");
            free(g);
            return;
        }
        char *out = prompt_path("  Save to [data/reads.txt]: ");
        if (!*out) {
            free(out);
            out = xstrdup("data/reads.txt");
        }
        if (seq_write_fragments(out, frags, k) == 0)
            printf("  Saved %d reads to %s\n", k, out);
        seq_free_list(frags, k);
        free(out);
        free(g);
    } else {
        printf("  Invalid choice.\n");
    }
}

/* ------------------------------------------------------------------ */
/* Main menu                                                           */
/* ------------------------------------------------------------------ */

int main(void) {
    printf("==============================================\n");
    printf("  DNA-GARS\n");
    printf("  DNA Genome Analysis and Reconstruction System\n");
    printf("==============================================\n");

    for (;;) {
        printf("\nMain menu\n");
        printf("  1. Pattern search        (Naive, KMP, Rabin-Karp, Boyer-Moore, Z)\n");
        printf("  2. Sequence alignment    (Needleman-Wunsch, Smith-Waterman)\n");
        printf("  3. Genome reconstruction (Brute force, DP bitmask, Greedy)\n");
        printf("  4. Benchmarks            (CSV for plotting)\n");
        printf("  5. Utilities             (generate / fragment genomes)\n");
        printf("  0. Exit\n");

        switch (prompt_int("Choice: ", -1)) {
        case 1: menu_search(); break;
        case 2: menu_align(); break;
        case 3: menu_reconstruct(); break;
        case 4: menu_benchmark(); break;
        case 5: menu_utilities(); break;
        case 0: return 0;
        default: printf("Invalid choice.\n");
        }
    }
}
