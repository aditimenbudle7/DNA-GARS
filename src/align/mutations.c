/*
 * mutations.c — turning an alignment into a mutation report, plus shared
 * alignment helpers (printing, rescoring, freeing).
 *
 * Reading an alignment column by column (reference on top, query below):
 *   A / A  match
 *   A / G  substitution   (SNP)
 *   A / -  deletion       (base missing from the query)
 *   - / G  insertion      (extra base in the query)
 */
#include "align.h"

/* ------------------------------------------------------------------ */
/* Shared helpers                                                      */
/* ------------------------------------------------------------------ */

void align_result_free(AlignResult *r) {
    free(r->aligned_a);
    free(r->aligned_b);
    free(r->matrix);
    r->aligned_a = r->aligned_b = NULL;
    r->matrix = NULL;
}

int align_rescore(const AlignResult *r, ScoreScheme s) {
    int score = 0;
    for (int k = 0; k < r->length; k++) {
        char x = r->aligned_a[k], y = r->aligned_b[k];
        if (x == '-' || y == '-')
            score += s.gap;
        else
            score += (x == y) ? s.match : s.mismatch;
    }
    return score;
}

void align_print(const AlignResult *r, FILE *out, int max_blocks) {
    const int width = 60;
    int pa = r->a_start, pb = r->b_start;   /* original positions at block start */
    int blocks = 0;
    char mid[61];

    for (int k = 0; k < r->length; k += width) {
        if (max_blocks > 0 && blocks == max_blocks) {
            fprintf(out, "... (%d more columns not shown)\n", r->length - k);
            break;
        }
        int w = MIN(width, r->length - k);
        int na = 0, nb = 0;
        for (int t = 0; t < w; t++) {
            char x = r->aligned_a[k + t], y = r->aligned_b[k + t];
            if (x == '-' || y == '-')
                mid[t] = ' ';
            else
                mid[t] = (x == y) ? '|' : '*';
            if (x != '-') na++;
            if (y != '-') nb++;
        }
        mid[w] = '\0';
        fprintf(out, "Ref %7d  %.*s\n", pa, w, r->aligned_a + k);
        fprintf(out, "             %s\n", mid);
        fprintf(out, "Qry %7d  %.*s\n\n", pb, w, r->aligned_b + k);
        pa += na;
        pb += nb;
        blocks++;
    }
}

void align_print_matrix(const AlignResult *r, const char *a, const char *b, FILE *out) {
    int cols = r->cols;
    fprintf(out, "%6s%5s", "", "-");
    for (int j = 1; j < cols; j++)
        fprintf(out, "%5c", b[j - 1]);
    fprintf(out, "\n");
    for (int i = 0; i < r->rows; i++) {
        fprintf(out, "%6c", i == 0 ? '-' : a[i - 1]);
        for (int j = 0; j < cols; j++)
            fprintf(out, "%5d", r->matrix[(size_t)i * (size_t)cols + (size_t)j]);
        fprintf(out, "\n");
    }
}

/* ------------------------------------------------------------------ */
/* Mutation report                                                     */
/* ------------------------------------------------------------------ */

static void add_mutation(MutationReport *m, int *cap, MutationType type, int pos,
                         char ref, char alt) {
    if (m->count == *cap) {
        *cap = *cap ? *cap * 2 : 16;
        m->items = xrealloc(m->items, (size_t)*cap * sizeof(Mutation));
    }
    Mutation mu = { type, pos, ref, alt };
    m->items[m->count++] = mu;
}

MutationReport mutations_detect(const AlignResult *r) {
    MutationReport m;
    memset(&m, 0, sizeof m);
    int cap = 0;
    int ref = r->a_start;   /* current position in the reference */

    for (int k = 0; k < r->length; k++) {
        char x = r->aligned_a[k], y = r->aligned_b[k];
        if (x != '-' && y != '-') {
            if (x == y) {
                m.matches++;
            } else {
                add_mutation(&m, &cap, MUT_SUBSTITUTION, ref, x, y);
                m.substitutions++;
            }
            ref++;
        } else if (x == '-') {
            add_mutation(&m, &cap, MUT_INSERTION, ref, '-', y);
            m.insertions++;
        } else {
            add_mutation(&m, &cap, MUT_DELETION, ref, x, '-');
            m.deletions++;
            ref++;
        }
    }
    m.identity = r->length ? 100.0 * m.matches / r->length : 0.0;
    return m;
}

void mutations_print(const MutationReport *m, FILE *out, int max_items) {
    fprintf(out, "Mutation report (reference positions, 0-based)\n");
    fprintf(out, "  Matches:       %d\n", m->matches);
    fprintf(out, "  Substitutions: %d\n", m->substitutions);
    fprintf(out, "  Insertions:    %d\n", m->insertions);
    fprintf(out, "  Deletions:     %d\n", m->deletions);
    fprintf(out, "  Identity:      %.2f%%\n", m->identity);

    if (m->count == 0) {
        fprintf(out, "  No mutations: the aligned region is identical.\n");
        return;
    }
    int shown = (max_items > 0) ? MIN(max_items, m->count) : m->count;
    for (int i = 0; i < shown; i++) {
        const Mutation *mu = &m->items[i];
        switch (mu->type) {
        case MUT_SUBSTITUTION:
            fprintf(out, "  Substitution at %d: %c -> %c\n", mu->ref_pos, mu->ref_base, mu->alt_base);
            break;
        case MUT_INSERTION:
            fprintf(out, "  Insertion before %d: +%c\n", mu->ref_pos, mu->alt_base);
            break;
        case MUT_DELETION:
            fprintf(out, "  Deletion at %d: -%c\n", mu->ref_pos, mu->ref_base);
            break;
        }
    }
    if (shown < m->count)
        fprintf(out, "  ... %d more not shown\n", m->count - shown);
}

void mutations_free(MutationReport *m) {
    free(m->items);
    m->items = NULL;
    m->count = 0;
}
