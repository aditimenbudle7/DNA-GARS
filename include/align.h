/*
 * align.h — Module 2: pairwise sequence alignment and mutation detection.
 *
 * Convention: `a` is the reference, `b` is the query (sample).
 * Mutations are reported relative to the reference, 0-based.
 */
#ifndef DNA_GARS_ALIGN_H
#define DNA_GARS_ALIGN_H

#include "common.h"

typedef struct {
    int match;      /* score for equal bases (positive) */
    int mismatch;   /* score for a substitution (negative) */
    int gap;        /* score per gap character (negative, linear gap penalty) */
} ScoreScheme;

static inline ScoreScheme score_default(void) {
    ScoreScheme s = { DEFAULT_MATCH, DEFAULT_MISMATCH, DEFAULT_GAP };
    return s;
}

typedef struct {
    char *aligned_a;        /* reference with '-' for gaps */
    char *aligned_b;        /* query with '-' for gaps */
    int length;             /* alignment length (columns) */
    int score;
    int a_start, a_end;     /* aligned region of a: [a_start, a_end) */
    int b_start, b_end;     /* aligned region of b: [b_start, b_end) */
    int *matrix;            /* DP score matrix, rows x cols, row-major */
    int rows, cols;         /* |a|+1, |b|+1 */
    long long cells;        /* DP cells filled = |a| * |b| */
    double time_ms;         /* filled by align_run() */
} AlignResult;

typedef AlignResult (*AlignFn)(const char *a, const char *b, ScoreScheme s);

/* Global alignment.  O(n*m) time and space. */
AlignResult needleman_wunsch(const char *a, const char *b, ScoreScheme s);

/* Local alignment.   O(n*m) time and space. */
AlignResult smith_waterman(const char *a, const char *b, ScoreScheme s);

void align_result_free(AlignResult *r);

/* Recompute the score of an alignment from its aligned strings (used by tests). */
int align_rescore(const AlignResult *r, ScoreScheme s);

/* Print the alignment in blocks of 60 columns, middle line:
 * '|' match, '*' substitution, ' ' gap. Shows at most max_blocks blocks (0 = all). */
void align_print(const AlignResult *r, FILE *out, int max_blocks);

/* Print the full DP matrix (intended for short inputs). */
void align_print_matrix(const AlignResult *r, const char *a, const char *b, FILE *out);

/* ---- Mutation report ---- */

typedef enum { MUT_SUBSTITUTION, MUT_INSERTION, MUT_DELETION } MutationType;

typedef struct {
    MutationType type;
    int ref_pos;        /* substitution/deletion: position in reference;
                           insertion: inserted before this reference position */
    char ref_base;      /* '-' for insertions */
    char alt_base;      /* '-' for deletions */
} Mutation;

typedef struct {
    Mutation *items;
    int count;
    int matches, substitutions, insertions, deletions;
    double identity;    /* matches / alignment length * 100 */
} MutationReport;

MutationReport mutations_detect(const AlignResult *r);
void mutations_print(const MutationReport *m, FILE *out, int max_items);
void mutations_free(MutationReport *m);

#endif /* DNA_GARS_ALIGN_H */
