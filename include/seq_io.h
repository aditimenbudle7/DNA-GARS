/*
 * seq_io.h — reading, validating, generating and fragmenting DNA sequences.
 */
#ifndef DNA_GARS_SEQ_IO_H
#define DNA_GARS_SEQ_IO_H

#include <stdio.h>

/* ---- Line input ---- */

/* Read one line of any length (newline removed). Returns NULL at EOF. Caller frees. */
char *seq_read_line(FILE *in);

/* ---- Validation ---- */

/* In place: remove whitespace and convert to upper case. */
void seq_normalize(char *s);

/* 1 if s is non-empty and contains only A, C, G, T; else 0. */
int seq_is_valid(const char *s);

/* ---- Files ---- */

/* Read a genome from a plain-text or FASTA file (first record only).
 * Returns a normalized, validated sequence, or NULL with an error printed. */
char *seq_read_file(const char *path);

/* Read fragments: one per line, or one per record if the file is FASTA.
 * Returns an array of *count strings, or NULL with an error printed. */
char **seq_read_fragments(const char *path, int *count);

/* Write a sequence as FASTA (60 bases per line). Returns 0 on success. */
int seq_write_fasta(const char *path, const char *header, const char *seq);

/* Write fragments, one per line. Returns 0 on success. */
int seq_write_fragments(const char *path, char **frags, int count);

/* ---- Random generation (deterministic for a given seed) ---- */

/* xorshift64* step. A state of 0 is replaced by a fixed non-zero constant. */
unsigned long long rng_next(unsigned long long *state);

/* Uniform double in [0, 1). */
double rng_uniform(unsigned long long *state);

/* Random genome of length n. Caller frees. */
char *seq_random(int n, unsigned long long *state);

/* Copy of src with random substitutions, insertions and deletions at the given
 * per-base rate (0..1). Counts are written to subs/ins/dels if not NULL. */
char *seq_mutate(const char *src, double rate, unsigned long long *state,
                 int *subs, int *ins, int *dels);

/* ---- Fragments ---- */

/* Cut genome into overlapping reads of frag_len, starting every `step` bases.
 * A final read is added so the end of the genome is always covered. */
char **seq_fragment(const char *genome, int frag_len, int step, int *count);

/* Append s to a growable string list (takes ownership of s). */
void seq_list_push(char ***list, int *count, int *capacity, char *s);

/* Free a list of strings and the list itself. */
void seq_free_list(char **list, int count);

#endif /* DNA_GARS_SEQ_IO_H */
