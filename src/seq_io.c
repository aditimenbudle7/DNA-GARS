/*
 * seq_io.c — sequence input/output, validation, random generation, fragmenting.
 */
#include <ctype.h>

#include "common.h"
#include "seq_io.h"

/* ------------------------------------------------------------------ */
/* Line input                                                          */
/* ------------------------------------------------------------------ */

char *seq_read_line(FILE *in) {
    size_t cap = 256, len = 0;
    char *buf = xmalloc(cap);
    int c;

    while ((c = fgetc(in)) != EOF && c != '\n') {
        if (len + 1 >= cap) {
            cap *= 2;
            buf = xrealloc(buf, cap);
        }
        buf[len++] = (char)c;
    }
    if (c == EOF && len == 0) {
        free(buf);
        return NULL;
    }
    if (len > 0 && buf[len - 1] == '\r') /* Windows line endings */
        len--;
    buf[len] = '\0';
    return buf;
}

/* ------------------------------------------------------------------ */
/* Validation                                                          */
/* ------------------------------------------------------------------ */

void seq_normalize(char *s) {
    char *w = s;
    for (; *s; s++)
        if (!isspace((unsigned char)*s))
            *w++ = (char)toupper((unsigned char)*s);
    *w = '\0';
}

int seq_is_valid(const char *s) {
    if (!s || !*s)
        return 0;
    for (; *s; s++)
        if (*s != 'A' && *s != 'C' && *s != 'G' && *s != 'T')
            return 0;
    return 1;
}

/* ------------------------------------------------------------------ */
/* String lists                                                        */
/* ------------------------------------------------------------------ */

void seq_list_push(char ***list, int *count, int *capacity, char *s) {
    if (*count == *capacity) {
        *capacity = *capacity ? *capacity * 2 : 16;
        *list = xrealloc(*list, (size_t)*capacity * sizeof **list);
    }
    (*list)[(*count)++] = s;
}

void seq_free_list(char **list, int count) {
    if (!list)
        return;
    for (int i = 0; i < count; i++)
        free(list[i]);
    free(list);
}

/* ------------------------------------------------------------------ */
/* Files                                                               */
/* ------------------------------------------------------------------ */

char *seq_read_file(const char *path) {
    FILE *f = fopen(path, "r");
    if (!f) {
        fprintf(stderr, "  Cannot open '%s'\n", path);
        return NULL;
    }

    size_t cap = 1024, len = 0;
    char *seq = xmalloc(cap);
    seq[0] = '\0';
    int headers = 0;
    char *line;

    while ((line = seq_read_line(f)) != NULL) {
        if (line[0] == '>') {           /* FASTA header: keep only the first record */
            free(line);
            if (++headers > 1)
                break;
            continue;
        }
        if (line[0] == ';') {           /* old-style FASTA comment */
            free(line);
            continue;
        }
        seq_normalize(line);
        size_t l = strlen(line);
        while (len + l + 1 > cap) {
            cap *= 2;
            seq = xrealloc(seq, cap);
        }
        memcpy(seq + len, line, l + 1);
        len += l;
        free(line);
    }
    fclose(f);

    if (!seq_is_valid(seq)) {
        fprintf(stderr, "  '%s' must contain a non-empty sequence of A, C, G, T only\n", path);
        free(seq);
        return NULL;
    }
    return seq;
}

char **seq_read_fragments(const char *path, int *count) {
    *count = 0;
    FILE *f = fopen(path, "r");
    if (!f) {
        fprintf(stderr, "  Cannot open '%s'\n", path);
        return NULL;
    }

    char **list = NULL;
    int k = 0, cap = 0;
    char *cur = NULL;       /* current FASTA record being built */
    size_t cur_len = 0;
    int fasta = 0;
    char *line;

    while ((line = seq_read_line(f)) != NULL) {
        if (line[0] == '>') {
            if (cur && cur_len)
                seq_list_push(&list, &k, &cap, cur);
            else
                free(cur);
            cur = xstrdup("");
            cur_len = 0;
            fasta = 1;
            free(line);
            continue;
        }
        seq_normalize(line);
        if (!*line) {
            free(line);
            continue;
        }
        if (fasta) {
            size_t l = strlen(line);
            cur = xrealloc(cur, cur_len + l + 1);
            memcpy(cur + cur_len, line, l + 1);
            cur_len += l;
            free(line);
        } else {
            seq_list_push(&list, &k, &cap, line);
        }
    }
    if (cur) {
        if (cur_len)
            seq_list_push(&list, &k, &cap, cur);
        else
            free(cur);
    }
    fclose(f);

    for (int i = 0; i < k; i++) {
        if (!seq_is_valid(list[i])) {
            fprintf(stderr, "  Fragment %d in '%s' contains characters other than A, C, G, T\n",
                    i + 1, path);
            seq_free_list(list, k);
            return NULL;
        }
    }
    if (k == 0) {
        fprintf(stderr, "  No fragments found in '%s'\n", path);
        free(list);
        return NULL;
    }
    *count = k;
    return list;
}

int seq_write_fasta(const char *path, const char *header, const char *seq) {
    FILE *f = fopen(path, "w");
    if (!f) {
        fprintf(stderr, "  Cannot write '%s'\n", path);
        return -1;
    }
    fprintf(f, ">%s\n", header);
    int n = (int)strlen(seq);
    for (int i = 0; i < n; i += 60)
        fprintf(f, "%.*s\n", MIN(60, n - i), seq + i);
    fclose(f);
    return 0;
}

int seq_write_fragments(const char *path, char **frags, int count) {
    FILE *f = fopen(path, "w");
    if (!f) {
        fprintf(stderr, "  Cannot write '%s'\n", path);
        return -1;
    }
    for (int i = 0; i < count; i++)
        fprintf(f, "%s\n", frags[i]);
    fclose(f);
    return 0;
}

/* ------------------------------------------------------------------ */
/* Random generation                                                   */
/* ------------------------------------------------------------------ */

unsigned long long rng_next(unsigned long long *state) {
    unsigned long long x = *state ? *state : 0x9E3779B97F4A7C15ULL;
    x ^= x >> 12;
    x ^= x << 25;
    x ^= x >> 27;
    *state = x;
    return x * 0x2545F4914F6CDD1DULL;
}

double rng_uniform(unsigned long long *state) {
    return (double)(rng_next(state) >> 11) * (1.0 / 9007199254740992.0); /* / 2^53 */
}

static char random_base(unsigned long long *state) {
    return DNA_ALPHABET[rng_next(state) >> 62];
}

char *seq_random(int n, unsigned long long *state) {
    if (n < 0)
        n = 0;
    char *s = xmalloc((size_t)n + 1);
    for (int i = 0; i < n; i++)
        s[i] = random_base(state);
    s[n] = '\0';
    return s;
}

char *seq_mutate(const char *src, double rate, unsigned long long *state,
                 int *subs, int *ins, int *dels) {
    int n = (int)strlen(src);
    int ns = 0, ni = 0, nd = 0, o = 0;
    char *out = xmalloc((size_t)2 * n + 2);

    for (int i = 0; i < n; i++) {
        if (rng_uniform(state) < rate) {
            switch (rng_next(state) % 3) {
            case 0: {                       /* substitution */
                char c;
                do {
                    c = random_base(state);
                } while (c == src[i]);
                out[o++] = c;
                ns++;
                break;
            }
            case 1:                         /* insertion before src[i] */
                out[o++] = random_base(state);
                out[o++] = src[i];
                ni++;
                break;
            default:                        /* deletion of src[i] */
                nd++;
                break;
            }
        } else {
            out[o++] = src[i];
        }
    }
    if (o == 0 && n > 0) {                  /* never return an empty sequence */
        out[o++] = src[0];
        nd--;
    }
    out[o] = '\0';
    if (subs) *subs = ns;
    if (ins)  *ins = ni;
    if (dels) *dels = nd;
    return out;
}

/* ------------------------------------------------------------------ */
/* Fragmenting                                                         */
/* ------------------------------------------------------------------ */

static char *substring(const char *s, int start, int len) {
    char *d = xmalloc((size_t)len + 1);
    memcpy(d, s + start, (size_t)len);
    d[len] = '\0';
    return d;
}

char **seq_fragment(const char *genome, int frag_len, int step, int *count) {
    int n = (int)strlen(genome);
    char **list = NULL;
    int k = 0, cap = 0;
    *count = 0;

    if (n == 0 || frag_len <= 0 || step <= 0)
        return NULL;

    if (frag_len >= n) {
        seq_list_push(&list, &k, &cap, xstrdup(genome));
    } else {
        int last = -1;
        for (int s = 0; s + frag_len <= n; s += step) {
            seq_list_push(&list, &k, &cap, substring(genome, s, frag_len));
            last = s;
        }
        if (last + frag_len < n) /* make sure the tail is covered */
            seq_list_push(&list, &k, &cap, substring(genome, n - frag_len, frag_len));
    }
    *count = k;
    return list;
}
