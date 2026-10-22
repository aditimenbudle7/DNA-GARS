/*
 * rabin_karp.c — Rabin-Karp pattern matching with a rolling hash.
 *
 * Bases are encoded A=0, C=1, G=2, T=3 and each window is read as a base-4
 * number modulo a prime q (RK_MOD). Sliding the window by one position updates
 * the hash in O(1):
 *     hash' = ((hash - code(out) * 4^(m-1)) * 4 + code(in)) mod q
 * Only windows whose hash equals the pattern hash are verified character by
 * character. A verified mismatch is a "spurious hit" (hash collision).
 *
 * With base 4 and q ~ 1e9, patterns up to 14 bases hash without collisions;
 * lower RK_MOD in common.h (e.g. to 101) to demonstrate spurious hits.
 *
 * Time:  O(n + m) expected, O(n * m) worst case (every window collides)
 * Space: O(1)
 */
#include "search.h"

SearchResult rabin_karp_search(const char *text, int n, const char *pat, int m) {
    SearchResult r = search_result_new();
    if (m <= 0 || m > n)
        return r;

    const unsigned long long q = RK_MOD;
    const unsigned long long d = RK_BASE;
    unsigned long long high = 1;   /* d^(m-1) mod q: weight of the outgoing base */
    unsigned long long hp = 0;     /* pattern hash */
    unsigned long long ht = 0;     /* current window hash */

    for (int i = 0; i < m - 1; i++)
        high = (high * d) % q;
    for (int i = 0; i < m; i++) {
        hp = (hp * d + (unsigned long long)dna_code(pat[i])) % q;
        ht = (ht * d + (unsigned long long)dna_code(text[i])) % q;
    }

    for (int s = 0; s <= n - m; s++) {
        r.hash_checks++;
        if (hp == ht) {
            int j = 0;
            while (j < m) {
                r.comparisons++;
                if (text[s + j] != pat[j])
                    break;
                j++;
            }
            if (j == m)
                search_result_add(&r, s);
            else
                r.spurious_hits++;
        }
        if (s < n - m) {
            unsigned long long out = ((unsigned long long)dna_code(text[s]) * high) % q;
            ht = (ht + q - out) % q;
            ht = (ht * d + (unsigned long long)dna_code(text[s + m])) % q;
        }
    }
    return r;
}
