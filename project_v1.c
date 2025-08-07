#include <stdio.h>
#include <stdlib.h>
#include <string.h>

/**
 * Computes the rank of a permutation in O(n) time using the rankK algorithm.
 *
 * This algorithm is based on a mixed-radix number system and an in-place
 * modification of a copy of the permutation. It is not the standard
 * lexicographical rank, but a different ordering.
 *
 * @param pi The permutation array (1-based values in a 0-based array).
 * @param n The size of the permutation.
 * @return The integer rank of the permutation.
 */
int rankK(int *pi, int n)
{
    int r = 0;

    for (int i = 1; i <= n; i++)
    {
        printf(" %d", pi[i]);
    }
    printf("\n");

    // The main loop iterates backward from n down to 2
    for (int k = n; k >= 2; k--)
    {
        int i = pi[k];

        // This is the core logic: follow the chain of pointers
        // until we find a value less than or equal to k
        while (i > k)
        {
            i = pi[i];
        }

        // This is the key "fixing" step to ensure O(n) amortized complexity
        pi[i] = i;

        // Update the rank using the mixed-radix formula
        r = (k - i) + k * r;
    }

    return r;
}

/**
 * Computes a permutation from its rank in O(n) time using the unrankK algorithm.
 *
 * This is the inverse function of rankK. It constructs a permutation
 * from its integer rank.
 *
 * @param r The integer rank of the permutation.
 * @param n The size of the permutation.
 * @return A new array containing the reconstructed permutation. The caller
 * is responsible for freeing this memory.
 */
int temp[20];
void unrankK(int r, int n, int *pi)
{
    pi[1] = 1;
    int i = 1;

    while (r > 0)
    {
        i++;
        pi[i] = i;
        temp[i] = i - (r % i);
        r = r / i;
    }

    for (int j = i; j >= 2; j--)
    {
        int w = pi[i];
        pi[i] = pi[temp[i]];
        pi[temp[i]] = w;
    }
}

// Main function to demonstrate the inverse property
int main()
{
    int n = 4;

    // Test cases for rankK and unrankK
    int test_perms[5] =
        //{0, 1, 4, 3, 2}; // Rank 12
        //{0, 3, 4, 1, 2}; // Rank 12
        {0, 1, 4, 2, 3}; // Rank 12

    //{0, 2, 3, 4, 1}; // Rank 23
    //{0, 3, 2, 4, 1}; // Rank 22

    for (int i = 1; i <= n; i++)
    {
        printf(" %d", test_perms[i]);
    }
    printf("\n");
    int r = rankK(test_perms, n);
    printf("Rank %d\n", r);

    int out[5];
    unrankK(r, n, out);
    for (int i = 1; i <= n; i++)
    {
        printf(" %d", out[i]);
    }
    printf("\n");

    return 0;
}
