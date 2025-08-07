// Rank1, rank2 from "Ranking and unranking permutations in linear time" by Wendy Myrvold and Frank Ruskey

#include <stdio.h>

void swap(int *a, int *b)
{
    int temp = *a;
    *a = *b;
    *b = temp;
}

// rank1 function: computes the lexicographic rank of a permutation
int rank1(int n, int pi[], int pi_inv[])
{
    if (n == 1)
        return 0;

    int s = pi[n - 1];

    swap(&pi[n - 1], &pi[pi_inv[n - 1]]);
    swap(&pi_inv[s], &pi_inv[n - 1]);

    return s + n * rank1(n - 1, pi, pi_inv);
}

// unrank1: Builds a permutation from a given rank
void unrank1(int n, int r, int pi[])
{
    if (n > 0)
    {
        swap(&pi[n - 1], &pi[r % n]);
        unrank1(n - 1, r / n, pi);
    }
}

// Helper to print a permutation
void print_array(int arr[], int n)
{
    for (int i = 0; i < n; i++)
        printf("%d ", arr[i]);
    printf("\n");
}

// Helper to compute inverse permutation
void compute_inverse(int pi[], int pi_inv[], int n)
{
    for (int i = 0; i < n; i++)
        pi_inv[pi[i]] = i;
}

int main()
{
    // Example permutation of 0..n-1
    // int pi[] = {1, 2, 3, 0}; // permutation π
    // int pi[] = {2, 1, 3, 0}; // permutation π
    int pi[] = {0, 1, 2, 3}; // permutation π

    int pi_inv[10]; // inverse permutation π⁻¹
    int n = 4;

    print_array(pi, n);
    compute_inverse(pi, pi_inv, n);
    print_array(pi_inv, n);
    int result = rank1(n, pi, pi_inv);
    printf("Rank: %d\n", result);

    unrank1(n, result, pi);
    print_array(pi, n);

    return 0;
}