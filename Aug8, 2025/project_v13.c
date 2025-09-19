#include "function.h"
#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

#define MAX_N 20  // max length of a permutation (increased from 10)
#define INF 99999 // Large number represent infinity

// Global variables - now using pointers instead of fixed arrays
int n;
long long FACT;
int *D;        // Dynamic array for distances
int *Q;        // Dynamic array for queue
bool *visited; // Array to track visited permutations

int main()
{
    clock_t start_time = clock();

    printf("Enter the value of n (permutation length): ");
    scanf("%d", &n);

    if (n <= 0 || n > MAX_N)
    {
        printf("Error: n must be between 1 and %d\n", MAX_N);
        return 1;
    }

    printf("Starting computation for n=%d\n", n);
    printf("This will process %lld permutations\n", factorial(n));

    int *pi = (int *)malloc(n * sizeof(int));

    initialize_identity_permutation(pi, n);

    // Load save distance array
    long long size = factorial(n);

    int *distance_array = load_D_from_file(n, &size);
    D = distance_array;

    // printf("-----Level 2 combined of n = %d---\n", n);
    // printBadTranslocationFromIdentityCombined_Level2(n, distance_array);

    printf("-----Level 1 combined of n = %d---\n", n);
    printBadTranslocationFromIdentityCombined_Level1(n, distance_array);

    //------------ End of print bad translocation combined level 1

    clock_t end_time = clock();
    double total_program_time = (double)(end_time - start_time) / CLOCKS_PER_SEC;

    printf("\nTotal program execution time: %.3f seconds\n", total_program_time);

    return 0;
}