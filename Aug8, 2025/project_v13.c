#include "function.h"
#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

#define MAX_N 20 // max length of a permutation (increased from 10)

// Global variables - now using pointers instead of fixed arrays
int n;
int *D; // Dynamic array for distances

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

    // Load save distance array
    long long FACT;
    long long size;

    int opt = 2;
    switch (opt)
    {
    case 1: // Lehmer Code ranking
    {
        // int *distance_array = ComputeTDistanceFromIdentity_LehmerCode(n);

        int *distance_array_LehmerCode = load_D_from_file_LehmerCode(n, &size);
        printf("T(n,d) using Lehmer Code Ranking:\n ");
        printf("\n");
        for (int d = 1; d < 8; d++)
        {
            long long result = T_LehmerCode(n, d, distance_array_LehmerCode);
            printf("T(%d,%d) = %lld\n", n, d, result);
        }

        break;
    }
    case 2: // Lexigraphical ranking
    {
        // int *distance_array = ComputeTDistanceFromIdentity_lex(n);

        int *distance_array_Lex = load_D_from_file_lex(n, &size);
        printf("T(n,d) using Lexigraphical ranking:");
        printf("\n");
        for (int d = 1; d < 8; d++)
        {
            long long result = T_lex(n, d, distance_array_Lex);
            printf("T(%d,%d) = %lld\n", n, d, result);
        }

        clock_t end_time = clock();
        double total_program_time = (double)(end_time - start_time) / CLOCKS_PER_SEC;

        printf("\nTotal program execution time: %.3f seconds\n", total_program_time);
        break;
        return 0;
    }
    }
}

// int main()
// {
//     clock_t start_time = clock();

//     printf("Enter the value of n (permutation length): ");
//     scanf("%d", &n);

//     if (n <= 0 || n > MAX_N)
//     {
//         printf("Error: n must be between 1 and %d\n", MAX_N);
//         return 1;
//     }

//     printf("Starting computation for n=%d\n", n);
//     printf("This will process %lld permutations\n", factorial(n));

//     int *pi = (int *)malloc(n * sizeof(int));

//     initialize_identity_permutation(pi, n);

//     // Load save distance array
//     long long size = factorial(n);

//     int *distance_array = load_D_from_file(n, &size);
//     D = distance_array;

//     //============ Case value depends on what you want
//     int x = 2;
//     switch (x)
//     {
//     case 1:
//         printf("-----Level 1 max consecutive value of n = %d---\n", n);
//         printBadTranslocationMaxConsecutiveValue_Level1(n, distance_array);
//         break;
//     case 2:
//         printf("-----Level 2 max consecutive value of n = %d---\n", n);
//         printBadTranslocationMaxConsecutiveValue_Level2(n, distance_array);
//         break;
//     case 3:
//         printf("-----Level 1 max odd cycle of n = %d---\n", n);
//         printBadTranslocationOddCycle_Level1(n, distance_array);
//         break;
//     case 4:
//         printf("-----Level 2 max odd cycle of n = %d---\n", n);
//         printBadTranslocationOddCycle_Level2(n, distance_array);
//         break;
//     case 5:
//         printf("-----Level 1 combined of n = %d---\n", n);
//         printBadTranslocationFromIdentityCombined_Level1(n, distance_array);
//         break;
//     case 6:
//         printf("-----Level 2 combined of n = %d---\n", n);
//         printBadTranslocationFromIdentityCombined_Level2(n, distance_array);
//         // [5 3 2 1 4 0 ] bad
//         break;
//     }

//     clock_t end_time = clock();
//     double total_program_time = (double)(end_time - start_time) / CLOCKS_PER_SEC;

//     printf("\nTotal program execution time: %.3f seconds\n", total_program_time);

//     return 0;
// }