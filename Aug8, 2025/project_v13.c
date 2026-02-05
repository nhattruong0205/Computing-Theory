#include "function.h"
#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

#define MAX_N 20 // max length of a permutation (increased from 10)

// Global variables - now using pointers instead of fixed arrays
int n;
int *D;

int main()
{
    int n = 12;                             // Start with small n for testing
    long long max_code_size = 100000000000; // Maximum permutations we want
    long long *selected = (long long *)malloc(max_code_size * sizeof(long long));

    // Allocate space for the permutation array
    // int (*perms)[n] = malloc(max_code_size * sizeof(*perms));

    // Compute the (n,2)-PA
    // int code_size = compute_n2_PA(n, perms, max_code_size, "LehmerAscendingRadix");
    int code_size = compute_n2_PA(n, selected, max_code_size, "LehmerAscendingRadix");
    printf("\n=== Final (n,2)-PA ===\n");
    printf("n = %d, code size = %d\n\n", n, code_size);

    // Print all permutations in the code
    for (int i = 0; i < code_size; i++)
    {
        // printf("Perm %3d: ", i + 1);
        // print_array(perms[i], n);
    }

    // Verify it's a valid (n,2)-PA
    // bool is_valid = verify_n2_PA(n, perms, code_size, "LehmerAscendingRadix");
    bool is_valid = verify_n2_PA(n, "LehmerAscendingRadix");

    printf("\nVerification: %s\n", is_valid ? "VALID (n,2)-PA" : "INVALID");

    // free(perms);
    return 0;
}

// int main()
// {
//     int n = 12;
//     int num_perms = 2;
//     int perms[][12] = {{0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11},
//                        {2, 4, 6, 7, 8, 3, 1, 0, 9, 5, 10, 11},
//                        {11, 10, 9, 8, 7, 6, 5, 4, 3, 2, 1, 0}};
//     // int perms[][4] = {{},
//     //                   {},
//     //                   {}};
//     int pi[] = {0, 1, 2, 3, 4, 5, 6, 7, 9, 10, 11, 8};

//     long long *selected = compute_ranks(n, perms, num_perms);

//     printf("Ranks: ");
//     for (int i = 0; i < num_perms; i++)
//     {
//         printf("%lld ", selected[i]);
//     }
//     printf("\n");

//     bool can_add = can_add_to_code_incremental_d2(n, pi, selected, num_perms, "LehmerAscendingRadix");

//     printf("Can add: %s\n", can_add ? "Yes" : "No");

//     bool is_valid = verify_n2_PA(n, perms, num_perms, "LehmerAscendingRadix");
//     printf("Is valid (n,2)-PA: %s\n", is_valid ? "Yes" : "No");
// }
// int main()
// {
//     int n = 7; // permutation length

//     // Your permutation (0-indexed)
//     int pi[] = {0, 5, 6, 2, 1, 4, 3};

//     int pi[] = {0, 5, 6, 2, 1, 4, 3};

//     // Shift to 1-indexed for cycle computation
//     int pi_shifted[MAX_N];
//     shift_permutation_by_one(pi, pi_shifted, n);

//     // Build breakpoint graph
//     int (*adj)[2] = calloc(2 * n + 1, sizeof *adj);
//     int *deg = calloc(2 * n + 1, sizeof *deg);

//     build_black_edges_from_perm(adj, deg, pi_shifted, n);
//     build_gray_edges_identity(adj, deg, n);

//     // Count odd cycles
//     int odd_cycles = count_odd_cycles(adj, n);
//     printf("Number of odd cycles: %d\n", odd_cycles);

//     // Or print all cycles with details
//     print_cycles(adj, n);

//     free(adj);
//     free(deg);

//     return 0;
// }

// int main()
//  {
//      int n = 7;
//      const char *rank_name = "LehmerAscendingRadix";

//     // Allocate memory
//     if (!allocate_memory(n))
//     {
//         printf("Memory allocation failed!\n");
//         return 1;
//     }

//     // Load or compute distance array
//     long long size_out;
//     int *distance_array = load_D_from_file(n, &size_out, rank_name);

//     if (!distance_array)
//     {
//         printf("Distance array not found. Computing it (this may take a while)...\n");
//         ComputeTDistanceFromIdentity(n, rank_name);
//     }
//     else
//     {
//         // Copy loaded array to global D
//         memcpy(D, distance_array, size_out * sizeof(int));
//         free(distance_array);
//     }

//     // Define your two permutations (0-indexed)
//     // Example: [1, 6, 5, 4, 7, 3, 2] -> [1, 2, 3, 4, 5, 6, 7]
//     int start[MAX_N] = {0, 5, 4, 3, 6, 2, 1};  // [1, 6, 5, 4, 7, 3, 2]
//     int target[MAX_N] = {0, 1, 2, 3, 4, 5, 6}; // [1, 2, 3, 4, 5, 6, 7]

//     printf("Start permutation (0-indexed):  ");
//     print_array(start, n);
//     printf("Target permutation (0-indexed): ");
//     print_array(target, n);

//     // Show the steps
//     show_translocation_steps(n, start, target, rank_name);

//     // Clean up
//     free_memory();

//     return 0;
// }

// Find min distance from 2 permutations
// int main()
// {
//     int n = 7;
//     const char *rank_name = "LehmerAscendingRadix"; // or "Lex", "Lehmer", "SJT", "ReverseColexOrder"

//     // Allocate memory
//     if (!allocate_memory(n))
//     {
//         printf("Memory allocation failed!\n");
//         return 1;
//     }

//     // Load or compute distance array
//     long long size_out;
//     int *distance_array = load_D_from_file(n, &size_out, rank_name);

//     if (!distance_array)
//     {
//         printf("Distance array not found. Computing it...\n");
//         distance_array = ComputeTDistanceFromIdentity(n, rank_name);
//         if (!distance_array)
//         {
//             printf("Failed to compute distance array!\n");
//             free_memory();
//             return 1;
//         }
//     }

//     // Define your two permutations (0-indexed)
//     int pi[MAX_N] = {0, 5, 4, 3, 6, 2, 1};    // [1, 6, 5, 4, 7, 3, 2] in 1-indexed
//     int sigma[MAX_N] = {0, 1, 2, 3, 4, 5, 6}; // [1, 2, 3, 4, 5, 6, 7] in 1-indexed

//     printf("Permutation 1: ");
//     print_array(pi, n);
//     printf("Permutation 2: ");
//     print_array(sigma, n);

//     // Calculate translocation distance
//     int dist = distance_between_2_permutations(n, pi, sigma, distance_array, rank_name);

//     printf("\nTranslocation distance: %d\n", dist);

//     // Clean up
//     free_memory();

//     return 0;
// }

// Calculate neighbor max odd cycles
// int main()
// {
//     // Set the permutation size
//     int n = 7;

//     // Allocate memory for the permutation
//     int *perm = malloc(n * sizeof(int));
//     if (!perm)
//     {
//         printf("Memory allocation failed!\n");
//         return 1;
//     }

//     // Initialize permutation as [1, 2, 3, 4, 5, 6, 7]
//     // Note: Your code uses 0-indexed permutations internally
//     // So we initialize as [1, 2, 3, 4, 5, 6, 7]
//     perm[0] = 1;
//     perm[1] = 6;
//     perm[2] = 5;
//     perm[3] = 4;
//     perm[4] = 7;
//     perm[5] = 3;
//     perm[6] = 2;

//     // Print the input permutation
//     printf("Input permutation (n=%d): ", n);
//     print_array(perm, n);

//     // Compute the maximum odd cycles among all neighbors
//     printf("\nComputing maximum odd cycles among neighbors...\n");

//     clock_t start = clock();
//     int max_odd = computeMostOddCycle(perm, n);
//     clock_t end = clock();

//     double cpu_time = ((double)(end - start)) / CLOCKS_PER_SEC;

//     // Print the result
//     printf("\nMaximum number of odd cycles found: %d\n", max_odd);
//     printf("Computation time: %.3f seconds\n", cpu_time);

//     // Clean up
//     free(perm);

//     return 0;
// }
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

//     // Load save distance array
//     long long FACT;
//     FACT = factorial(n);
//     long long size;

//     // Allocate memory
//     if (!allocate_memory(n))
//     {
//         printf("Memory allocation failed. Exiting.\n");
//         return 1;
//     }

//     int *pi = (int *)malloc(n * sizeof(int));
//     if (!pi)
//     {
//         printf("Failed to allocate pi array\n");
//         free_memory();
//         return 1;
//     }
//     initialize_identity_permutation(pi, n);

//     //========================== OPTION ==============
//     int opt = 2;
//     //========================== OPTION ==============
//     const char *rank_name = NULL;

//     switch (opt)
//     {
//     case 1: // Lehmer Ascending Radix ranking
//     {
//         rank_name = "LehmerAscendingRadix";

//         // int *distance_array_LehmerAscendingRadix = ComputeTDistanceFromIdentity(n, rank_name);
//         int *distance_array_LehmerAscendingRadix = load_D_from_file(n, &size, rank_name);
//         printf("T(n,d) using Lehmer Ascending Radix Ranking:");
//         printf("\n");
//         for (int d = 1; d < 8; d++)
//         {
//             long long result = T_n_d(n, d, distance_array_LehmerAscendingRadix, rank_name);
//             printf("T(%d,%d) = %lld\n", n, d, result);
//         }

//         break;
//     }
//     case 2: // Lexigraphical ranking
//     {
//         rank_name = "Lex";
//         int *distance_array_Lex = ComputeTDistanceFromIdentity(n, rank_name);

//         // int *distance_array_Lex = load_D_from_file(n, &size, rank_name);
//         printf("T(n,d) using Lexigraphical ranking:");
//         printf("\n");
//         for (int d = 1; d < 8; d++)
//         {
//             long long result = T_n_d(n, d, distance_array_Lex, rank_name);
//             printf("T(%d,%d) = %lld\n", n, d, result);
//         }

//         // clock_t end_time = clock();
//         // double total_program_time = (double)(end_time - start_time) / CLOCKS_PER_SEC;

//         // printf("\nTotal program execution time: %.3f seconds\n", total_program_time);
//         break;
//     }
//     case 3: // Lehmer Ranking
//     {
//         rank_name = "Lehmer";

//         // int *distance_array_Lehmer = ComputeTDistanceFromIdentity(n, rank_name);
//         int *distance_array_Lehmer = load_D_from_file(n, &size, rank_name);
//         printf("T(n,d) using Lehmer ranking:");
//         printf("\n");
//         for (int d = 2; d < 8; d++)
//         {
//             long long result = T_n_d(n, d, distance_array_Lehmer, rank_name);
//             printf("T(%d,%d) = %lld\n", n, d, result);
//         }
//         break;
//     }
//     case 4:
//     {
//         rank_name = "SJT";

//         // Allocate arrays for testing
//         int *dir = (int *)malloc(n * sizeof(int));
//         int *test_pi = (int *)malloc(n * sizeof(int));

//         // for (int index = 0; index < FACT; index++)
//         // {
//         //     unrankSJT(index, n, test_pi, dir);
//         //     int computed_rank = rankSJT(n, test_pi);
//         //     printf("Current rank: %d", computed_rank);
//         // }
//         // int *distance_array_SJT = ComputeTDistanceFromIdentity(n, rank_name);
//         int *distance_array_SJT = load_D_from_file(n, &size, rank_name);

//         printf("T(n,d) using SJT ranking:");
//         printf("\n");
//         for (int d = 2; d < 8; ++d)
//         {
//             long long result = T_n_d(n, d, distance_array_SJT, rank_name);
//             printf("T(%d,%d) = %lld\n", n, d, result);
//         }

//         break;
//     }
//     case 5:
//     {
//         rank_name = "ReverseColexOrder";
//         // int *distance_array_ReverseColexOrder = ComputeTDistanceFromIdentity(n, rank_name);
//         int *distance_array_ReverseColexOrder = load_D_from_file(n, &size, rank_name);
//         printf("T(n,d) using Reverse Colex Order ranking:");
//         printf("\n");

//         for (int d = 4; d < 8; ++d)
//         {
//             long long result = T_n_d(n, d, distance_array_ReverseColexOrder, rank_name);
//             printf("T(%d,%d) = %lld\n", n, d, result);
//         }

//         break;
//     }
//     }
//     return 0;
// }

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

// int *distance_array = load_D_from_file(n, &size);
// D = distance_array;

// //============ Case value depends on what you want
// int x = 2;
// switch (x)
// {
// case 1:
//     printf("-----Level 1 max consecutive value of n = %d---\n", n);
//     printBadTranslocationMaxConsecutiveValue_Level1(n, distance_array);
//     break;
// case 2:
//     printf("-----Level 2 max consecutive value of n = %d---\n", n);
//     printBadTranslocationMaxConsecutiveValue_Level2(n, distance_array);
//     break;
// case 3:
//     printf("-----Level 1 max odd cycle of n = %d---\n", n);
//     printBadTranslocationOddCycle_Level1(n, distance_array);
//     break;
// case 4:
//     printf("-----Level 2 max odd cycle of n = %d---\n", n);
//     printBadTranslocationOddCycle_Level2(n, distance_array);
//     break;
// case 5:
//     printf("-----Level 1 combined of n = %d---\n", n);
//     printBadTranslocationFromIdentityCombined_Level1(n, distance_array);
//     break;
// case 6:
//     printf("-----Level 2 combined of n = %d---\n", n);
//     printBadTranslocationFromIdentityCombined_Level2(n, distance_array);
//     // [5 3 2 1 4 0 ] bad
//     break;
// }

//     clock_t end_time = clock();
//     double total_program_time = (double)(end_time - start_time) / CLOCKS_PER_SEC;

//     printf("\nTotal program execution time: %.3f seconds\n", total_program_time);

//     return 0;
// }