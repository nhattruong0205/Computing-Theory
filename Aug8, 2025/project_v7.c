// Updated from project_v6.c where I added function to calculate cycle and maximize the number of odd cycles

#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

//------------- Begin of helper function-------------
// Helper to print a permutation
void print_array(int arr[], int n)
{
    printf("[");
    for (int i = 0; i < n; i++)
    {
        printf("%d ", arr[i]);
    }
    printf("]");
    printf("\n");
}

// Function to create a breakpoint graph from a given permutation
int *creatingBreakpointGraph(int arr[], int size)
{
    // Allocate memory for result array of double size
    int *result = (int *)malloc(2 * size * sizeof(int));

    // Check if memory allocation was successful
    if (result == NULL)
    {
        printf("Memory allocation failed!\n");
        return NULL;
    }

    // For each element in the original permutation
    for (int i = 0; i < size; i++)
    {
        // Map element arr[i] to positions 2*arr[i]-1 and 2*arr[i]
        result[2 * i] = 2 * arr[i] - 1; // 2i-1
        result[2 * i + 1] = 2 * arr[i]; // 2i
    }

    return result;
}
// ------------------------End of help function-----

int main()
{
    int size = 7; // Size of the input array

    clock_t start_time = clock();

    int input[] = {1, 6, 5, 4, 7, 3, 2}; // Input permutation

    printf("\nInput permutation: ");
    print_array(input, size);

    int *result = creatingBreakpointGraph(input, size);
    if (result != NULL)
    {
        printf("Breakpoint graph representation: ");
        print_array(result, size * 2);
        free(result); // Free the allocated memory
    }

    clock_t end_time = clock();
    double total_program_time = (double)(end_time - start_time) / CLOCKS_PER_SEC;

    printf("\nTotal program execution time: %.3f seconds\n", total_program_time);
}