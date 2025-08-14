// Updated from project_v5.c where I'll not recompute visiting permutation where I checked if it's exist, not adding it to the queue meaning that its translocation permutation will not be there as well.

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
int *D;                        // Dynamic array for distances
int *Q;                        // Dynamic array for queue
long long front = 0, rear = 0; // Changed to long long for large arrays
bool *visited;                 // Array to track visited permutations

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

// Helper to compute inverse permutation
void compute_inverse(int pi[], int pi_inv[], int n)
{
    for (int i = 0; i < n; i++)
    {
        pi_inv[pi[i]] = i;
    }
}

// Helper function to check the longest increasing subsequence
int longest_increasing_contiguous_subsequence(int *arr, int n)
{
    int max_len = 1, curr_len = 1;

    for (int i = 1; i < n; i++)
    {
        if (arr[i] == arr[i - 1] + 1)
        {
            curr_len++;
            if (curr_len > max_len)
            {
                max_len = curr_len;
            }
        }
        else
        {
            curr_len = 1;
        }
    }
    return max_len;
}

int computeMaxLen(int pi[])
{
    int tmp[MAX_N];
    int tmp_inv[MAX_N];
    int max_len = 0;
    for (int i = 0; i < n; ++i)
        for (int j = i + 1; j < n; ++j)
            for (int k = j; k < n; ++k)
            {
                /* Build translocated permutation */
                int idx = 0;

                /* 1. Prefix: [0..i-1] */
                for (int x = 0; x < i; ++x)
                    tmp[idx++] = pi[x];

                /* 2. Block: [j..k] */
                for (int x = j; x <= k; ++x)
                    tmp[idx++] = pi[x];

                /* 3. Middle: [i..j-1] */
                for (int x = i; x < j; ++x)
                    tmp[idx++] = pi[x];

                /* 4. Suffix: [k+1..n-1] */
                for (int x = k + 1; x < n; ++x)
                    tmp[idx++] = pi[x];

                // printf("Its neighbors");
                // print_array(tmp, n);
                // printf("\n");

                int neighbor_sub_len = longest_increasing_contiguous_subsequence(tmp, n);
                if (neighbor_sub_len > max_len)
                {
                    max_len = neighbor_sub_len;
                }
            }
    return max_len;
}

// ------------------------End of help function-----

//------------- Begin of queue functions----------------
// Function to initialize the queue
void initQueue()
{
    front = 0;
    rear = 0;
}

// Function to check if the queue is full
bool isFull()
{
    return rear >= FACT; // Use FACT instead of MAX_SIZE
}

// Function to check if the queue is empty
bool isEmpty()
{
    return front == rear;
}

// Add a rank (integer) to the queue
void enqueue(int rank)
{
    if (isFull())
    {
        printf("Queue overflow! Cannot add rank %d\n", rank);
        printf("Current queue size: %lld, FACT: %lld\n", rear - front, FACT);
        return;
    }
    Q[rear] = rank;
    rear++;
}

// Remove and return a rank from the queue
int dequeue()
{
    if (isEmpty())
    {
        printf("Queue underflow! Queue is empty\n");
        return -1;
    }
    int rank = Q[front];
    front++;
    return rank;
}

// Get the current size of the queue
long long queueSize()
{
    return rear - front;
}

// Function to print the current queue contents (limited output for large queues)
void printQueue()
{
    if (isEmpty())
    {
        printf("Queue is empty\n");
        return;
    }

    printf("Queue contents (size=%lld): [", queueSize());
    long long print_limit = queueSize() < 20 ? queueSize() : 20;
    for (long long i = front; i < front + print_limit; i++)
    {
        printf("%d", Q[i]);
        if (i < front + print_limit - 1)
        {
            printf(", ");
        }
    }
    if (queueSize() > 20)
    {
        printf("...");
    }
    printf("]\n");
}

//------------- End of queue functions----------------

//------------- Memory Management Functions-------------
// Function to calculate factorial
long long factorial(int n)
{
    long long fact = 1;
    for (int i = 1; i <= n; i++)
    {
        fact *= i;
    }
    return fact;
}

// Memory allocation function
bool allocate_memory(int n)
{
    FACT = factorial(n);

    printf("Allocating memory for n=%d (n! = %lld)...\n", n, FACT);

    // Check if factorial is reasonable (warning for very large values)
    if (FACT > 1000000000LL)
    { // 1 billion
        printf("WARNING: Very large factorial (%lld). This will require ~%.1f GB of RAM.\n",
               FACT, (2.0 * FACT * sizeof(int)) / (1024.0 * 1024.0 * 1024.0));
        printf("Continue? (y/n): ");
        char response;
        scanf(" %c", &response);
        if (response != 'y' && response != 'Y')
        {
            return false;
        }
    }

    // Allocate memory for distance array D
    D = (int *)malloc(FACT * sizeof(int));
    if (!D)
    {
        printf("Failed to allocate memory for D array (%lld integers)\n", FACT);
        return false;
    }

    // Initialize D array to INF
    for (long long i = 0; i < FACT; i++)
    {
        D[i] = INF;
    }

    // Allocate memory for queue Q
    Q = (int *)malloc(FACT * sizeof(int));
    if (!Q)
    {
        printf("Failed to allocate memory for Q array (%lld integers)\n", FACT);
        free(D);
        D = NULL;
        return false;
    }

    // Allocate memory for visited array
    visited = (bool *)calloc(FACT, sizeof(bool));
    if (!visited)
    {
        printf("Failed to allocate memory for visited array (%lld integers)\n", FACT);
        free(D);
        free(Q);
        D = NULL;
        Q = NULL;
        return false;
    }

    printf("Memory allocated successfully.\n");
    printf("D array: %lld integers (%.2f MB)\n", FACT, (FACT * sizeof(int)) / (1024.0 * 1024.0));
    printf("Q array: %lld integers (%.2f MB)\n", FACT, (FACT * sizeof(int)) / (1024.0 * 1024.0));
    printf("Visited array: %lld integers (%.2f MB)\n", FACT, (FACT * sizeof(bool)) / (1024.0 * 1024.0));
    printf("Total memory: %.2f MB\n", (2.0 * FACT * sizeof(int) + FACT * sizeof(bool)) / (1024.0 * 1024.0));

    return true;
}

// Memory cleanup function
void free_memory()
{
    if (D)
    {
        free(D);
        D = NULL;
    }
    if (Q)
    {
        free(Q);
        Q = NULL;
    }
    if (visited)
    {
        free(visited);
        visited = NULL;
    }
    printf("Memory freed.\n");
}

//------------- End Memory Management Functions-------------

// ------------------------Rank 1 (ORIGINAL ALGORITHMS)----------------

void swap(int *a, int *b)
{
    int temp = *a;
    *a = *b;
    *b = temp;
}

// Original recursive rank1 function: computes the lexicographic rank of a permutation
int rank1(int n, int pi[], int pi_inv[])
{
    if (n == 1)
        return 0;

    int s = pi[n - 1];

    swap(&pi[n - 1], &pi[pi_inv[n - 1]]);
    swap(&pi_inv[s], &pi_inv[n - 1]);

    return s + n * rank1(n - 1, pi, pi_inv);
}

int rank_safe(int n, const int src[], int *inv_buf)
{
    int tmp[MAX_N], tmp_inv[MAX_N];
    memcpy(tmp, src, n * sizeof(int));         // work on a copy
    memcpy(tmp_inv, inv_buf, n * sizeof(int)); // inv must start correct
    return rank1(n, tmp, tmp_inv);             // rank1 can now swap freely
}

// Original recursive unrank1: Builds a permutation from a given rank
void unrank1(int n, int r, int pi[])
{
    if (n > 0)
    {
        swap(&pi[n - 1], &pi[r % n]);
        unrank1(n - 1, r / n, pi);
    }
}

//--------------End of Rank1-------------------------

// ------------ Begin of main functions-------------

void initialize_identity_permutation(int *pid, int n)
{
    for (int i = 0; i < n; i++)
    {
        pid[i] = i; // Initialize with the identity permutation (0, 1, 2, ..., n-1)
    }
}

void print_D(int n, long long size)
{
    int pi[MAX_N];
    printf("Translocation distances from identity permutation:\n");
    printf("===============================================\n");
    long long print_limit = size < 50 ? size : 50; // Limit output for large n
    for (long long i = 0; i < print_limit; i++)
    {
        initialize_identity_permutation(pi, n); // reset
        unrank1(n, (int)i, pi);
        printf("Rank %lld: ", i);
        print_array(pi, n);
        if (D[i] == INF)
        {
            printf("distance = INF (unreachable)\n");
        }
        else
        {
            printf("distance = %d\n", D[i]);
        }
    }
    if (size > 50)
    {
        printf("... (%lld more permutations)\n", size - 50);
    }
    printf("===============================================\n");
}

int get_max_distance(long long size)
{
    int max_val = -1;
    for (long long i = 0; i < size; i++)
    {
        if (D[i] != INF && D[i] > max_val)
        {
            max_val = D[i];
        }
    }
    return max_val;
}

int *ComputeTDistanceFromIdentity(int n)
{
    int *pi = (int *)malloc(n * sizeof(int));
    int *pi_inv = (int *)malloc(n * sizeof(int));

    if (!pi || !pi_inv)
    {
        printf("Failed to allocate memory for pi or pi_inv\n");
        if (pi)
            free(pi);
        if (pi_inv)
            free(pi_inv);
        return NULL;
    }

    initialize_identity_permutation(pi, n);
    compute_inverse(pi, pi_inv, n);

    int pid = rank_safe(n, pi, pi_inv);
    D[pid] = 0;
    visited[pid] = true;

    initQueue();
    enqueue(pid);

    long long processed = 0;
    int result[MAX_N];
    while (!isEmpty())
    {
        processed++;
        if (processed % 100000 == 0)
        {
            printf("Processed %lld permutations, queue size: %lld\n", processed, queueSize());
        }

        int current_rank = dequeue();
        // Convert rank back to permutation
        initialize_identity_permutation(result, n);
        unrank1(n, current_rank, result);

        int tmp[MAX_N];
        int tmp_inv[MAX_N];
        for (int i = 0; i < n; ++i)
            for (int j = i + 1; j < n; ++j)
                for (int k = j; k < n; ++k)
                {
                    /* Build translocated permutation */
                    int idx = 0;

                    /* 1. Prefix: [0..i-1] */
                    for (int x = 0; x < i; ++x)
                        tmp[idx++] = result[x];

                    /* 2. Block: [j..k] */
                    for (int x = j; x <= k; ++x)
                        tmp[idx++] = result[x];

                    /* 3. Middle: [i..j-1] */
                    for (int x = i; x < j; ++x)
                        tmp[idx++] = result[x];

                    /* 4. Suffix: [k+1..n-1] */
                    for (int x = k + 1; x < n; ++x)
                        tmp[idx++] = result[x];

                    // Get rank of translocated permutation
                    compute_inverse(tmp, tmp_inv, n);
                    int rank_tmp = rank_safe(n, tmp, tmp_inv);

                    if (rank_tmp < 0 || rank_tmp >= FACT)
                    {
                        printf("Error: Invalid rank %d (FACT=%lld)\n", rank_tmp, FACT);
                        continue;
                    }

                    // Check if the permutation has already been visited
                    if (!visited[rank_tmp])
                    {
                        visited[rank_tmp] = true;
                        if (D[rank_tmp] > D[current_rank] + 1)
                        {
                            D[rank_tmp] = D[current_rank] + 1;
                            enqueue(rank_tmp);
                        }
                    }
                }
    }

    printf("Total processed: %lld permutations\n", processed);
    free(pi);
    free(pi_inv);
    return D;
}

void printBadTranslocationFromIdentity(int n, int *distance_array)
{
    int *pi = (int *)malloc(n * sizeof(int));
    int *pi_inv = (int *)malloc(n * sizeof(int));

    long long size = factorial(n);

    initialize_identity_permutation(pi, n);
    print_array(pi, n);
    compute_inverse(pi, pi_inv, n);
    // print_array(pi_inv, n);

    int pid = rank_safe(n, pi, pi_inv);

    printf("PiD: %d\n", pid);
    int count = 0;
    for (int index = 0; index < size; ++index)
    // for (int index = 3; index < size - 1; ++index)
    {
        int current_distance;
        int current_sub_len;
        int neighbor_distance;
        int neighbor_sub_len;
        int result[MAX_N];
        int longestSubSequence_len = 0;
        int longestSubSequence_ind;

        initialize_identity_permutation(result, n);

        if (distance_array[index] == 0)
            continue;
        // if (index != pid)
        {
            initialize_identity_permutation(pi, n);

            unrank1(n, index, pi);
            // printf("Current permutation");
            // print_array(pi, n);
            // printf("\n");

            current_distance = distance_array[index];
            current_sub_len = longest_increasing_contiguous_subsequence(pi, n);
            int tmp[MAX_N];
            int tmp_inv[MAX_N];
            int max_len = 0, max_rank, max_len_2 = 0;

            for (int i = 0; i < n; ++i)
                for (int j = i + 1; j < n; ++j)
                    for (int k = j; k < n; ++k)
                    {
                        /* Build translocated permutation */
                        int idx = 0;

                        /* 1. Prefix: [0..i-1] */
                        for (int x = 0; x < i; ++x)
                            tmp[idx++] = pi[x];

                        /* 2. Block: [j..k] */
                        for (int x = j; x <= k; ++x)
                            tmp[idx++] = pi[x];

                        /* 3. Middle: [i..j-1] */
                        for (int x = i; x < j; ++x)
                            tmp[idx++] = pi[x];

                        /* 4. Suffix: [k+1..n-1] */
                        for (int x = k + 1; x < n; ++x)
                            tmp[idx++] = pi[x];

                        // printf("Its neighbors");
                        // print_array(tmp, n);
                        // printf("\n");

                        neighbor_sub_len = longest_increasing_contiguous_subsequence(tmp, n);
                        if (neighbor_sub_len > max_len)
                        {
                            max_len = neighbor_sub_len;
                            // printf("Max len : %d\n", max_len);
                            // print_array(tmp, n);
                            compute_inverse(tmp, tmp_inv, n);
                            max_rank = rank_safe(n, tmp, tmp_inv);
                            max_len_2 = computeMaxLen(tmp);
                        }
                        else
                        {
                            if (neighbor_sub_len == max_len)
                            {
                                int neighbor_2_len = computeMaxLen(tmp);

                                if (neighbor_2_len > max_len_2)
                                {
                                    compute_inverse(tmp, tmp_inv, n);
                                    max_rank = rank_safe(n, tmp, tmp_inv);
                                    max_len_2 = neighbor_2_len;
                                }
                            }
                            {
                            }
                        }

                        // printf("Neighbor distance: %d", neighbor_distance);
                        // printf("Current distance: %d\n", current_distance);

                        // printf("Neighbor sublen: %d", neighbor_sub_len);
                        // printf("Current sublen: %d\n", current_sub_len);
                    }
            if (distance_array[max_rank] != (distance_array[index] - 1))
            {
                count++;

                initialize_identity_permutation(pi, n);

                unrank1(n, index, pi);

                // printf("Bad index: %d, %d, %d, %d, %d\n", index, max_len, max_len_2, distance_array[index], distance_array[max_rank]);
                // print_array(pi, n);
                neighbor_distance = distance_array[max_rank];
            }
        }
    }
    printf("Number of bad permutation: %d\n", count);
}

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

    // Allocate memory
    if (!allocate_memory(n))
    {
        printf("Memory allocation failed. Exiting.\n");
        return 1;
    }

    int *pi = (int *)malloc(n * sizeof(int));
    if (!pi)
    {
        printf("Failed to allocate pi array\n");
        free_memory();
        return 1;
    }
    initialize_identity_permutation(pi, n);

    // Run the main algorithm
    int *distance_array = ComputeTDistanceFromIdentity(n);
    if (!distance_array)
    {
        printf("Failed to compute distance array\n");
        free(pi);
        free_memory();
        return 1;
    }

    // Get results
    int max_dist = get_max_distance(FACT);
    printf("Maximum reachable distance = %d\n", max_dist);

    printBadTranslocationFromIdentity(n, distance_array);

    clock_t end_time = clock();
    double total_program_time = (double)(end_time - start_time) / CLOCKS_PER_SEC;

    printf("\nTotal program execution time: %.3f seconds\n", total_program_time);

    // Cleanup
    free(pi);
    free_memory();
    return 0;
}