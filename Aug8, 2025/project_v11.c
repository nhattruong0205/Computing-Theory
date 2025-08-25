// Find the lower bound of T(n,d) using the greedy algorithm.
// Using standard ranking
// Updated from project_v5.c where I'll not recompute visiting permutation where I checked if it's exist, not adding it to the queue meaning that its translocation permutation will not be there as well.

#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

#define MAX_N 20  // max length of a permutation (increased from 10)
#define INF 99999 // Large number represent infinity

// Global variables - now using pointers instead of fixed arrays
int n, d;
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
void translocate(const int *src, int *dst, int n, int i, int j, int k)
{
    int idx = 0;

    // 1. Prefix: [0..i-1]
    for (int x = 0; x < i; ++x)
        dst[idx++] = src[x];

    // 2. Block: [j..k]
    for (int x = j; x <= k; ++x)
        dst[idx++] = src[x];

    // 3. Middle: [i..j-1]
    for (int x = i; x < j; ++x)
        dst[idx++] = src[x];

    // 4. Suffix: [k+1..n-1]
    for (int x = k + 1; x < n; ++x)
        dst[idx++] = src[x];
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

// Original recursive rank lex function: computes the lexicographic rank of a permutation
int rank_lex(int pi[], int n)
{
    int rank = 0;
    int fact = 1;
    for (int i = 2; i <= n; i++)
        fact *= i; // fact = n!

    bool used[MAX_N] = {false};

    for (int i = 0; i < n; i++)
    {
        fact /= (n - i); // fact = (n-i-1)!
        int smaller = 0;
        for (int j = 0; j < pi[i]; j++)
        {
            if (!used[j])
                smaller++;
        }
        rank += smaller * fact;
        used[pi[i]] = true;
    }
    return rank;
}

// Build the permutation corresponding to rank r in lexicographic order
void unrank_lex(int n, int r, int pi[])
{
    int fact = 1;
    for (int i = 2; i <= n; i++)
        fact *= i;

    int elems[MAX_N];
    for (int i = 0; i < n; i++)
        elems[i] = i;

    for (int i = 0; i < n; i++)
    {
        fact /= (n - i);
        int idx = r / fact;
        r = r % fact;
        pi[i] = elems[idx];
        // remove elems[idx]
        for (int j = idx; j < n - i - 1; j++)
            elems[j] = elems[j + 1];
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
        unrank_lex(n, (int)i, pi); // directly get permutation of rank i
        printf("Rank %lld: ", i);
        print_array(pi, n);

        if (D[i] == INF)
            printf("distance = INF (unreachable)\n");
        else
            printf("distance = %d\n", D[i]);
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

// Function to compute counts
int *get_distance_counts(long long size, int *max_val_out)
{
    int max_val = -1;

    // Find maximum valid distance
    for (long long i = 0; i < size; i++)
    {
        if (D[i] != INF && D[i] > max_val)
        {
            max_val = D[i];
        }
    }

    if (max_val < 0)
    {
        *max_val_out = -1;
        return NULL; // No valid distances
    }

    // Allocate array for counts
    int *count = (int *)calloc(max_val + 1, sizeof(int));
    if (!count)
    {
        *max_val_out = -1;
        return NULL; // Allocation failed
    }

    // Count frequencies
    for (long long i = 0; i < size; i++)
    {
        if (D[i] != INF)
        {
            count[D[i]]++;
        }
    }

    *max_val_out = max_val;
    return count;
}

// Function to print counts
void print_distance_counts(int *count, int max_val)
{
    if (!count || max_val < 0)
    {
        printf("No distance counts available.\n");
        return;
    }

    printf("Distance counts:\n");
    for (int d = 0; d <= max_val; d++)
    {
        if (count[d] > 0)
        {
            printf("Distance %d → %d times\n", d, count[d]);
        }
    }
}

void save_D_to_file(const char *filename, int *D, long long size)
{
    FILE *fp = fopen(filename, "w");
    if (!fp)
    {
        printf("Error: cannot open %s for writing.\n", filename);
        return;
    }
    for (long long i = 0; i < size; i++)
    {
        fprintf(fp, "%d\n", D[i]);
    }
    fclose(fp);
    printf("Saved D array to %s (%lld elements)\n", filename, size);
}

int *ComputeTDistanceFromIdentity(int n)
{
    int *pi = (int *)malloc(n * sizeof(int));
    if (!pi)
    {
        printf("Failed to allocate memory for pi\n");
        return NULL;
    }

    // Identity permutation
    initialize_identity_permutation(pi, n);

    int pid = rank_lex(pi, n); // rank of identity
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
        unrank_lex(n, current_rank, result);

        int tmp[MAX_N];
        for (int i = 0; i < n; ++i)
            for (int j = i + 1; j < n; ++j)
                for (int k = j; k < n; ++k)
                {
                    /* Build translocated permutation */
                    int idx = 0;

                    // 1. Prefix: [0..i-1]
                    for (int x = 0; x < i; ++x)
                        tmp[idx++] = result[x];

                    // 2. Block: [j..k]
                    for (int x = j; x <= k; ++x)
                        tmp[idx++] = result[x];

                    // 3. Middle: [i..j-1]
                    for (int x = i; x < j; ++x)
                        tmp[idx++] = result[x];

                    // 4. Suffix: [k+1..n-1]
                    for (int x = k + 1; x < n; ++x)
                        tmp[idx++] = result[x];

                    // Get rank of translocated permutation
                    int rank_tmp = rank_lex(tmp, n);

                    if (rank_tmp < 0 || rank_tmp >= FACT)
                    {
                        printf("Error: Invalid rank %d (FACT=%lld)\n", rank_tmp, FACT);
                        continue;
                    }

                    // If not visited yet, update distance and enqueue
                    if (!visited[rank_tmp])
                    {
                        visited[rank_tmp] = true;
                        D[rank_tmp] = D[current_rank] + 1;
                        enqueue(rank_tmp);
                    }
                }
    }

    printf("Total processed: %lld permutations\n", processed);

    char filename[512];
    snprintf(filename, sizeof(filename),
             "/Users/nhattruong/Documents/ComputingTheoryDArraydistances_n%d.txt", n);
    save_D_to_file(filename, D, FACT);

    free(pi);
    return D;
}

// Load D array from the fixed directory path
int *load_D_from_file(int n, long long *size_out)
{
    char filepath[512];
    snprintf(filepath, sizeof(filepath),
             "/Users/nhattruong/Documents/ComputingTheoryDArraydistances_n%d.txt", n);

    FILE *f = fopen(filepath, "r");
    if (!f)
    {
        perror("Failed to open file for reading");
        return NULL;
    }

    // Count number of lines to determine size
    long long count = 0;
    int temp;
    while (fscanf(f, "%d", &temp) == 1)
    {
        count++;
    }
    rewind(f); // go back to beginning

    int *D_loaded = (int *)malloc(count * sizeof(int));
    if (!D_loaded)
    {
        fclose(f);
        perror("Memory allocation failed");
        return NULL;
    }

    for (long long i = 0; i < count; i++)
    {
        if (fscanf(f, "%d", &D_loaded[i]) != 1)
        {
            printf("Error reading element %lld\n", i);
            free(D_loaded);
            fclose(f);
            return NULL;
        }
    }

    fclose(f);
    *size_out = count;
    printf("Loaded D array from %s (%lld elements)\n", filepath, count);
    return D_loaded;
}

// Compute distance between two permutations pi and sigma
int distance_between_2_permutations(int n, int *pi, int *sigma, int *D)
{
    // Compute composition: pi⁻¹ ∘ sigma
    int pi_inv[MAX_N];
    for (int i = 0; i < n; i++)
    {
        pi_inv[pi[i]] = i; // inverse of pi
    }

    int composed[MAX_N];
    for (int i = 0; i < n; i++)
    {
        composed[i] = pi_inv[sigma[i]];
    }

    // Rank the composed permutation in lex order
    int r = rank_lex(composed, n);

    return D[r]; // lookup precomputed distance
}

// Function to compute T(n,d) using greedy algorithm
// T(n,d) = maximum size of a code where all pairs have distance >= d
long long T(int n, int d)
{
    // Create forbidden array to track which permutations are too close to chosen ones
    bool *forbidden = (bool *)calloc(FACT, sizeof(bool));
    if (!forbidden)
    {
        printf("Failed to allocate forbidden array\n");
        return -1;
    }

    long long code_size = 0;
    int pi[MAX_N], sigma[MAX_N];

    // Greedy algorithm: keep selecting permutations until none remain
    for (long long i = 0; i < FACT; i++)
    {
        if (!forbidden[i])
        {
            // Select this permutation as a codeword
            code_size++;

            // Convert rank i to permutation pi
            unrank_lex(n, (int)i, pi);

            // print_array(pi, n);

            // Forbid all permutations within distance < d of pi
            for (long long j = 0; j < FACT; j++)
            {
                if (!forbidden[j])
                {
                    // Convert rank j to permutation sigma
                    unrank_lex(n, (int)j, sigma);

                    // Compute distance between pi and sigma
                    int dist = distance_between_2_permutations(n, pi, sigma, D);

                    if (dist < d)
                    {
                        forbidden[j] = true;
                    }
                }
            }
        }
    }

    free(forbidden);
    return code_size;
}

long long V(int n, int r, int *D)
{
    if (r < 0)
        return 0;

    long long count = 0;
    int tmp[MAX_N];

    for (long long rank = 0; rank < FACT; ++rank)
    {
        if (D[rank] <= r)
        {
            ++count;
            //  initialize_identity_permutation(tmp, n);
            //     unrank1(n, (int)rank, tmp);
            //     printf("Rank %lld: ", rank);
            //     print_array(tmp, n);
            //     printf("distance = %d\n", D[rank]);
        }
    }

    printf("\nV(n=%d, r=%d) = %lld\n", n, r, count);
    return count;
}

// Gilbert–Varshamov lower bound using precomputed ball size
long long GV_lower_bound(int n, int d, int *D)
{
    int r = d - 1;
    long long ball = V(n, r, D);
    if (ball == 0)
        return 0; // safety
    long long bound = FACT / ball;
    printf("GV lower bound: T(%d,%d) >= %lld (since %lld!/%lld)\n", n, d, bound, (long long)n, ball);
    return bound;
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

//     // printf("Enter the value of distance d: ");
//     // scanf("%d", &d);

//     printf("Starting computation for n=%d\n", n);
//     printf("This will process %lld permutations\n", factorial(n));

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

//     int *distance_array = ComputeTDistanceFromIdentity(n);
//     if (!distance_array)
//     {
//         printf("Failed to compute distance array\n");
//         free(pi);
//         free_memory();
//         return 1;
//     }

//     // long long size;
//     // int *distance_array = load_D_from_file(n, &size);
//     // Get results
//     int max_dist = get_max_distance(FACT);

//     for (int d_try = 1; d_try <= 7; ++d_try)
//     {
//         V(n, d_try, distance_array);
//     }

//     clock_t end_time = clock();
//     double total_program_time = (double)(end_time - start_time) / CLOCKS_PER_SEC;

//     printf("\nTotal program execution time: %.3f seconds\n", total_program_time);

//     // Cleanup
//     free(pi);
//     free_memory();
//     return 0;
// }

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

    // printf("Enter the value of distance d: ");
    // scanf("%d", &d);

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

    // int perm1[MAX_N] = {2, 4, 1, 3};
    // int perm2[MAX_N] = {2, 3, 1, 4};
    // int dis = distance_between_2_permutations(n, perm1, perm2, distance_array);
    // printf("Distance between pi and sigma: %d", dis);
    // print_D(n, FACT);

    // Get results
    // int max_dist = get_max_distance(FACT);
    // printf("Maximum reachable distance = %d\n", max_dist);

    for (int d = 4; d < 8; d++)
    {
        GV_lower_bound(n, d, distance_array);
    }

    // int *counts = get_distance_counts(FACT, &max_dist);

    // print_distance_counts(counts, max_dist);

    // long long result = T(n, d);
    // printf("T(%d,%d) = %lld\n", n, d, result);

    clock_t end_time = clock();
    double total_program_time = (double)(end_time - start_time) / CLOCKS_PER_SEC;

    printf("\nTotal program execution time: %.3f seconds\n", total_program_time);

    // Cleanup
    free(pi);
    free_memory();
    return 0;
}