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
int *D;                        // Dynamic array for distances
int *Q;                        // Dynamic array for queue
long long front = 0, rear = 0; // Changed to long long for large arrays
bool *visited;                 // Array to track visited permutations

// ===================Helper functions====================
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
void shift_permutation_by_one(int *src, int *dest, int n)
{
    for (int i = 0; i < n; i++)
    {
        dest[i] = src[i] + 1;
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

long long factorial(int n)
{
    long long fact = 1;
    for (int i = 1; i <= n; i++)
    {
        fact *= i;
    }
    return fact;
}

void initialize_identity_permutation(int *pid, int n)
{
    for (int i = 0; i < n; i++)
    {
        pid[i] = i; // Initialize with the identity permutation (0, 1, 2, ..., n-1)
    }
}

//===============Calculating values (max_len, max_cycle)=============
int longest_increasing_consecutive_values(int *arr, int n)
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

int longest_increasing_subsequence(int *arr, int n)
{
    if (n <= 0)
        return 0;

    int *dp = (int *)malloc(n * sizeof(int));

    // Initialize all lengths as 1
    for (int i = 0; i < n; i++)
        dp[i] = 1;

    // Compute LIS values for all indices
    for (int i = 1; i < n; i++)
    {
        for (int j = 0; j < i; j++)
        {
            if (arr[j] < arr[i] && dp[i] < dp[j] + 1)
                dp[i] = dp[j] + 1;
        }
    }

    // Find maximum value in dp array
    int max_len = dp[0];
    for (int i = 1; i < n; i++)
        if (dp[i] > max_len)
            max_len = dp[i];

    free(dp);
    return max_len;
}

int computeMaxLen(int pi[], int n)
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

                int neighbor_sub_len = longest_increasing_consecutive_values(tmp, n);
                if (neighbor_sub_len > max_len)
                {
                    max_len = neighbor_sub_len;
                }
            }
    return max_len;
}

int computeMaxLen_v2(int pi[], int n)
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

                int neighbor_sub_len = longest_increasing_subsequence(tmp, n);
                if (neighbor_sub_len > max_len)
                {
                    max_len = neighbor_sub_len;
                }
            }
    return max_len;
}

// ================= Read distances array from files ============
// Load D array from the fixed directory path
int *load_D_from_file_LehmerAscendingRadix(int n, long long *size_out)
{
    char filepath[512];
    snprintf(filepath, sizeof(filepath),
             "/Users/nhattruong/Documents/ComputingTheoryDArraydistanceLehmerAscendingRadixRank/distances_n%d.txt", n);

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
    // printf("Loaded D array from %s (%lld elements)\n", filepath, count);
    return D_loaded;
}

// Load D array from the fixed directory path
int *load_D_from_file_lex(int n, long long *size_out)
{
    char filepath[512];
    snprintf(filepath, sizeof(filepath),
             "/Users/nhattruong/Documents/ComputingTheoryDArraydistanceStdRank/distances_n%d.txt", n);

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
    // printf("Loaded D array from %s (%lld elements)\n", filepath, count);
    return D_loaded;
}

// Load D array from the fixed directory path
int *load_D_from_file_Lehmer(int n, long long *size_out)
{
    char filepath[512];
    snprintf(filepath, sizeof(filepath),
             "/Users/nhattruong/Documents/ComputingTheoryDArraydistanceLehmerRank/distances_n%d.txt", n);

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
    // printf("Loaded D array from %s (%lld elements)\n", filepath, count);
    return D_loaded;
}

// -----------------Ranking permutation into index.--------------

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

// Rank 2
// Original recursive rank1 function: computes the lexicographic rank of a permutation
int rank2(int n, int pi[], int pi_inv[])
{
    if (n == 1)
        return 0;

    int s = pi[n - 1];

    swap(&pi[n - 1], &pi[pi_inv[n - 1]]);
    swap(&pi_inv[s], &pi_inv[n - 1]);

    return s * factorial(n - 1) + rank2(n - 1, pi, pi_inv);
}

int rank2_safe(int n, const int src[], int *inv_buf)
{
    int tmp[MAX_N], tmp_inv[MAX_N];
    memcpy(tmp, src, n * sizeof(int));         // work on a copy
    memcpy(tmp_inv, inv_buf, n * sizeof(int)); // inv must start correct
    return rank2(n, tmp, tmp_inv);             // rank1 can now swap freely
}

// Original recursive unrank1: Builds a permutation from a given rank
void unrank2(int n, int r, int pi[])
{
    if (n > 0)
    {
        int s = r / factorial(n - 1);
        swap(&pi[n - 1], &pi[s]);
        unrank2(n - 1, r % factorial(n - 1), pi);
    }
}

//==================== Find odd cycle ====================
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

//----------------Edges
static inline int L(int x) { return 2 * x - 1; } // left endpoint
static inline int R(int x) { return 2 * x; }     // right endpoint

static void add_edge(int (*adj)[2], int *deg, int u, int v, const char *label)
{
    if (deg[u] >= 2 || deg[v] >= 2)
    {
        fprintf(stderr, "Degree overflow when adding %s edge %d--%d\n", label, u, v);
        exit(1);
    }
    adj[u][deg[u]++] = v;
    adj[v][deg[v]++] = u;
    // printf("%s: %d -- %d\n", label, u, v);
}

static void build_black_edges_from_perm(int (*adj)[2], int *deg, int *pi, int n)
{
    for (int i = 0; i < n; i++)
    {
        int a = pi[i];
        int b = pi[(i + 1) % n]; // circular adjacency
        add_edge(adj, deg, R(a), L(b), "black");
    }
}

static void build_gray_edges_identity(int (*adj)[2], int *deg, int n)
{
    for (int i = 1; i <= n; i++)
    {
        int j = (i < n) ? (i + 1) : 1; // circular identity adjacency
        add_edge(adj, deg, R(i), L(j), "gray");
    }
}

static void count_cycles_colored(int (*adj)[2], int n)
{
    int V = 2 * n;
    int *vis = calloc(V + 1, sizeof(int));
    int total = 0, odd = 0, even = 0;

    for (int start = 1; start <= V; start++)
    {
        if (vis[start])
            continue;
        int len = 0;
        int curr = start, prev = 0;

        // walk alternates automatically (each vertex has 1 black + 1 gray)
        do
        {
            vis[curr] = 1;
            len++;
            int nxt = (adj[curr][0] != prev) ? adj[curr][0] : adj[curr][1];
            prev = curr;
            curr = nxt;
        } while (!vis[curr]);

        total++;
        // len counts vertices; edges = len. #black edges per cycle = len/2
        if (((len / 2) % 2) == 0)
            even++;
        else
            odd++;
    }

    free(vis);
    printf("Total cycles: %d\nEven cycles: %d\nOdd cycles: %d\n", total, even, odd);
}

int count_odd_cycles(int (*adj)[2], int n)
{
    int V = 2 * n;
    int *vis = calloc(V + 1, sizeof(int));
    int odd = 0;

    for (int start = 1; start <= V; start++)
    {
        if (vis[start])
            continue;
        int len = 0;
        int curr = start, prev = 0;

        do
        {
            vis[curr] = 1;
            len++;
            int nxt = (adj[curr][0] != prev) ? adj[curr][0] : adj[curr][1];
            prev = curr;
            curr = nxt;
        } while (!vis[curr]);

        if (((len / 2) % 2) != 0) // odd cycle
            odd++;
    }

    free(vis);
    return odd;
}

void print_cycles(int (*adj)[2], int n)
{
    int V = 2 * n;
    int *vis = calloc(V + 1, sizeof(int));
    int cycle_num = 0;

    for (int start = 1; start <= V; start++)
    {
        if (vis[start])
            continue;

        cycle_num++;
        printf("Cycle %d: ", cycle_num);

        int curr = start, prev = 0;
        bool first = true;

        do
        {
            vis[curr] = 1;
            if (!first)
                printf(" -> ");
            printf("%d", curr);
            first = false;

            int nxt = (adj[curr][0] != prev) ? adj[curr][0] : adj[curr][1];
            prev = curr;
            curr = nxt;

        } while (!vis[curr]);

        printf("\n");
    }

    free(vis);
}

void reset_graph(int n, int adj[][2], int *deg)
{
    for (int v = 0; v <= 2 * n; v++)
    {
        deg[v] = 0;
        adj[v][0] = 0;
        adj[v][1] = 0;
    }
}
int computeMostOddCycle(int *perm, int n)
{
    int *tmp = malloc(n * sizeof(int));
    int *best_perm = malloc(n * sizeof(int));
    int best_i = -1, best_j = -1, best_k = -1;
    int max_odd = -1;

    int (*adj)[2] = calloc(2 * n + 1, sizeof *adj);
    int *deg = calloc(2 * n + 1, sizeof *deg);

    for (int i = 0; i < n; ++i)
    {
        for (int j = i + 1; j < n; ++j)
        {
            for (int k = j; k < n; ++k)
            {
                translocate(perm, tmp, n, i, j, k);

                memset(adj, 0, (2 * n + 1) * 2 * sizeof(int));
                memset(deg, 0, (2 * n + 1) * sizeof(int));

                build_black_edges_from_perm(adj, deg, tmp, n);
                build_gray_edges_identity(adj, deg, n);

                int odd = count_odd_cycles(adj, n);

                if (odd > max_odd)
                {
                    max_odd = odd;
                    best_i = i;
                    best_j = j;
                    best_k = k;
                    memcpy(best_perm, tmp, n * sizeof(int));
                }
            }
        }
    }

    return max_odd;
}

//==================== End of find max cycle ====================

// ============== Find best neighbor pairs =============

// Function that finds the single best neighbor and returns both its metrics
BestNeighborMetrics findBestNeighborMetrics(int *perm, int n)
{
    int tmp[MAX_N];
    int tmp_shift[MAX_N];
    BestNeighborMetrics result = {-1, -1, -1, -1, -1};

    // Allocate graph structures
    int (*adj)[2] = calloc(2 * n + 1, sizeof *adj);
    int *deg = calloc(2 * n + 1, sizeof *deg);

    for (int i = 0; i < n; ++i)
    {
        for (int j = i + 1; j < n; ++j)
        {
            for (int k = j; k < n; ++k)
            {
                /* Build translocated permutation */
                int idx = 0;

                /* 1. Prefix: [0..i-1] */
                for (int x = 0; x < i; ++x)
                    tmp[idx++] = perm[x];

                /* 2. Block: [j..k] */
                for (int x = j; x <= k; ++x)
                    tmp[idx++] = perm[x];

                /* 3. Middle: [i..j-1] */
                for (int x = i; x < j; ++x)
                    tmp[idx++] = perm[x];

                /* 4. Suffix: [k+1..n-1] */
                for (int x = k + 1; x < n; ++x)
                    tmp[idx++] = perm[x];

                // Compute length for this neighbor
                int neighbor_len = longest_increasing_consecutive_values(tmp, n);

                // Compute cycles for this neighbor
                shift_permutation_by_one(tmp, tmp_shift, n);
                reset_graph(n, adj, deg);
                build_black_edges_from_perm(adj, deg, tmp_shift, n);
                build_gray_edges_identity(adj, deg, n);
                int neighbor_cycles = count_odd_cycles(adj, n);

                // Check if this neighbor is better using the same tie-breaking logic
                // Priority: length first, then cycles as tie-breaker
                if (neighbor_len > result.max_len ||
                    (neighbor_len == result.max_len && neighbor_cycles > result.max_cycles))
                {
                    result.max_len = neighbor_len;
                    result.max_cycles = neighbor_cycles;
                    result.best_i = i;
                    result.best_j = j;
                    result.best_k = k;
                }
            }
        }
    }

    free(adj);
    free(deg);

    return result;
}

// ================= Queue functions ===========================
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

// =============== Dealing with memory ========================
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

    // printf("Memory allocated successfully.\n");
    // printf("D array: %lld integers (%.2f MB)\n", FACT, (FACT * sizeof(int)) / (1024.0 * 1024.0));
    // printf("Q array: %lld integers (%.2f MB)\n", FACT, (FACT * sizeof(int)) / (1024.0 * 1024.0));
    // printf("Visited array: %lld integers (%.2f MB)\n", FACT, (FACT * sizeof(bool)) / (1024.0 * 1024.0));
    // printf("Total memory: %.2f MB\n", (2.0 * FACT * sizeof(int) + FACT * sizeof(bool)) / (1024.0 * 1024.0));

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

// =============== Saving files function ======================
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

//================== Compute distance array ====================
int *ComputeTDistanceFromIdentity_lex(int n)
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
    // snprintf(filename, sizeof(filename),
    //          "/Users/nhattruong/Documents/ComputingTheoryDArraydistanceStdRank/distances_n%d.txt", n);
    // save_D_to_file(filename, D, FACT);

    free(pi);
    return D;
}

int *ComputeTDistanceFromIdentity_LehmerAscendingRadix(int n)
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

    // char filename[512];
    // snprintf(filename, sizeof(filename),
    //          "/Users/nhattruong/Documents/ComputingTheoryDArraydistanceLehmerAscendingRadixRank/distances_n%d.txt", n);
    // save_D_to_file(filename, D, FACT);
    return D;
}

int *ComputeTDistanceFromIdentity_Lehmer(int n)
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

    int pid = rank2_safe(n, pi, pi_inv);
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
        unrank2(n, current_rank, result);

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
                    int rank_tmp = rank2_safe(n, tmp, tmp_inv);

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

    char filename[512];
    snprintf(filename, sizeof(filename),
             "/Users/nhattruong/Documents/ComputingTheoryDArraydistanceLehmerRank/distances_n%d.txt", n);
    save_D_to_file(filename, D, FACT);
    return D;
}

// ================== Computing PAs =========================
// Compute distance between two permutations pi and sigma
int distance_between_2_permutations_LehmerAscendingRadix(int n, int *pi, int *sigma, int *D)
{
    int pi_inv[MAX_N];
    compute_inverse(pi, pi_inv, n);

    // Compute composition: pi_inv ∘ sigma
    int composed[MAX_N];
    for (int i = 0; i < n; i++)
    {
        composed[i] = pi_inv[sigma[i]];
    }

    // Rank the composed permutation
    int composed_inv[MAX_N];
    compute_inverse(composed, composed_inv, n);
    int r = rank_safe(n, composed, composed_inv);

    return D[r]; // lookup precomputed distance
}

// Computing T(n,d) PA - an array A of permutation on [1..n] with dt(A) >= d.
long long T_LehmerAscendingRadix(int n, int d, int *D)
{
    // Create forbidden array to track which permutations are too close to chosen ones
    bool *forbidden = (bool *)calloc(factorial(n), sizeof(bool));
    if (!forbidden)
    {
        printf("Failed to allocate forbidden array\n");
        return -1;
    }

    long long code_size = 0;
    int pi[MAX_N], sigma[MAX_N];
    int pi_inv[MAX_N], sigma_inv[MAX_N];

    // Greedy algorithm: keep selecting permutations until none remain
    for (long long i = 0; i < factorial(n); i++)
    {
        if (!forbidden[i])
        {
            // Select this permutation as a codeword
            code_size++;

            // Convert rank i to permutation pi
            initialize_identity_permutation(pi, n);
            unrank1(n, (int)i, pi);

            // print_array(pi, n);

            // Forbid all permutations within distance d-1 of pi
            for (long long j = 0; j < factorial(n); j++)
            {
                if (!forbidden[j])
                {
                    // Convert rank j to permutation sigma
                    initialize_identity_permutation(sigma, n);
                    unrank1(n, (int)j, sigma);

                    // Compute distance between pi and sigma
                    int dist = distance_between_2_permutations_LehmerAscendingRadix(n, pi, sigma, D);

                    // If distance < d, forbid this permutation
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

// Compute distance between two permutations pi and sigma
int distance_between_2_permutations_lex(int n, int *pi, int *sigma, int *D)
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
long long T_lex(int n, int d, int *D)
{
    // Create forbidden array to track which permutations are too close to chosen ones
    bool *forbidden = (bool *)calloc(factorial(n), sizeof(bool));
    if (!forbidden)
    {
        printf("Failed to allocate forbidden array\n");
        return -1;
    }

    long long code_size = 0;
    int pi[MAX_N], sigma[MAX_N];

    // Greedy algorithm: keep selecting permutations until none remain
    for (long long i = 0; i < factorial(n); i++)
    {
        if (!forbidden[i])
        {
            // Select this permutation as a codeword
            code_size++;

            // Convert rank i to permutation pi
            unrank_lex(n, (int)i, pi);

            // print_array(pi, n);

            // Forbid all permutations within distance < d of pi
            for (long long j = 0; j < factorial(n); j++)
            {
                if (!forbidden[j])
                {
                    // Convert rank j to permutation sigma
                    unrank_lex(n, (int)j, sigma);

                    // Compute distance between pi and sigma
                    int dist = distance_between_2_permutations_lex(n, pi, sigma, D);

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

// Compute distance between two permutations pi and sigma
int distance_between_2_permutations_Lehmer(int n, int *pi, int *sigma, int *D)
{
    int pi_inv[MAX_N];
    compute_inverse(pi, pi_inv, n);

    // Compute composition: pi_inv ∘ sigma
    int composed[MAX_N];
    for (int i = 0; i < n; i++)
    {
        composed[i] = pi_inv[sigma[i]];
    }

    // Rank the composed permutation
    int composed_inv[MAX_N];
    compute_inverse(composed, composed_inv, n);
    int r = rank2_safe(n, composed, composed_inv);

    return D[r]; // lookup precomputed distance
}

// Computing T(n,d) PA - an array A of permutation on [1..n] with dt(A) >= d.
long long T_Lehmer(int n, int d, int *D)
{
    // Create forbidden array to track which permutations are too close to chosen ones
    bool *forbidden = (bool *)calloc(factorial(n), sizeof(bool));
    if (!forbidden)
    {
        printf("Failed to allocate forbidden array\n");
        return -1;
    }

    long long code_size = 0;
    int pi[MAX_N], sigma[MAX_N];
    int pi_inv[MAX_N], sigma_inv[MAX_N];

    // Greedy algorithm: keep selecting permutations until none remain
    for (long long i = 0; i < factorial(n); i++)
    {
        if (!forbidden[i])
        {
            // Select this permutation as a codeword
            code_size++;

            // Convert rank i to permutation pi
            initialize_identity_permutation(pi, n);
            unrank2(n, (int)i, pi);

            // print_array(pi, n);

            // Forbid all permutations within distance d-1 of pi
            for (long long j = 0; j < factorial(n); j++)
            {
                if (!forbidden[j])
                {
                    // Convert rank j to permutation sigma
                    initialize_identity_permutation(sigma, n);
                    unrank2(n, (int)j, sigma);

                    // Compute distance between pi and sigma
                    int dist = distance_between_2_permutations_Lehmer(n, pi, sigma, D);

                    // If distance < d, forbid this permutation
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

// ============== Print bad permutation ==================
void printBadTranslocationMaxConsecutiveValue_Level1(int n, int *distance_array)
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
    long long last_progress = -1;

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

            long long progress = (index * 100) / size;
            if (progress != last_progress)
            {
                printf("\rProgress: %lld%% (%d/%lld)", progress, index, size);
                fflush(stdout);
                last_progress = progress;
            }

            initialize_identity_permutation(pi, n);

            unrank1(n, index, pi);
            // printf("Current permutation");
            // print_array(pi, n);
            // printf("\n");

            current_distance = distance_array[index];
            current_sub_len = longest_increasing_consecutive_values(pi, n);
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

                        neighbor_sub_len = longest_increasing_consecutive_values(tmp, n);
                        if (neighbor_sub_len > max_len)
                        {
                            max_len = neighbor_sub_len;
                            // printf("Max len : %d\n", max_len);
                            // print_array(tmp, n);
                            compute_inverse(tmp, tmp_inv, n);
                            max_rank = rank_safe(n, tmp, tmp_inv);
                        }
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
    printf("\nNumber of bad permutation: %d\n", count);
}

void printBadTranslocationMaxConsecutiveValue_Level2(int n, int *distance_array)
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
    long long last_progress = -1;

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

            long long progress = (index * 100) / size;
            if (progress != last_progress)
            {
                printf("\rProgress: %lld%% (%d/%lld)", progress, index, size);
                fflush(stdout);
                last_progress = progress;
            }

            initialize_identity_permutation(pi, n);

            unrank1(n, index, pi);
            // printf("Current permutation");
            // print_array(pi, n);
            // printf("\n");

            current_distance = distance_array[index];
            current_sub_len = longest_increasing_consecutive_values(pi, n);
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

                        neighbor_sub_len = longest_increasing_consecutive_values(tmp, n);
                        if (neighbor_sub_len > max_len)
                        {
                            max_len = neighbor_sub_len;
                            // printf("Max len : %d\n", max_len);
                            // print_array(tmp, n);
                            compute_inverse(tmp, tmp_inv, n);
                            max_rank = rank_safe(n, tmp, tmp_inv);
                            max_len_2 = computeMaxLen(tmp, n);
                        }
                        else
                        {
                            if (neighbor_sub_len == max_len)
                            {
                                int neighbor_2_len = computeMaxLen(tmp, n);

                                if (neighbor_2_len > max_len_2)
                                {
                                    compute_inverse(tmp, tmp_inv, n);
                                    max_rank = rank_safe(n, tmp, tmp_inv);
                                    max_len_2 = neighbor_2_len;
                                }
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
    printf("\nNumber of bad permutation: %d\n", count);
}

void printBadTranslocationOddCycle_Level1(int n, int *distance_array)
{
    int *pi = (int *)malloc(n * sizeof(int));
    int *pi_inv = (int *)malloc(n * sizeof(int));

    long long size = factorial(n);

    initialize_identity_permutation(pi, n);
    // print_array(pi, n);
    compute_inverse(pi, pi_inv, n);
    // print_array(pi_inv, n);

    int pid = rank_safe(n, pi, pi_inv);

    printf("PiD: %d\n", pid);
    int count = 0;
    long long last_progress = -1;

    for (int index = 0; index < size; ++index)
    // for (int index = 3; index < size - 1; ++index)
    {
        int current_distance;
        int current_maxCycle;
        int neighbor_distance;
        int neighbor_maxCycle;
        int result[MAX_N];
        int pi_shifted[MAX_N];
        int longestOddCycle_len = 0;
        int longestOddCycle_ind;

        long long progress = (index * 100) / size;
        if (progress != last_progress)
        {
            printf("\rProgress: %lld%% (%d/%lld)", progress, index, size);
            fflush(stdout);
            last_progress = progress;
        }

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
            // printf("Current distance: %d", current_distance);
            //------------ Added
            int (*adj)[2] = calloc(2 * n + 1, sizeof *adj);
            int *deg = calloc(2 * n + 1, sizeof *deg);

            shift_permutation_by_one(pi, pi_shifted, n);
            // print_array(pi_shifted, n);
            //  printf("Building edges...\n");
            build_black_edges_from_perm(adj, deg, pi_shifted, n); // e.g., 2--11, 12--9, ...
            build_gray_edges_identity(adj, deg, n);               // e.g., 2--3, 4--5, ..., 14--1
            current_maxCycle = count_odd_cycles(adj, n);
            // printf("\n");

            // printf("Current max cycle\n: %d", current_maxCycle);
            // printf("\n");

            //------------ Added

            int tmp[MAX_N];
            int tmp_shift[MAX_N];
            int tmp_inv[MAX_N];
            int max_cycle = 0, max_rank, max_cycle_2 = 0;

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

                        shift_permutation_by_one(tmp, tmp_shift, n);

                        // printf("Building edges...\n");
                        reset_graph(n, adj, deg);
                        build_black_edges_from_perm(adj, deg, tmp_shift, n); // e.g., 2--11, 12--9, ...
                        build_gray_edges_identity(adj, deg, n);              // e.g., 2--3, 4--5, ..., 14--1
                        neighbor_maxCycle = count_odd_cycles(adj, n);

                        if (neighbor_maxCycle > max_cycle)
                        {
                            max_cycle = neighbor_maxCycle;
                            // printf("Max len : %d\n", max_cycle);
                            //  print_array(tmp, n);
                            compute_inverse(tmp, tmp_inv, n);
                            max_rank = rank_safe(n, tmp, tmp_inv);
                        }
                    }
            if (distance_array[max_rank] != (distance_array[index] - 1))
            {
                count++;

                initialize_identity_permutation(pi, n);

                unrank1(n, index, pi);

                // printf("Bad index: %d, %d, %d, %d, %d\n", index, max_cycle, max_cycle_2, distance_array[index], distance_array[max_rank]);
                // print_array(pi, n);
                neighbor_distance = distance_array[max_rank];
            }
        }
    }
    printf("Number of bad permutation: %d\n", count);
}

void printBadTranslocationOddCycle_Level2(int n, int *distance_array)
{
    int *pi = (int *)malloc(n * sizeof(int));
    int *pi_inv = (int *)malloc(n * sizeof(int));

    long long size = factorial(n);

    initialize_identity_permutation(pi, n);
    // print_array(pi, n);
    compute_inverse(pi, pi_inv, n);
    // print_array(pi_inv, n);

    int pid = rank_safe(n, pi, pi_inv);

    printf("PiD: %d\n", pid);
    int count = 0;
    long long last_progress = -1;

    for (int index = 0; index < size; ++index)
    // for (int index = 3; index < size - 1; ++index)
    {
        int current_distance;
        int current_maxCycle;
        int neighbor_distance;
        int neighbor_maxCycle;
        int result[MAX_N];
        int pi_shifted[MAX_N];

        long long progress = (index * 100) / size;
        if (progress != last_progress)
        {
            printf("\rProgress: %lld%% (%d/%lld)", progress, index, size);
            fflush(stdout);
            last_progress = progress;
        }

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
            // printf("Current distance: %d", current_distance);
            //------------ Added
            int (*adj)[2] = calloc(2 * n + 1, sizeof *adj);
            int *deg = calloc(2 * n + 1, sizeof *deg);

            shift_permutation_by_one(pi, pi_shifted, n);
            // print_array(pi_shifted, n);
            //  printf("Building edges...\n");
            build_black_edges_from_perm(adj, deg, pi_shifted, n); // e.g., 2--11, 12--9, ...
            build_gray_edges_identity(adj, deg, n);               // e.g., 2--3, 4--5, ..., 14--1
            current_maxCycle = count_odd_cycles(adj, n);
            // printf("Current max cycle\n: %d", current_maxCycle);
            //------------ Added

            int tmp[MAX_N];
            int tmp_shift[MAX_N];
            int tmp_inv[MAX_N];
            int max_cycle = 0, max_rank, max_cycle_2 = 0;

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

                        shift_permutation_by_one(tmp, tmp_shift, n);

                        // printf("Building edges...\n");
                        reset_graph(n, adj, deg);
                        build_black_edges_from_perm(adj, deg, tmp_shift, n); // e.g., 2--11, 12--9, ...
                        build_gray_edges_identity(adj, deg, n);              // e.g., 2--3, 4--5, ..., 14--1
                        neighbor_maxCycle = count_odd_cycles(adj, n);

                        if (neighbor_maxCycle > max_cycle)
                        {
                            max_cycle = neighbor_maxCycle;
                            // printf("Max len : %d\n", max_len);
                            // print_array(tmp, n);
                            compute_inverse(tmp, tmp_inv, n);
                            max_rank = rank_safe(n, tmp, tmp_inv);

                            shift_permutation_by_one(tmp, tmp_shift, n);

                            reset_graph(n, adj, deg);

                            // printf("Building edges...\n");
                            build_black_edges_from_perm(adj, deg, tmp_shift, n); // e.g., 2--11, 12--9, ...
                            build_gray_edges_identity(adj, deg, n);              // e.g., 2--3, 4--5, ..., 14--1

                            max_cycle_2 = computeMostOddCycle(tmp_shift, n);
                        }
                        else
                        {
                            if (neighbor_maxCycle == max_cycle)
                            {
                                int neighbor_2_maxCycle = computeMostOddCycle(tmp_shift, n);

                                if (neighbor_2_maxCycle > max_cycle_2)
                                {
                                    compute_inverse(tmp, tmp_inv, n);
                                    max_rank = rank_safe(n, tmp, tmp_inv);
                                    max_cycle_2 = neighbor_2_maxCycle;
                                }
                            }
                        }
                    }
            if (distance_array[max_rank] != (distance_array[index] - 1))
            {
                count++;

                initialize_identity_permutation(pi, n);

                unrank1(n, index, pi);

                // printf("Bad index: %d, %d, %d, %d, %d\n", index, max_cycle, max_cycle_2, distance_array[index], distance_array[max_rank]);
                // print_array(pi, n);
                neighbor_distance = distance_array[max_rank];
            }
        }
    }
    printf("\n");
    printf("Odd cycles level-2 with n = %d\n", n);
    printf("Number of bad permutation: %d\n", count);
}

void printBadTranslocationFromIdentityCombined_Level1(int n, int *distance_array)
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
    long long last_progress = -1;

    for (int index = 0; index < size; ++index)
    // for (int index = 3; index < size - 1; ++index)
    {
        int current_distance;
        int current_sub_len;
        int current_maxCycle;
        int neighbor_maxCycle;
        int neighbor_distance;
        int neighbor_sub_len;
        int result[MAX_N];
        int pi_shifted[MAX_N];

        initialize_identity_permutation(result, n);

        if (distance_array[index] == 0)
            continue;
        // if (index != pid)
        {

            long long progress = (index * 100) / size;
            if (progress != last_progress)
            {
                printf("\rProgress: %lld%% (%d/%lld)", progress, index, size);
                fflush(stdout);
                last_progress = progress;
            }

            initialize_identity_permutation(pi, n);

            unrank1(n, index, pi);
            // printf("Current permutation");
            // print_array(pi, n);
            // printf("\n");

            current_distance = distance_array[index];
            current_sub_len = longest_increasing_consecutive_values(pi, n);

            //------------ Added
            int (*adj)[2] = calloc(2 * n + 1, sizeof *adj);
            int *deg = calloc(2 * n + 1, sizeof *deg);

            // shift_permutation_by_one(pi, pi_shifted, n);
            // // print_array(pi_shifted, n);
            // //  printf("Building edges...\n");
            // build_black_edges_from_perm(adj, deg, pi_shifted, n); // e.g., 2--11, 12--9, ...
            // build_gray_edges_identity(adj, deg, n);               // e.g., 2--3, 4--5, ..., 14--1
            // current_maxCycle = count_odd_cycles(adj, n);
            // // printf("Current max cycle\n: %d", current_maxCycle);
            //------------ Added

            int tmp[MAX_N];
            int tmp_inv[MAX_N];
            int max_len = -1, max_rank = -1, max_cycle = -1;
            int tmp_shift[MAX_N];

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

                        neighbor_sub_len = longest_increasing_consecutive_values(tmp, n);

                        if (neighbor_sub_len > max_len)
                        {
                            max_len = neighbor_sub_len;

                            compute_inverse(tmp, tmp_inv, n);
                            max_rank = rank_safe(n, tmp, tmp_inv);
                        }
                        else if (neighbor_sub_len == max_len)
                        {
                            // ------------Calculate current max_cycle ----------
                            int temp_perm[MAX_N];
                            int temp_perm_shift[MAX_N];

                            initialize_identity_permutation(temp_perm, n);
                            unrank1(n, max_rank, temp_perm);
                            // compute cycle now so it's in sync
                            reset_graph(n, adj, deg);

                            shift_permutation_by_one(temp_perm, temp_perm_shift, n);
                            reset_graph(n, adj, deg);
                            build_black_edges_from_perm(adj, deg, temp_perm_shift, n);
                            build_gray_edges_identity(adj, deg, n);
                            max_cycle = count_odd_cycles(adj, n);

                            // ------------End of Calculate current max_cycle ----------

                            // ------------Calculate current neighbor_maxCycle ----------

                            shift_permutation_by_one(tmp, tmp_shift, n);

                            reset_graph(n, adj, deg);
                            build_black_edges_from_perm(adj, deg, tmp_shift, n); // e.g., 2--11, 12--9, ...
                            build_gray_edges_identity(adj, deg, n);              // e.g., 2--3, 4--5, ..., 14--1
                            neighbor_maxCycle = count_odd_cycles(adj, n);

                            // ------------End of calculate current neighbor_maxCycle ----------

                            if (neighbor_maxCycle > max_cycle)
                            {
                                max_cycle = neighbor_maxCycle;

                                compute_inverse(tmp, tmp_inv, n);
                                max_rank = rank_safe(n, tmp, tmp_inv);
                            }
                        }
                    }
            if (distance_array[max_rank] != (distance_array[index] - 1))
            {
                count++;

                initialize_identity_permutation(pi, n);

                unrank1(n, index, pi);

                printf("Bad index: %d, %d, %d\n", max_len, max_cycle, index);

                print_array(pi, n);

                neighbor_distance = distance_array[max_rank];
            }
        }
    }
    printf("\n");

    printf("Max_len + Odd Cycle level-1 with n = %d\n", n);

    printf("Number of bad permutation: %d\n", count);
}

void printBadTranslocationFromIdentityCombined_Level2(int n, int *distance_array)
{
    int *pi = (int *)malloc(n * sizeof(int));
    int *pi_inv = (int *)malloc(n * sizeof(int));

    long long size = factorial(n);

    initialize_identity_permutation(pi, n);
    // print_array(pi, n);
    compute_inverse(pi, pi_inv, n);
    // print_array(pi_inv, n);

    int pid = rank_safe(n, pi, pi_inv);

    // printf("PiD: %d\n", pid);
    int count = 0;
    long long last_progress = -1;

    for (int index = 0; index < size; ++index)
    // for (int index = 3; index < size - 1; ++index)
    {
        int current_distance;
        int current_sub_len;
        int current_maxCycle;
        int neighbor_distance;
        int neighbor_maxCycle;
        int neighbor_sub_len;
        int result[MAX_N];
        int pi_shifted[MAX_N];

        initialize_identity_permutation(result, n);

        if (distance_array[index] == 0)
            continue;
        // if (index != pid)
        {

            long long progress = (index * 100) / size;
            if (progress != last_progress)
            {
                printf("\rProgress: %lld%% (%d/%lld)", progress, index, size);
                fflush(stdout);
                last_progress = progress;
            }

            initialize_identity_permutation(pi, n);

            unrank1(n, index, pi);
            // printf("Current permutation");
            // print_array(pi, n);
            // printf("\n");

            //------------ Added
            int (*adj)[2] = calloc(2 * n + 1, sizeof *adj);
            int *deg = calloc(2 * n + 1, sizeof *deg);

            int tmp[MAX_N];
            int tmp_inv[MAX_N];
            int tmp_shift[MAX_N];
            int max_rank = -1;
            int max_len = -1, max_len_2 = -1;
            int max_cycle = -1, max_cycle_2 = -1;

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

                        neighbor_sub_len = longest_increasing_consecutive_values(tmp, n);
                        // neighbor_sub_len = longest_increasing_subsequence(tmp, n);

                        // compute cycle now so it's in sync
                        reset_graph(n, adj, deg);

                        shift_permutation_by_one(tmp, tmp_shift, n);
                        reset_graph(n, adj, deg);
                        build_black_edges_from_perm(adj, deg, tmp_shift, n);
                        build_gray_edges_identity(adj, deg, n);
                        neighbor_maxCycle = count_odd_cycles(adj, n);

                        // ------------End of Update max_cycle ----------

                        if (neighbor_sub_len > max_len || (neighbor_sub_len == max_len && neighbor_maxCycle > max_cycle))
                        {
                            // ----------------Update max_rank
                            compute_inverse(tmp, tmp_inv, n);
                            max_rank = rank_safe(n, tmp, tmp_inv);

                            // -----------------Update max_len
                            max_len = neighbor_sub_len;

                            // -----------------Update max_Cycle
                            max_cycle = neighbor_maxCycle;

                            BestNeighborMetrics level2_metrics = findBestNeighborMetrics(tmp, n);
                            max_len_2 = level2_metrics.max_len;
                            max_cycle_2 = level2_metrics.max_cycles;
                        }

                        else if (neighbor_sub_len == max_len && neighbor_maxCycle == max_cycle)
                        {

                            BestNeighborMetrics level2_metrics = findBestNeighborMetrics(tmp, n);
                            int neighbor_sub_len_2 = level2_metrics.max_len;
                            int neighbor_maxCycle_2 = level2_metrics.max_cycles;

                            // ------------End of calculate current neighbor_maxCycle ----------

                            if (neighbor_sub_len_2 > max_len_2 || (neighbor_sub_len_2 == max_len_2 && neighbor_maxCycle_2 > max_cycle_2))
                            {
                                // -----------------Update max_rank----------
                                compute_inverse(tmp, tmp_inv, n);
                                max_rank = rank_safe(n, tmp, tmp_inv);
                                // ----------------Update max_len_2
                                max_len_2 = neighbor_sub_len_2;

                                // -----------------Update max_cycle_2 ----------
                                max_cycle_2 = neighbor_maxCycle_2;
                                // -----------------Update max_cycle_2 ----------
                            }
                            else if (neighbor_sub_len_2 == max_len_2 && neighbor_maxCycle_2 == max_cycle_2)
                            {
                                // TO be design later
                            }
                        }
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
    printf("\n");
    printf("Max_len + Odd Cycle level-2 with n = %d\n", n);

    printf("Number of bad permutation: %d\n", count);
}
