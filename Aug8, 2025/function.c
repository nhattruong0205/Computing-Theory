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
int compare_long_long(const void *a, const void *b)
{
    long long arg1 = *(const long long *)a;
    long long arg2 = *(const long long *)b;

    if (arg1 < arg2)
        return -1;
    if (arg1 > arg2)
        return 1;
    return 0;
}

// Helper function to generate next permutation in lexicographic order
bool next_permutation(int *arr, int n)
{
    // Find the largest index i such that arr[i] < arr[i + 1]
    int i = n - 2;
    while (i >= 0 && arr[i] >= arr[i + 1])
        i--;

    // If no such index exists, we're at the last permutation
    if (i < 0)
        return false;

    // Find the largest index j greater than i such that arr[i] < arr[j]
    int j = n - 1;
    while (arr[j] <= arr[i])
        j--;

    // Swap arr[i] and arr[j]
    int temp = arr[i];
    arr[i] = arr[j];
    arr[j] = temp;

    // Reverse the suffix starting at arr[i + 1]
    int left = i + 1;
    int right = n - 1;
    while (left < right)
    {
        temp = arr[left];
        arr[left] = arr[right];
        arr[right] = temp;
        left++;
        right--;
    }

    return true;
}

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
int *load_D_from_file(int n, long long *size_out, const char *rank_name)
{
    char filepath[512];
    snprintf(filepath, sizeof(filepath),
             "/Users/nhattruong/Documents/ComputingTheoryDArraydistance%sRank/distances_n%d.txt", rank_name, n);

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

// ===================Ranking permutation into index=====================
// Rank_name
// === 1. LehmerAscendingRadix ===
// === 2. Lex ===
// === 3. Lehmer ===
// === 5. ReverseColexOrder===
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

void unrank2(int n, int r, int pi[])
{
    if (n > 0)
    {
        int s = r / factorial(n - 1);
        swap(&pi[n - 1], &pi[s]);
        unrank2(n - 1, r % factorial(n - 1), pi);
    }
}

int rankSJT(int n, int *pi)
{
    int j, k, r;

    if (n == 1)
        return 0;

    // Find position of largest element (n-1 in 0-based)
    j = 0;
    while (pi[j] != n - 1)
    {
        j++;
    }

    // Count elements to the left that are smaller than n-1
    k = 0;
    for (int i = 0; i < j; i++)
    {
        if (pi[i] < n - 1)
        {
            k++;
        }
    }

    // Create reduced permutation by removing the largest element
    // and reducing all elements > pi[i] by 1
    int *reduced_pi = (int *)malloc((n - 1) * sizeof(int));
    if (!reduced_pi)
    {
        printf("Memory allocation failed in rankSJT\n");
        return -1;
    }

    int idx = 0;
    for (int i = 0; i < n; i++)
    {
        if (pi[i] != n - 1)
        { // Skip the largest element
            // Reduce elements larger than current by 1
            reduced_pi[idx] = pi[i] > (n - 1) ? pi[i] - 1 : pi[i];
            idx++;
        }
    }

    // Recursive call
    r = rankSJT(n - 1, reduced_pi);
    free(reduced_pi);

    // Apply the SJT ranking formula
    if (r % 2 == 1)
    { // r is odd
        return n * r + k;
    }
    else
    { // r is even
        return n * r + n - 1 - k;
    }
}

// Algorithm 5.14: Unranking algorithm for SJT algorithm
void unrankSJT(int n, int r, int *pi, int *dir)
{
    int j, k, rem, c;

    // Initialize pi array to 0
    for (j = 0; j < n; j++)
    {
        pi[j] = 0;
    }

    // Process from n down to 1
    for (j = n; j >= 1; j--)
    {
        rem = r % j;
        r = r / j;

        if (r % 2 == 1)
        {                   // r is odd
            k = -1;         // Start at -1 so first increment gives 0
            dir[j - 1] = 1; // Moving right
        }
        else
        {
            k = n;           // Start at n so first increment gives n-1
            dir[j - 1] = -1; // Moving left
        }

        c = -1;

        do
        {
            k = k + dir[j - 1];

            // Add bounds checking to prevent segmentation fault
            if (k < 0 || k >= n)
            {
                printf("Error: k out of bounds in unrankSJT: k=%d, n=%d, j=%d\n", k, n, j);
                return;
            }

            if (pi[k] == 0)
                c = c + 1;
        } while (c != rem);

        pi[k] = j - 1; // Convert to 0-based: store j-1 instead of j
    }
}
// Reverse Colex order
int rankReverseColexOrder(int n, int pi[])
{
    int r = 0;
    int f = 1;
    for (int i = 1; i < n; i++)
    {
        int c = 0;
        for (int j = 0; j < i; j++)
        {
            if (pi[j] > pi[i])
                c++;
        }
        r += c * f;
        f *= (i + 1);
    }
    return r;
}

void unrankReverseColexOrder(int n, int r, int pi[])
{
    int fact[32];
    fact[0] = 1;
    for (int i = 1; i <= n; i++)
        fact[i] = fact[i - 1] * i;

    if (r < 0 || r >= fact[n])
    {
        printf("rank out of range (0..%d)\n", fact[n] - 1);
        return;
    }

    int c[32];
    c[0] = 0;
    for (int i = 1; i < n; i++)
        c[i] = (r / fact[i]) % (i + 1);

    int avail[32];
    for (int i = 0; i < n; i++)
        avail[i] = i; // ascending
    int avail_size = n;

    for (int i = n - 1; i >= 0; i--)
    {
        // pick (c[i]+1)-th largest from avail
        int idx = avail_size - (c[i] + 1); // index in ascending array
        pi[i] = avail[idx];

        // remove element at idx
        for (int j = idx; j < avail_size - 1; j++)
            avail[j] = avail[j + 1];
        avail_size--;
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
int L(int x) { return 2 * x - 1; } // left endpoint
int R(int x) { return 2 * x; }     // right endpoint

void add_edge(int (*adj)[2], int *deg, int u, int v, const char *label)
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

void build_black_edges_from_perm(int (*adj)[2], int *deg, int *pi, int n)
{
    for (int i = 0; i < n; i++)
    {
        int a = pi[i];
        int b = pi[(i + 1) % n]; // circular adjacency
        add_edge(adj, deg, R(a), L(b), "black");
    }
}

void build_gray_edges_identity(int (*adj)[2], int *deg, int n)
{
    for (int i = 1; i <= n; i++)
    {
        int j = (i < n) ? (i + 1) : 1; // circular identity adjacency
        add_edge(adj, deg, R(i), L(j), "gray");
    }
}

void count_cycles_colored(int (*adj)[2], int n)
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
int *ComputeTDistanceFromIdentity(int n, const char *rank_name)
{
    int *pi = (int *)malloc(n * sizeof(int));

    if (!pi)
    {
        printf("Failed to allocate memory for pi\n");
        return NULL;
    }

    int *pi_inv = malloc(n * sizeof(int));
    if (!pi_inv)
    {
        printf("Failed to allocate pi_inv\n");
        free(pi);
        return NULL;
    }

    long long FACT = factorial(n);

    printf("Clearing arrays for %s computation...\n", rank_name);
    for (long long i = 0; i < FACT; i++)
    {
        D[i] = -1; // Use -1 to indicate uncomputed distances
        visited[i] = false;
    }

    int *result_dir = (int *)malloc(n * sizeof(int));

    // Identity permutation
    initialize_identity_permutation(pi, n);
    compute_inverse(pi, pi_inv, n);

    int pid;
    if (strcmp(rank_name, "Lex") == 0)
        pid = rank_lex(pi, n);
    else if (strcmp(rank_name, "Lehmer") == 0)
        pid = rank2_safe(n, pi, pi_inv);
    else if (strcmp(rank_name, "LehmerAscendingRadix") == 0)
        pid = rank_safe(n, pi, pi_inv);
    else if (strcmp(rank_name, "SJT") == 0)
        pid = rankSJT(n, pi);
    else if (strcmp(rank_name, "ReverseColexOrder") == 0)
        pid = rankReverseColexOrder(n, pi);

    D[pid] = 0;
    visited[pid] = true;

    initQueue();
    enqueue(pid);

    long long processed = 0;
    int *result = (int *)malloc(n * sizeof(int));
    int *tmp = (int *)malloc(n * sizeof(int));
    int *tmp_inv = (int *)malloc(n * sizeof(int));
    while (!isEmpty())
    {
        processed++;
        if (processed % 100000 == 0)
        {
            printf("Processed %lld permutations, queue size: %lld\n", processed, queueSize());
        }

        int current_rank = dequeue();

        // Convert rank back to permutation
        if (strcmp(rank_name, "Lex") == 0)
            unrank_lex(n, current_rank, result);
        else if (strcmp(rank_name, "Lehmer") == 0)
            unrank2(n, current_rank, result);
        else if (strcmp(rank_name, "LehmerAscendingRadix") == 0)
            unrank1(n, current_rank, result);
        else if (strcmp(rank_name, "SJT") == 0)
            unrankSJT(n, current_rank, result, result_dir);
        else if (strcmp(rank_name, "ReverseColexOrder") == 0)
            unrankReverseColexOrder(n, current_rank, result);

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

                    compute_inverse(tmp, tmp_inv, n);
                    // Get rank of translocated permutation
                    int rank_tmp;
                    // Convert rank back to permutation
                    if (strcmp(rank_name, "Lex") == 0)
                        rank_tmp = rank_lex(tmp, n);
                    else if (strcmp(rank_name, "Lehmer") == 0)
                        rank_tmp = rank2_safe(n, tmp, tmp_inv);
                    else if (strcmp(rank_name, "LehmerAscendingRadix") == 0)
                        rank_tmp = rank_safe(n, tmp, tmp_inv);
                    else if (strcmp(rank_name, "SJT") == 0)
                        rank_tmp = rankSJT(n, tmp);
                    else if (strcmp(rank_name, "ReverseColexOrder") == 0)
                        rank_tmp = rankReverseColexOrder(n, tmp);

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
                        // printf("Rank: %d\n", rank_tmp);
                        // printf("Distance: %d\n", D[rank_tmp]);
                        enqueue(rank_tmp);
                    }
                }
    }

    printf("Total processed: %lld permutations\n", processed);

    // Comment after
    char filename[512];
    snprintf(filename, sizeof(filename),
             "/Users/nhattruong/Documents/ComputingTheoryDArraydistance%sRank/distances_n%d.txt", rank_name, n);
    save_D_to_file(filename, D, FACT);

    return D;
}

// ================== Computing PAs =========================
// Compute distance between two permutations pi and sigma
int distance_between_2_permutations(int n, int *pi, int *sigma, int *D, const char *rank_name)
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

    int r;
    if (strcmp(rank_name, "Lex") == 0)
        r = rank_lex(composed, n);
    else if (strcmp(rank_name, "Lehmer") == 0)
        r = rank2_safe(n, composed, composed_inv);
    else if (strcmp(rank_name, "LehmerAscendingRadix") == 0)
        r = rank_safe(n, composed, composed_inv);
    else if (strcmp(rank_name, "SJT") == 0)
        r = rankSJT(n, composed);
    else if (strcmp(rank_name, "ReverseColexOrder") == 0)
        r = rankReverseColexOrder(n, composed);

    return D[r]; // lookup precomputed distance
}

// Computing T(n,d) PA - an array A of permutation on [1..n] with dt(A) >= d.
long long T_n_d(int n, int d, int *D, const char *rank_name)
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
    int *result_dir = (int *)malloc(n * sizeof(int));

    // Greedy algorithm: keep selecting permutations until none remain
    for (long long i = 0; i < factorial(n); i++)
    {
        if (!forbidden[i])
        {
            // Select this permutation as a codeword
            code_size++;

            // Convert rank i to permutation pi
            initialize_identity_permutation(pi, n);

            if (strcmp(rank_name, "Lex") == 0)
                unrank_lex(n, (int)i, pi);
            else if (strcmp(rank_name, "Lehmer") == 0)
                unrank2(n, (int)i, pi);
            else if (strcmp(rank_name, "LehmerAscendingRadix") == 0)
                unrank1(n, (int)i, pi);
            else if (strcmp(rank_name, "SJT") == 0)
                unrankSJT(n, (int)i, pi, result_dir);
            else if (strcmp(rank_name, "ReverseColexOrder") == 0)
                unrankReverseColexOrder(n, (int)i, pi);
            // print_array(pi, n);

            // Forbid all permutations within distance d-1 of pi
            for (long long j = 0; j < factorial(n); j++)
            {
                if (!forbidden[j])
                {
                    int dist;
                    // Convert rank j to permutation sigma
                    initialize_identity_permutation(sigma, n);

                    if (strcmp(rank_name, "Lex") == 0)
                    {
                        unrank_lex(n, (int)j, sigma);
                        dist = distance_between_2_permutations(n, pi, sigma, D, rank_name);
                    }
                    else if (strcmp(rank_name, "Lehmer") == 0)
                    {
                        unrank2(n, (int)j, pi);
                        dist = distance_between_2_permutations(n, pi, sigma, D, rank_name);
                    }
                    else if (strcmp(rank_name, "LehmerAscendingRadix") == 0)
                    {
                        unrank1(n, (int)j, sigma);
                        // Compute distance between pi and sigma
                        dist = distance_between_2_permutations(n, pi, sigma, D, rank_name);
                    }
                    else if (strcmp(rank_name, "SJT") == 0)
                    {
                        unrankSJT(n, (int)j, sigma, result_dir);
                        // Compute distance between pi and sigma
                        dist = distance_between_2_permutations(n, pi, sigma, D, rank_name);
                    }
                    else if (strcmp(rank_name, "ReverseColexOrder") == 0)
                    {
                        unrankReverseColexOrder(n, (int)j, sigma);
                        // Compute distance between pi and sigma
                        dist = distance_between_2_permutations(n, pi, sigma, D, rank_name);
                    }

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

// Check if rank in selected
bool is_rank_in_selected(long long *selected, int size, long long rank)
{
    // LINEAR SEARCH - works on unsorted arrays
    for (int i = 0; i < size; i++)
    {
        if (selected[i] == rank)
        {
            return true; // Found
        }
    }
    return false; // Not found
}

long long *compute_ranks(int n, int perms[][n], int num_perms)
{
    long long *selected = (long long *)malloc(num_perms * sizeof(long long));

    for (int i = 0; i < num_perms; i++)
    {
        int pi_working[n];
        int pi_inv[n];

        // Copy permutation
        for (int j = 0; j < n; j++)
        {
            pi_working[j] = perms[i][j];
        }

        // Compute inverse
        compute_inverse(pi_working, pi_inv, n);

        // Compute and store rank
        selected[i] = rank_safe(n, pi_working, pi_inv);
    }

    return selected;
}

bool can_add_to_code_incremental_d2(int n, int *pi, long long *selected, int code_size, const char *rank_name)
{
    int neighbor[MAX_N];
    int neighbor_inv[MAX_N]; // Add this for inverse computation
    int *result_dir = (int *)malloc(n * sizeof(int));
    int pi_inv[n];
    compute_inverse(pi, pi_inv, n);

    // First check if pi itself is already in A
    long long pi_rank;
    if (strcmp(rank_name, "Lex") == 0)
        pi_rank = rank_lex(pi, n);
    else if (strcmp(rank_name, "Lehmer") == 0)
        pi_rank = rank2_safe(n, pi, pi_inv);
    else if (strcmp(rank_name, "LehmerAscendingRadix") == 0)
        pi_rank = rank_safe(n, pi, pi_inv);
    else if (strcmp(rank_name, "SJT") == 0)
        pi_rank = rankSJT(n, pi);
    else if (strcmp(rank_name, "ReverseColexOrder") == 0)
        pi_rank = rankReverseColexOrder(n, pi);

    if (is_rank_in_selected(selected, code_size, pi_rank))
    {
        free(result_dir);
        return false; // pi is already in A
    }

    for (int i = 0; i < n; ++i)
    {
        for (int j = i + 1; j < n; ++j)
        {
            for (int k = j; k < n; ++k)
            {
                /* Build translocated permutation */
                int idx = 0;

                // 1. Prefix: [0..i-1]
                for (int x = 0; x < i; ++x)
                    neighbor[idx++] = pi[x];

                // 2. Block: [j..k]
                for (int x = j; x <= k; ++x)
                    neighbor[idx++] = pi[x];

                // 3. Middle: [i..j-1]
                for (int x = i; x < j; ++x)
                    neighbor[idx++] = pi[x];

                // 4. Suffix: [k+1..n-1]
                for (int x = k + 1; x < n; ++x)
                    neighbor[idx++] = pi[x];

                // Compute inverse for neighbor if needed
                compute_inverse(neighbor, neighbor_inv, n);

                // Get rank of this neighbor - FIXED: use 'neighbor' not 'tmp'
                long long neighbor_rank;
                if (strcmp(rank_name, "Lex") == 0)
                    neighbor_rank = rank_lex(neighbor, n);
                else if (strcmp(rank_name, "Lehmer") == 0)
                    neighbor_rank = rank2_safe(n, neighbor, neighbor_inv);
                else if (strcmp(rank_name, "LehmerAscendingRadix") == 0)
                    neighbor_rank = rank_safe(n, neighbor, neighbor_inv);
                else if (strcmp(rank_name, "SJT") == 0)
                    neighbor_rank = rankSJT(n, neighbor);
                else if (strcmp(rank_name, "ReverseColexOrder") == 0)
                    neighbor_rank = rankReverseColexOrder(n, neighbor);

                // Check if this neighbor is in A
                if (is_rank_in_selected(selected, code_size, neighbor_rank))
                {
                    free(result_dir);
                    return false; // A 1-neighbor of pi is already in A
                }
            }
        }
    }

    free(result_dir);
    return true; // Can safely add pi to A
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

// Function to compute a (n,2)-PA using greedy construction with checkpointing
// Returns the code size, fills selected array with ranks only
int compute_n2_PA(int n, long long *selected, int max_size, const char *rank_name)
{
    int code_size = 0;

    char checkpoint_file[256];
    snprintf(checkpoint_file, sizeof(checkpoint_file), "checkpoint_n%d.txt", n);

    // Calculate total permutations
    long long total_perms = 1;
    for (int i = 1; i <= n; i++)
        total_perms *= i;

    // Try to resume from checkpoint
    FILE *fp = fopen(checkpoint_file, "r");
    int current_perm[MAX_N];
    long long last_checked = -1;

    if (fp != NULL)
    {
        printf("Found checkpoint file, resuming...\n");

        // Read code size
        fscanf(fp, "%d\n", &code_size);

        // Read selected ranks
        for (int i = 0; i < code_size; i++)
        {
            fscanf(fp, "%lld\n", &selected[i]);
        }

        // Read last checked permutation rank
        fscanf(fp, "%lld\n", &last_checked);

        fclose(fp);

        printf("Resumed: code_size=%d, last_checked=%lld, total_perms=%lld\n",
               code_size, last_checked, total_perms);

        // VALIDATE BEFORE UNRANKING
        if (last_checked >= total_perms)
        {
            printf("Error: Invalid last_checked rank %lld (max is %lld)\n",
                   last_checked, total_perms - 1);
            printf("Checkpoint file is corrupted. Please delete checkpoint_n%d.txt and restart.\n", n);
            return -1;
        }

        // Validate all selected ranks
        for (int i = 0; i < code_size; i++)
        {
            if (selected[i] >= total_perms || selected[i] < 0)
            {
                printf("Error: Invalid rank at position %d: %lld (max is %lld)\n",
                       i, selected[i], total_perms - 1);
                printf("Checkpoint file is corrupted. Please delete checkpoint_n%d.txt and restart.\n", n);
                return -1;
            }
        }

        // Check if we've already completed
        if (last_checked >= total_perms - 1)
        {
            printf("Already completed! All permutations checked.\n");
            printf("Final code size: %d\n", code_size);
            return code_size;
        }

        // Unrank to get the starting permutation
        if (last_checked >= 0)
        {
            // Initialize with identity permutation FIRST
            for (int i = 0; i < n; i++)
            {
                current_perm[i] = i;
            }
            // Then unrank
            unrank1(n, last_checked, current_perm);

            // Advance to next permutation
            if (!next_permutation(current_perm, n))
            {
                printf("Already completed!\n");
                printf("Final code size: %d\n", code_size);
                return code_size;
            }
        }
        else
        {
            for (int i = 0; i < n; i++)
                current_perm[i] = i;
        }
    }
    else
    {
        printf("No checkpoint found, starting fresh...\n");

        // Start with identity permutation
        int pi[MAX_N];
        for (int i = 0; i < n; i++)
        {
            pi[i] = i;
            current_perm[i] = i;
        }

        // Compute rank of identity and add to selected
        int pi_inv[MAX_N];
        compute_inverse(pi, pi_inv, n);
        selected[0] = rank_safe(n, pi, pi_inv);
        code_size = 1;
    }

    printf("Starting greedy construction...\n");

    long long checked = (last_checked >= 0) ? last_checked + 1 : 0;
    int checkpoint_interval = 10000; // Save every 10000 permutations

    do
    {
        checked++;

        // Check if this permutation can be added
        if (can_add_to_code_incremental_d2(n, current_perm, selected, code_size, rank_name))
        {
            // Compute and store its rank
            int temp_inv[MAX_N];
            compute_inverse(current_perm, temp_inv, n);
            selected[code_size] = rank_safe(n, current_perm, temp_inv);
            code_size++;

            // Check if we've reached max size
            if (code_size >= max_size)
            {
                printf("Reached maximum code size: %d\n", max_size);
                break;
            }
        }

        // Save checkpoint periodically
        if (checked % checkpoint_interval == 0)
        {
            printf("Progress: checked %lld/%lld permutations, code size = %d\n",
                   checked, total_perms, code_size);

            // Save checkpoint
            fp = fopen(checkpoint_file, "w");
            if (fp != NULL)
            {
                fprintf(fp, "%d\n", code_size);
                for (int i = 0; i < code_size; i++)
                {
                    fprintf(fp, "%lld\n", selected[i]);
                }
                // Save current permutation rank (NOT checked counter!)
                int temp_inv[MAX_N];
                compute_inverse(current_perm, temp_inv, n);
                long long current_rank = rank_safe(n, current_perm, temp_inv);
                fprintf(fp, "%lld\n", current_rank);
                fclose(fp);
            }
        }
    } while (next_permutation(current_perm, n) && code_size < max_size);

    printf("\nFinal code size: %d\n", code_size);
    printf("Checked %lld permutations\n", checked);

    // Save final checkpoint
    fp = fopen(checkpoint_file, "w");
    if (fp != NULL)
    {
        fprintf(fp, "%d\n", code_size);
        for (int i = 0; i < code_size; i++)
        {
            fprintf(fp, "%lld\n", selected[i]);
        }
        // Save the rank of the last permutation checked
        fprintf(fp, "%lld\n", total_perms - 1); // Mark as complete with max rank
        fclose(fp);
    }

    return code_size;
}

// // Function to compute a (n,2)-PA using greedy construction with checkpointing
// int compute_n2_PA(int n, int perms[][n], int max_size, const char *rank_name)
// {
//     int code_size = 0;
//     long long *selected = (long long *)malloc(max_size * sizeof(long long));
//     int D[MAX_N] = {0}; // Distance array (not strictly needed for d=2)

//     char checkpoint_file[256];
//     snprintf(checkpoint_file, sizeof(checkpoint_file), "checkpoint_n%d.txt", n);

//     // Calculate total permutations
//     long long total_perms = 1;
//     for (int i = 1; i <= n; i++)
//         total_perms *= i;

//     // Try to resume from checkpoint
//     FILE *fp = fopen(checkpoint_file, "r");
//     int current_perm[MAX_N];
//     long long last_checked = -1;

//     if (fp != NULL)
//     {
//         printf("Found checkpoint file, resuming...\n");

//         // Read code size
//         fscanf(fp, "%d\n", &code_size);

//         // Read selected ranks (but don't unrank yet)
//         for (int i = 0; i < code_size; i++)
//         {
//             fscanf(fp, "%lld\n", &selected[i]);
//         }

//         // Read last checked permutation rank
//         fscanf(fp, "%lld\n", &last_checked);

//         fclose(fp);

//         printf("Resumed: code_size=%d, last_checked=%lld, total_perms=%lld\n",
//                code_size, last_checked, total_perms);

//         // VALIDATE BEFORE UNRANKING
//         if (last_checked >= total_perms)
//         {
//             printf("Error: Invalid last_checked rank %lld (max is %lld)\n",
//                    last_checked, total_perms - 1);
//             printf("Checkpoint file is corrupted. Please delete checkpoint_n%d.txt and restart.\n", n);
//             free(selected);
//             return -1;
//         }

//         // Validate all selected ranks
//         for (int i = 0; i < code_size; i++)
//         {
//             if (selected[i] >= total_perms || selected[i] < 0)
//             {
//                 printf("Error: Invalid rank at position %d: %lld (max is %lld)\n",
//                        i, selected[i], total_perms - 1);
//                 printf("Checkpoint file is corrupted. Please delete checkpoint_n%d.txt and restart.\n", n);
//                 free(selected);
//                 return -1;
//             }
//         }

//         printf("Unranking %d permutations...\n", code_size);

//         // NOW it's safe to unrank
//         for (int i = 0; i < code_size; i++)
//         {
//             // Initialize with identity permutation FIRST
//             for (int j = 0; j < n; j++)
//             {
//                 perms[i][j] = j;
//             }
//             // Then unrank
//             unrank1(n, selected[i], perms[i]);
//         }

//         printf("Unranking complete.\n");

//         // Check if we've already completed
//         if (last_checked >= total_perms - 1)
//         {
//             printf("Already completed! All permutations checked.\n");
//             printf("Final code size: %d\n", code_size);
//             return code_size;
//         }

//         // Unrank to get the starting permutation
//         if (last_checked >= 0)
//         {
//             // Initialize with identity permutation FIRST
//             for (int i = 0; i < n; i++)
//             {
//                 current_perm[i] = i;
//             }
//             // Then unrank
//             unrank1(n, last_checked, current_perm);

//             // Advance to next permutation
//             if (!next_permutation(current_perm, n))
//             {
//                 printf("Already completed!\n");
//                 printf("Final code size: %d\n", code_size);
//                 return code_size;
//             }
//         }
//         else
//         {
//             for (int i = 0; i < n; i++)
//                 current_perm[i] = i;
//         }
//     }
//     else
//     {
//         printf("No checkpoint found, starting fresh...\n");

//         // Start with identity permutation
//         int pi[MAX_N];
//         for (int i = 0; i < n; i++)
//         {
//             pi[i] = i;
//             perms[code_size][i] = i;
//             current_perm[i] = i;
//         }

//         // Compute rank of identity and add to selected
//         int pi_inv[MAX_N];
//         compute_inverse(pi, pi_inv, n);
//         selected[0] = rank_safe(n, pi, pi_inv);
//         code_size = 1;
//     }

//     printf("Starting greedy construction...\n");

//     long long checked = (last_checked >= 0) ? last_checked + 1 : 0;
//     int checkpoint_interval = 10000; // Save every 10000 permutations

//     do
//     {
//         checked++;

//         // Check if this permutation can be added
//         if (can_add_to_code_incremental_d2(n, current_perm, selected, code_size, rank_name))
//         {
//             // Add this permutation to the code
//             for (int i = 0; i < n; i++)
//             {
//                 perms[code_size][i] = current_perm[i];
//             }

//             // Compute and store its rank
//             int temp_inv[MAX_N];
//             compute_inverse(current_perm, temp_inv, n);
//             selected[code_size] = rank_safe(n, current_perm, temp_inv);
//             code_size++;

//             // Check if we've reached max size
//             if (code_size >= max_size)
//             {
//                 printf("Reached maximum code size: %d\n", max_size);
//                 break;
//             }
//         }

//         // Save checkpoint periodically
//         if (checked % checkpoint_interval == 0)
//         {
//             printf("Progress: checked %lld/%lld permutations, code size = %d\n",
//                    checked, total_perms, code_size);

//             // Save checkpoint
//             fp = fopen(checkpoint_file, "w");
//             if (fp != NULL)
//             {
//                 fprintf(fp, "%d\n", code_size);
//                 for (int i = 0; i < code_size; i++)
//                 {
//                     fprintf(fp, "%lld\n", selected[i]);
//                 }
//                 // Save current permutation rank (NOT checked counter!)
//                 int temp_inv[MAX_N];
//                 compute_inverse(current_perm, temp_inv, n);
//                 long long current_rank = rank_safe(n, current_perm, temp_inv);
//                 fprintf(fp, "%lld\n", current_rank);
//                 fclose(fp);
//             }
//         }
//     } while (next_permutation(current_perm, n) && code_size < max_size);

//     printf("\nFinal code size: %d\n", code_size);
//     printf("Checked %lld permutations\n", checked);

//     // Save final checkpoint
//     fp = fopen(checkpoint_file, "w");
//     if (fp != NULL)
//     {
//         fprintf(fp, "%d\n", code_size);
//         for (int i = 0; i < code_size; i++)
//         {
//             fprintf(fp, "%lld\n", selected[i]);
//         }
//         // Save the rank of the last permutation checked
//         fprintf(fp, "%lld\n", total_perms - 1); // Mark as complete with max rank
//         fclose(fp);
//     }

//     free(selected);
//     return code_size;
// }

// // Function to compute a (n,2)-PA using greedy construction
// int compute_n2_PA(int n, int perms[][n], int max_size, const char *rank_name)
// {
//     int code_size = 0;
//     long long *selected = (long long *)malloc(max_size * sizeof(long long));
//     int D[MAX_N] = {0}; // Distance array (not strictly needed for d=2)

//     // Start with identity permutation
//     int pi[MAX_N];
//     for (int i = 0; i < n; i++)
//     {
//         pi[i] = i;
//         perms[code_size][i] = i;
//     }

//     // Compute rank of identity and add to selected
//     int pi_inv[MAX_N];
//     compute_inverse(pi, pi_inv, n);
//     selected[0] = rank_safe(n, pi, pi_inv);
//     code_size = 1;

//     printf("Starting greedy construction...\n");
//     // printf("Added permutation %d: ", code_size);
//     // print_array(pi, n);

//     // Generate permutations and try to add them
//     long long total_perms = 1;
//     for (int i = 1; i <= n; i++)
//         total_perms *= i;

//     int current_perm[MAX_N];
//     for (int i = 0; i < n; i++)
//         current_perm[i] = i;

//     // Try all permutations in lexicographic order
//     long long checked = 0;
//     do
//     {
//         checked++;

//         // Check if this permutation can be added
//         if (can_add_to_code_incremental_d2(n, current_perm, selected, code_size, "LehmerAscendingRadix"))
//         {
//             // Add this permutation to the code
//             for (int i = 0; i < n; i++)
//             {
//                 perms[code_size][i] = current_perm[i];
//             }

//             // Compute and store its rank
//             int temp_inv[MAX_N];
//             compute_inverse(current_perm, temp_inv, n);
//             selected[code_size] = rank_safe(n, current_perm, temp_inv);
//             code_size++;

//             // printf("Added permutation %d: ", code_size);
//             // print_array(current_perm, n);

//             // Check if we've reached max size
//             if (code_size >= max_size)
//             {
//                 printf("Reached maximum code size: %d\n", max_size);
//                 break;
//             }
//         }

//         // Progress indicator
//         if (checked % 10000 == 0)
//         {
//             printf("Progress: checked %lld/%lld permutations, code size = %d\n",
//                    checked, total_perms, code_size);
//         }

//     } while (next_permutation(current_perm, n) && code_size < max_size);

//     printf("\nFinal code size: %d\n", code_size);
//     printf("Checked %lld permutations\n", checked);

//     free(selected);
//     return code_size;
// }

bool verify_n2_PA(int n, const char *rank_name)
{
    char checkpoint_file[256];
    snprintf(checkpoint_file, sizeof(checkpoint_file), "checkpoint_n%d.txt", n);

    // Calculate total permutations
    long long total_perms = 1;
    for (int i = 1; i <= n; i++)
        total_perms *= i;

    // Open checkpoint file
    FILE *fp = fopen(checkpoint_file, "r");
    if (fp == NULL)
    {
        printf("Error: Checkpoint file not found!\n");
        return false;
    }

    // Read code size
    int num_perms;
    fscanf(fp, "%d\n", &num_perms);

    // Allocate memory for ranks
    long long *selected = (long long *)malloc(num_perms * sizeof(long long));
    if (!selected)
    {
        printf("Error: Memory allocation failed!\n");
        fclose(fp);
        return false;
    }

    // Read ONLY the code_size ranks (NOT the last_checked marker)
    for (int i = 0; i < num_perms; i++)
    {
        fscanf(fp, "%lld\n", &selected[i]);
    }
    // printf("Loaded codeword ranks:\n");
    // for (int i = 0; i < num_perms; i++)
    // {
    //     printf("  A[%d] = %lld\n", i, selected[i]);
    // }
    // We don't read the last line (progress marker) - it's not a codeword
    fclose(fp);

    printf("Verifying %d permutations...\n", num_perms);

    // Validate all selected ranks
    for (int i = 0; i < num_perms; i++)
    {
        if (selected[i] >= total_perms || selected[i] < 0)
        {
            printf("Error: Invalid rank at position %d: %lld (max is %lld)\n",
                   i, selected[i], total_perms - 1);
            free(selected);
            return false;
        }
    }

    // Sort the ranks array for binary search
    qsort(selected, num_perms, sizeof(long long), compare_long_long);

    // Step 2: Check for duplicates
    for (int i = 0; i < num_perms - 1; i++)
    {
        if (selected[i] == selected[i + 1])
        {
            printf("Error: Found duplicate rank %lld at positions %d and %d\n",
                   selected[i], i, i + 1);
            free(selected);
            return false;
        }
    }

    printf("No duplicates found. Checking transpositions...\n");

    // Step 5-9: For each permutation, check all transpositions
    for (int perm_idx = 0; perm_idx < num_perms; perm_idx++)
    {
        // Unrank to get the permutation
        int pi[MAX_N];
        for (int i = 0; i < n; i++)
        {
            pi[i] = i;
        }
        unrank1(n, selected[perm_idx], pi);

        // Progress indicator
        if ((perm_idx + 1) % 100 == 0 || perm_idx == 0)
        {
            printf("Checking permutation %d/%d (rank %lld)...\n",
                   perm_idx + 1, num_perms, selected[perm_idx]);
        }

        // Step 6: Generate all (i,j,k)-transpositions
        for (int i = 0; i < n; i++)
        {
            for (int j = i + 1; j < n; j++)
            {
                for (int k = j; k < n; k++)
                {
                    // Step 7: Build transposed permutation τ
                    int tau[MAX_N];
                    int idx = 0;

                    // 1. Prefix: [0..i-1]
                    for (int x = 0; x < i; x++)
                        tau[idx++] = pi[x];

                    // 2. Block: [j..k]
                    for (int x = j; x <= k; x++)
                        tau[idx++] = pi[x];

                    // 3. Middle: [i..j-1]
                    for (int x = i; x < j; x++)
                        tau[idx++] = pi[x];

                    // 4. Suffix: [k+1..n-1]
                    for (int x = k + 1; x < n; x++)
                        tau[idx++] = pi[x];

                    // Compute rank of τ
                    int tau_inv[MAX_N];
                    compute_inverse(tau, tau_inv, n);

                    long long tau_rank;
                    if (strcmp(rank_name, "Lex") == 0)
                        tau_rank = rank_lex(tau, n);
                    else if (strcmp(rank_name, "Lehmer") == 0)
                        tau_rank = rank2_safe(n, tau, tau_inv);
                    else if (strcmp(rank_name, "LehmerAscendingRadix") == 0)
                        tau_rank = rank_safe(n, tau, tau_inv);
                    else if (strcmp(rank_name, "SJT") == 0)
                        tau_rank = rankSJT(n, tau);
                    else if (strcmp(rank_name, "ReverseColexOrder") == 0)
                        tau_rank = rankReverseColexOrder(n, tau);

                    // Step 8: Binary search for τ in A
                    // If found, it means a transposition is in the code - INVALID!
                    if (is_rank_in_selected(selected, num_perms, tau_rank))
                    {
                        printf("Error: Found transposition in code!\n");
                        printf("Permutation at index %d (rank %lld): ", perm_idx, selected[perm_idx]);
                        print_array(pi, n);
                        printf("Transposition (i=%d, j=%d, k=%d) with rank %lld: ", i, j, k, tau_rank);
                        print_array(tau, n);
                        free(selected);
                        return false;
                    }
                }
            }
        }
    }

    // Step 10: All checks passed
    printf("All checks passed!\n");
    free(selected);
    return true;
}

// bool verify_n2_PA(int n, int perms[][n], int num_perms, const char *rank_name)
// {
//     // Step 1: Sort A in lexicographical order (using ranks)
//     long long *selected = compute_ranks(n, perms, num_perms);

//     // Sort the ranks array
//     qsort(selected, num_perms, sizeof(long long), compare_long_long);

//     // Step 2: Check for duplicates
//     for (int i = 0; i < num_perms - 1; i++)
//     {
//         if (selected[i] == selected[i + 1])
//         {
//             free(selected);
//             return false; // Found duplicate
//         }
//     }

//     // Step 5-9: For each permutation, check all transpositions
//     for (int perm_idx = 0; perm_idx < num_perms; perm_idx++)
//     {
//         int *pi = perms[perm_idx];

//         // Step 6: Generate all (i,j,k)-transpositions
//         for (int i = 0; i < n; i++)
//         {
//             for (int j = i + 1; j < n; j++)
//             {
//                 for (int k = j; k < n; k++)
//                 {
//                     // Step 7: Build transposed permutation τ
//                     int tau[MAX_N];
//                     int idx = 0;

//                     // 1. Prefix: [0..i-1]
//                     for (int x = 0; x < i; x++)
//                         tau[idx++] = pi[x];

//                     // 2. Block: [j..k]
//                     for (int x = j; x <= k; x++)
//                         tau[idx++] = pi[x];

//                     // 3. Middle: [i..j-1]
//                     for (int x = i; x < j; x++)
//                         tau[idx++] = pi[x];

//                     // 4. Suffix: [k+1..n-1]
//                     for (int x = k + 1; x < n; x++)
//                         tau[idx++] = pi[x];

//                     // Compute rank of τ
//                     int tau_inv[MAX_N];
//                     compute_inverse(tau, tau_inv, n);

//                     long long tau_rank;
//                     if (strcmp(rank_name, "Lex") == 0)
//                         tau_rank = rank_lex(tau, n);
//                     else if (strcmp(rank_name, "Lehmer") == 0)
//                         tau_rank = rank2_safe(n, tau, tau_inv);
//                     else if (strcmp(rank_name, "LehmerAscendingRadix") == 0)
//                         tau_rank = rank_safe(n, tau, tau_inv);
//                     else if (strcmp(rank_name, "SJT") == 0)
//                         tau_rank = rankSJT(n, tau);
//                     else if (strcmp(rank_name, "ReverseColexOrder") == 0)
//                         tau_rank = rankReverseColexOrder(n, tau);

//                     // Step 8: Binary search for τ in A
//                     if (is_rank_in_selected(selected, num_perms, tau_rank))
//                     {
//                         free(selected);
//                         return false; // τ not found in A
//                     }
//                 }
//             }
//         }
//     }

//     // Step 10: All checks passed
//     free(selected);
//     return true;
// }

// bool verify_n2_PA(int n, const char *rank_name)
// {
//     char checkpoint_file[256];
//     snprintf(checkpoint_file, sizeof(checkpoint_file), "checkpoint_n%d.txt", n);

//     FILE *fp = fopen(checkpoint_file, "r");
//     if (fp == NULL)
//     {
//         printf("Error: Could not open checkpoint file %s\n", checkpoint_file);
//         return false;
//     }

//     // Read code size
//     int num_perms;
//     if (fscanf(fp, "%d\n", &num_perms) != 1)
//     {
//         printf("Error: Failed to read code size\n");
//         fclose(fp);
//         return false;
//     }

//     printf("Reading %d permutations from checkpoint...\n", num_perms);

//     // Allocate memory for selected ranks
//     long long *selected = (long long *)malloc(num_perms * sizeof(long long));
//     if (selected == NULL)
//     {
//         printf("Error: Memory allocation failed\n");
//         fclose(fp);
//         return false;
//     }

//     // Read all ranks
//     for (int i = 0; i < num_perms; i++)
//     {
//         if (fscanf(fp, "%lld\n", &selected[i]) != 1)
//         {
//             printf("Error: Failed to read rank %d\n", i);
//             free(selected);
//             fclose(fp);
//             return false;
//         }
//     }
//     fclose(fp);

//     printf("Successfully read %d ranks from checkpoint\n", num_perms);

//     // Sort the ranks array
//     qsort(selected, num_perms, sizeof(long long), compare_long_long);

//     // Step 2: Check for duplicates
//     printf("Checking for duplicates...\n");
//     for (int i = 0; i < num_perms - 1; i++)
//     {
//         if (selected[i] == selected[i + 1])
//         {
//             printf("Error: Found duplicate rank at positions %d and %d: %lld\n",
//                    i, i + 1, selected[i]);
//             free(selected);
//             return false;
//         }
//     }
//     printf("No duplicates found.\n");

//     // Step 5-9: For each permutation, check all transpositions
//     printf("Verifying (n,2)-PA property...\n");

//     for (int perm_idx = 0; perm_idx < num_perms; perm_idx++)
//     {
//         // Unrank to get the permutation
//         int pi[MAX_N];
//         unrank1(n, selected[perm_idx], pi);

//         // Progress indicator
//         if (perm_idx % 100 == 0)
//         {
//             printf("Verified %d/%d permutations...\n", perm_idx, num_perms);
//         }

//         // Step 6: Generate all (i,j,k)-transpositions
//         for (int i = 0; i < n; i++)
//         {
//             for (int j = i + 1; j < n; j++)
//             {
//                 for (int k = j; k < n; k++)
//                 {
//                     // Step 7: Build transposed permutation τ
//                     int tau[MAX_N];
//                     int idx = 0;

//                     // 1. Prefix: [0..i-1]
//                     for (int x = 0; x < i; x++)
//                         tau[idx++] = pi[x];

//                     // 2. Block: [j..k]
//                     for (int x = j; x <= k; x++)
//                         tau[idx++] = pi[x];

//                     // 3. Middle: [i..j-1]
//                     for (int x = i; x < j; x++)
//                         tau[idx++] = pi[x];

//                     // 4. Suffix: [k+1..n-1]
//                     for (int x = k + 1; x < n; x++)
//                         tau[idx++] = pi[x];

//                     // Compute rank of τ
//                     int tau_inv[MAX_N];
//                     compute_inverse(tau, tau_inv, n);

//                     long long tau_rank;
//                     if (strcmp(rank_name, "Lex") == 0)
//                         tau_rank = rank_lex(tau, n);
//                     else if (strcmp(rank_name, "Lehmer") == 0)
//                         tau_rank = rank2_safe(n, tau, tau_inv);
//                     else if (strcmp(rank_name, "LehmerAscendingRadix") == 0)
//                         tau_rank = rank_safe(n, tau, tau_inv);
//                     else if (strcmp(rank_name, "SJT") == 0)
//                         tau_rank = rankSJT(n, tau);
//                     else if (strcmp(rank_name, "ReverseColexOrder") == 0)
//                         tau_rank = rankReverseColexOrder(n, tau);
//                     else
//                     {
//                         printf("Error: Unknown rank name %s\n", rank_name);
//                         free(selected);
//                         return false;
//                     }

//                     // Step 8: Binary search for τ in A
//                     if (is_rank_in_selected(selected, num_perms, tau_rank))
//                     {
//                         printf("Error: Found transposition in code!\n");
//                         printf("Permutation %d: ", perm_idx);
//                         print_array(pi, n);
//                         printf("Transposition (%d,%d,%d): ", i, j, k);
//                         print_array(tau, n);
//                         printf("Rank: %lld\n", tau_rank);
//                         free(selected);
//                         return false;
//                     }
//                 }
//             }
//         }
//     }

//     // Step 10: All checks passed
//     printf("\nVerification successful!\n");
//     printf("Code size: %d\n", num_perms);
//     printf("All permutations satisfy the (n,2)-PA property.\n");

//     free(selected);
//     return true;
// }
