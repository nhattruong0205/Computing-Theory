// project_v7 where I re-implement the figure 2

#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

#define MAX_N 20          // max length of a permutation (increased from 10)
#define MAX_QUEUE 4000000 // adjust if needed

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
    printf("%s: %d -- %d\n", label, u, v);
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

// --------------- BFS-------------

// Map permutation to unique integer
long long factorial[MAX_N + 1];

void init_factorial(int n)
{
    factorial[0] = 1;
    for (int i = 1; i <= n; i++)
        factorial[i] = factorial[i - 1] * i;
}

// Convert permutation to unique index
long long perm_to_index(int *perm, int n)
{
    int used[MAX_N + 1] = {0};
    long long idx = 0;
    for (int i = 0; i < n; i++)
    {
        int cnt = 0;
        for (int j = 1; j < perm[i]; j++)
            if (!used[j])
                cnt++;
        idx += cnt * factorial[n - i - 1];
        used[perm[i]] = 1;
    }
    return idx;
}

// BFS node structurew
typedef struct Node
{
    int perm[MAX_N];
    int trans_count;
    int last_i, last_j, last_k; // the last translocation that created this node
    struct Node *parent;        // pointer to parent node
} Node;

void print_translocation_path(Node *node, int n)
{
    if (node->parent != NULL)
    {
        print_translocation_path(node->parent, n); // recurse to original
        printf("Translocation: i=%d, j=%d, k=%d\n", node->last_i, node->last_j, node->last_k);
        print_array(node->perm, n);
    }
}

// Standard bfs queue
Node *queue;
int front = 0, rear = 0;

// BFS function
// ------------------- BFS function -------------------
void computeMostOddCycle_BFS(int *start_perm, int n)
{
    init_factorial(n);

    bool *visited = calloc(factorial[n], sizeof(bool));
    queue = malloc(MAX_QUEUE * sizeof(Node));
    if (!visited || !queue)
    {
        fprintf(stderr, "Memory allocation failed!\n");
        return;
    }

    Node best;
    int max_odd = -1;
    bool stop = false;

    // Initialize BFS
    Node root;
    memcpy(root.perm, start_perm, n * sizeof(int));
    root.trans_count = 0;
    queue[rear++] = root;
    visited[perm_to_index(root.perm, n)] = true;

    // adjacency arrays
    int (*adj)[2] = calloc(2 * n + 1, sizeof *adj);
    int *deg = calloc(2 * n + 1, sizeof *deg);

    while (front < rear && !stop)
    {
        Node cur = queue[front++];

        memset(adj, 0, (2 * n + 1) * 2 * sizeof(int));
        memset(deg, 0, (2 * n + 1) * sizeof(int));

        build_black_edges_from_perm(adj, deg, cur.perm, n);
        build_gray_edges_identity(adj, deg, n);

        int odd = count_odd_cycles(adj, n);

        if (odd > max_odd)
        {
            max_odd = odd;
            best = cur;
        }

        // Generate all possible translocations
        for (int i = 0; i < n && !stop; i++)
        {
            for (int j = i + 1; j < n && !stop; j++)
            {
                for (int k = j; k < n && !stop; k++)
                {
                    Node next;
                    translocate(cur.perm, next.perm, n, i, j, k);
                    next.trans_count = cur.trans_count + 1;
                    next.last_i = i;
                    next.last_j = j;
                    next.last_k = k;
                    next.parent = &queue[front - 1]; // point to current node
                    long long idx = perm_to_index(next.perm, n);
                    if (!visited[idx])
                    {
                        if (rear >= MAX_QUEUE)
                        {
                            printf("Queue overflow!\n");
                            stop = true;
                            break;
                        }
                        visited[idx] = true;
                        queue[rear++] = next;
                    }
                }
            }
        }
    }

    printf("Best permutation (max odd cycles = %d) reached with %d translocations:\n",
           max_odd, best.trans_count);
    print_array(best.perm, n);

    printf("Sequence of translocations from original:\n");
    print_translocation_path(&best, n);

    free(visited);
    free(queue);
    free(adj);
    free(deg);
}

//-------------- Main --------------------
// int main()
// {
//     int size;

//     int *dst = (int *)malloc(size * sizeof(int));
//     // Input the size of the OG permutation
//     printf("Enter the size of the permutation: ");
//     scanf("%d", &size);

//     if (size <= 0 || size > MAX_N)
//     {
//         printf("Error: Size must be between 1 and %d\n", MAX_N);
//         return 1;
//     }

//     int *input = (int *)malloc(size * sizeof(int));
//     if (!input)
//     {
//         printf("Failed to allocate memory for input array\n");
//         return 1;
//     }
//     clock_t start_time = clock();

//     // Input the permutation
//     printf("Enter the permutation (space-separated): ");
//     for (int i = 0; i < size; i++)
//         scanf("%d", &input[i]);

//     // translocate(input, dst, size, 1, 2, 4);
//     print_array(input, size);

//     computeMostOddCycle_BFS(input, size);

//     clock_t end_time = clock();
//     double total_time = (double)(end_time - start_time) / CLOCKS_PER_SEC;
//     printf("\nExecution time: %.3f seconds\n", total_time);

//     return 0;
// }

int main()
{
    int size;

    clock_t start_time = clock();

    printf("Enter the size of the permutation: ");
    scanf("%d", &size);

    if (size <= 0 || size > MAX_N)
    {
        printf("Error: Size must be between 1 and %d\n", MAX_N);
        return 1;
    }

    int *input = (int *)malloc(size * sizeof(int));
    if (!input)
    {
        printf("Failed to allocate memory for input array\n");
        return 1;
    }

    int *dst = (int *)malloc(size * sizeof(int));

    printf("Enter the permutation (space-separated): ");
    for (int i = 0; i < size; i++)
        scanf("%d", &input[i]);

    // adjacency: index 1..2n, each vertex has degree 2
    int (*adj)[2] = calloc(2 * size + 1, sizeof *adj);
    int *deg = calloc(2 * size + 1, sizeof *deg);
    if (!adj || !deg)
    {
        fprintf(stderr, "OOM\n");
        return 1;
    }

    // printf("Building edges...\n");
    build_black_edges_from_perm(adj, deg, input, size); // e.g., 2--11, 12--9, ...
    build_gray_edges_identity(adj, deg, size);          // e.g., 2--3, 4--5, ..., 14--1

    printf("\nCounting cycles...\n");
    count_cycles_colored(adj, size);

    translocate(input, dst, size, 1, 2, 4);

    free(adj);
    free(deg);

    clock_t end_time = clock();
    double total_time = (double)(end_time - start_time) / CLOCKS_PER_SEC;
    printf("\nExecution time: %.3f seconds\n", total_time);

    return 0;
}
