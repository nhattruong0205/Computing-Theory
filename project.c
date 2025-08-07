#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#define N 10
#define MAX_QUEUE 10000000
#define MAX_PERMS 3628800 // 6! permutations

// Queue for BFS
int queue[MAX_QUEUE][N];
int front = 0, rear = 0;

// Visited array and parent mapping
int visited[MAX_PERMS] = {0};
int parent[MAX_PERMS];

// Factorial lookup for permutation indexing
int factorial[] = {1, 1, 2, 6, 24, 120, 720};

// Encode a permutation to a unique index using Lehmer code
int perm_to_index(int perm[N])
{
    int idx = 0;
    int used[N] = {0};
    for (int i = 0; i < N; i++)
    {
        int count = 0;
        for (int j = 0; j < perm[i] - 1; j++)
            if (!used[j])
                count++;
        idx += count * factorial[N - 1 - i];
        used[perm[i] - 1] = 1;
    }
    return idx;
}

// Copy a permutation
void copy_perm(int dest[N], int src[N])
{
    for (int i = 0; i < N; i++)
        dest[i] = src[i];
}

// Print a permutation
void print_perm(int perm[N])
{
    printf("(");
    for (int i = 0; i < N; i++)
    {
        printf("%d", perm[i]);
        if (i < N - 1)
            printf(", ");
    }
    printf(")\n");
}

// Add to BFS queue
void enqueue(int perm[N])
{
    for (int i = 0; i < N; i++)
        queue[rear][i] = perm[i];
    rear++;
}

// Remove from BFS queue
void dequeue(int perm[N])
{
    for (int i = 0; i < N; i++)
        perm[i] = queue[front][i];
    front++;
}

// Perform BFS and store optimal path
void bfs(int start[N], int target[N])
{
    enqueue(start);
    int start_idx = perm_to_index(start);
    visited[start_idx] = 1;
    parent[start_idx] = -1;

    int current[N];

    while (front < rear)
    {
        dequeue(current);
        int current_idx = perm_to_index(current);

        if (memcmp(current, target, sizeof(int) * N) == 0)
            break;

        for (int i = 0; i < N; i++)
        {
            for (int j = i; j < N; j++)
            {
                int interval_len = j - i + 1;
                int interval[interval_len];
                for (int x = 0; x < interval_len; x++)
                    interval[x] = current[i + x];

                int rest[N - interval_len];
                int r = 0;
                for (int x = 0; x < N; x++)
                {
                    if (x < i || x > j)
                        rest[r++] = current[x];
                }

                for (int k = 0; k <= r; k++)
                {
                    if (k == i)
                        continue;

                    int new_perm[N], pos = 0;
                    for (int x = 0; x < k; x++)
                        new_perm[pos++] = rest[x];
                    for (int x = 0; x < interval_len; x++)
                        new_perm[pos++] = interval[x];
                    for (int x = k; x < r; x++)
                        new_perm[pos++] = rest[x];

                    int new_idx = perm_to_index(new_perm);
                    if (!visited[new_idx])
                    {
                        visited[new_idx] = 1;
                        parent[new_idx] = current_idx;
                        enqueue(new_perm);
                    }
                }
            }
        }
    }

    // Reconstruct and print path
    int path[MAX_PERMS][N];
    int path_len = 0;
    int idx = perm_to_index(target);
    while (idx != -1)
    {
        for (int i = 0; i < N; i++)
            path[path_len][i] = queue[idx][i];
        idx = parent[idx];
        path_len++;
    }

    printf("\nOptimal path of translocations:\n");
    for (int i = path_len - 1; i >= 0; i--)
    {
        printf("Step %d: ", path_len - 1 - i);
        print_perm(path[i]);
    }
    printf("\nMinimum translocations: %d\n", path_len - 1);
}

int main()
{
    int start[N] = {1, 2, 3, 4, 5, 6, 7, 8, 9, 10};
    int target[N] = {10, 9, 8, 7, 6, 5, 4, 3, 2, 1};

    bfs(start, target);
    return 0;
}
