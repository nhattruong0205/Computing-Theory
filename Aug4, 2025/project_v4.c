#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

#define MAX_SIZE 4000000 // How to fix this max size !!!!!
#define MAX_N 10         // max length of a permutation
#define INF 99999        // Large number represent infinity
#define MAX_FACT 4000000 // How to fix this max size !!!!!

// Global variables
int n, D[MAX_FACT], Q[MAX_FACT];
int front = 0, rear = 0;

//------------- Begin of helper function-------------
// Helper to print a permutation
void print_array(int arr[], int n)
{
    printf("[");
    for (int i = 0; i < n; i++)
        printf("%d ", arr[i]);
    printf("]");

    printf("\n");
}

// Helper to compute inverse permutation
void compute_inverse(int pi[], int pi_inv[], int n)
{
    for (int i = 0; i < n; i++)
        pi_inv[pi[i]] = i;
}

// Helper fucntion to check the longest increasing subsequence
int longest_increasing_contiguous_subsequence(int *arr, int n)
{
    int max_len = 1, curr_len = 1;

    for (int i = 1; i < n; i++)
    {
        if (arr[i] == arr[i - 1] + 1)
        {
            curr_len++;
            if (curr_len > max_len)
                max_len = curr_len;
        }
        else
        {
            curr_len = 1;
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
    return rear >= MAX_SIZE;
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
int queueSize()
{
    return rear - front;
}

// Function to print the current queue contents
void printQueue()
{
    if (isEmpty())
    {
        printf("Queue is empty\n");
        return;
    }

    printf("Queue contents (size=%d): [", queueSize());
    for (int i = front; i < rear; i++)
    {
        printf("%d", Q[i]);
        if (i < rear - 1)
            printf(", ");
    }
    printf("]\n");
}

//------------- End of queue functions----------------

// ------------------------Rank 1-----------------

void swap(int *a, int *b)
{
    int temp = *a;
    *a = *b;
    *b = temp;
}

// rank1 function: computes the lexicographic rank of a permutation
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

// unrank1: Builds a permutation from a given rank
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

void initialize_identity_permutation(int *pid, int n)
{
    for (int i = 0; i < n; i++)
    {
        pid[i] = i; // Initialize with the identity permutation (1, 2, ..., n)
    }
}

/* reverse a[left..right] inclusive */
void print_adjacent_translocations(const int p[], int n)
{
    int tmp[n];
    for (int i = 0; i < n; ++i)
        for (int j = i + 1; j < n; ++j)
            for (int k = j; k < n; ++k)
            {
                /* Build translocated permutation */
                int idx = 0;

                /* 1. Prefix: [0..i-1] */
                for (int x = 0; x < i; ++x)
                    tmp[idx++] = p[x];

                /* 2. Block: [j..k] */
                for (int x = j; x <= k; ++x)
                    tmp[idx++] = p[x];

                /* 3. Middle: [i..j-1] */
                for (int x = i; x < j; ++x)
                    tmp[idx++] = p[x];

                /* 4. Suffix: [k+1..n-1] */
                for (int x = k + 1; x < n; ++x)
                    tmp[idx++] = p[x];

                print_array(tmp, n);
            }
}

void print_D(int n, int *D, long long size)
{
    int pi[MAX_N];
    printf("Translocation distances from identity permutation:\n");
    printf("===============================================\n");
    for (int i = 0; i < size; i++)
    {
        initialize_identity_permutation(pi, n); // reset
        unrank1(n, i, pi);
        printf("Rank %d: ", i);
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
    printf("===============================================\n");
}

int get_max_distance(int *D, long long size)
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

void translocate(int *p, int n, int i, int j, int k)
{
    int *tmp = malloc(n * sizeof(int));
    memcpy(tmp, p, n * sizeof(int));

    int idx = 0;
    for (int x = 0; x < i; ++x)
        tmp[idx++] = p[x];
    for (int x = j; x <= k; ++x)
        tmp[idx++] = p[x];
    for (int x = i; x < j; ++x)
        tmp[idx++] = p[x];
    for (int x = k + 1; x < n; ++x)
        tmp[idx++] = p[x];

    memcpy(p, tmp, n * sizeof(int));
    free(tmp);
}

int *ComputeTDistanceFromIdentity(int n)
{
    int *pi = (int *)malloc(n * sizeof(int));
    int *pi_inv = (int *)malloc(n * sizeof(int));
    long long size = factorial(n);              // Calculate n!
    int *D = (int *)malloc(size * sizeof(int)); // Allocate memory for D

    // Initialize all elements of D to infinity
    for (long long i = 0; i < size; i++)
    {
        D[i] = INF;
    }
    initialize_identity_permutation(pi, n);
    // print_array(pi, n);

    compute_inverse(pi, pi_inv, n);
    // print_array(pi_inv, n);

    int pid = rank_safe(n, pi, pi_inv);
    D[pid] = 0;

    initQueue();
    enqueue(pid);

    // printQueue(&Q);

    int result[MAX_N];
    while (!isEmpty())
    {
        int current_rank = dequeue();
        // Convert rank back to permutation
        int result[MAX_N];
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

                    // printf("temp");

                    // print_array(tmp, n);

                    // Pop rank
                    compute_inverse(tmp, tmp_inv, n);
                    int rank_tmp = rank_safe(n, tmp, tmp_inv);

                    // printf("Distance rank temp: %d", D[rank_tmp]);
                    // printf("Distance rank temp: %d", D[current_rank]);

                    if (D[rank_tmp] > D[current_rank] + 1)
                    {
                        D[rank_tmp] = D[current_rank] + 1;
                        // printf("Distance rank temp: %d", D[rank_tmp]);

                        // printf("Before enqueue: ");
                        // printQueue();
                        enqueue(rank_tmp);
                        // printf("Before enqueue: ");
                        // printQueue();
                    }
                }
    }
    // print_D(n, D, pi, size);
    free(pi);
    free(pi_inv);
    return D;
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
// Run ComputeTDistanceFromIdentity
int main()
{

    printf("Enter the value of n (permutation length): ");
    scanf("%d", &n);

    if (n <= 0 || n > MAX_N)
    {
        printf("Error: n must be between 1 and %d\n", MAX_N);
        return 1;
    }

    printf("Starting computation for n=%d\n", n);
    int *pi = (int *)malloc(n * sizeof(int));
    initialize_identity_permutation(pi, n);

    int *distance_array = ComputeTDistanceFromIdentity(n);

    long long size = factorial(n);
    printBadTranslocationFromIdentity(n, distance_array);

    // unrank1(n,3, pi);
    // print_array(pi, n);
    // printf("Distance: %d\n", distance_array[23]);
}

// int main()
// {
//     clock_t start_time = clock();

//     n = 7;
//     int *pi = (int *)malloc(n * sizeof(int));
//     initialize_identity_permutation(pi, n);

//     int *distance_array = ComputeTDistanceFromIdentity(n);

//     long long size = factorial(n);
//     // print_D(n, distance_array, size);

//     int max_dist = get_max_distance(distance_array, size);
//     printf("Maximum reachable distance = %d\n", max_dist);

//     // unrank1(n, 4, pi);
//     // print_array(pi, n);
//     // printf("Distance: %d", distance_array[4]);

//     clock_t end_time = clock();
//     double total_program_time = (double)(end_time - start_time) / CLOCKS_PER_SEC;

//     printf("\nTotal program execution time: %.3f seconds\n", total_program_time);

//     return 0;
// }
