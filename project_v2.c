#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#define MAX_SIZE 10000 // How to fix this max size !!!!!
#define MAX_N 10       // max length of a permutation
#define INF 99999      // Large number represent infinity

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

// ------------------------End of help function-----

//------------- Begin queue functions----------------
// Defining the Queue structure
typedef struct
{
    int items[MAX_SIZE][MAX_N];
    int front;
    int rear;
    int n;
} Queue;

// Function to initialize the queue
void initializeQueue(Queue *q, int n)
{
    q->front = 0;
    q->rear = 0;
    q->n = n;
}

// Function to check if the queue is full
// Function to check if the queue is full
bool isFull(Queue *q)
{
    return q->rear >= MAX_SIZE; // Check if the rear index has reached the maximum size
}

// Function to check if the queue is empty
bool isEmpty(Queue *q)
{
    return q->front == q->rear; // Check if the front and rear indices are the same
}

void enqueue(Queue *q, int perm[])
{
    if (isFull(q))
    {
        printf("Queue is full\n");
        return;
    }
    for (int i = 0; i < q->n; i++)
    {
        q->items[q->rear][i] = perm[i];
    }
    q->rear++;
}

void pop(Queue *q, int dest[])
{
    if (isEmpty(q))
    {
        printf("Queue is empty\n");
        return;
    }
    for (int i = 0; i < q->n; i++)
    {
        dest[i] = q->items[q->front][i];
    }
    q->front++;
}

// Function to print the current queue
void printQueue(Queue *q)
{
    if (isEmpty(q))
    {
        printf("Queue is empty\n");
        return;
    }

    printf("Current Queue:\n");
    for (int i = q->front; i < q->rear; i++)
    {
        print_array(q->items[i], q->n); // Print each permutation in the queue
    }
    printf("\n");
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

void print_D(int n, const int *D, const int pi[], long long size)
{
    for (int i = 0; i < size; i++)
    {
        initialize_identity_permutation(pi, n); // reset
        unrank1(n, i, pi);
        printf("Index %d: ", i);
        print_array(pi, n);
        printf("distance = %d\n", D[i]);
    }
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

    Queue Q;
    initializeQueue(&Q, n);

    enqueue(&Q, pi);
    // printQueue(&Q);

    int result[MAX_N];
    while (!isEmpty(&Q))
    {
        pop(&Q, result);
        // printf("Popped permutation: ");
        // print_array(result, n);
        // printf("Q: ");
        // printQueue(&Q);
        int *tmp = (int *)malloc(n * sizeof(int));
        int *tmp_inv = (int *)malloc(n * sizeof(int));
        for (int i = 0; i < n; ++i)
            for (int j = i + 1; j < n; ++j)
                for (int k = j; k < n; ++k)
                {
                    // // ----- Added-----

                    // memcpy(tmp, pi, n * sizeof(int)); // Copy original permutation
                    // translocate(tmp, n, i, j, k);

                    // // ----- Added-----

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
                    compute_inverse(result, tmp_inv, n);
                    int rank_cur = rank_safe(n, result, tmp_inv);

                    // Neighbour rank
                    compute_inverse(tmp, tmp_inv, n);
                    int rank_tmp = rank_safe(n, tmp, tmp_inv);
                    // printf("Rank temp: %d\n", rank_tmp);

                    if (D[rank_tmp] > D[rank_cur] + 1)
                    {
                        D[rank_tmp] = D[rank_cur] + 1;

                        enqueue(&Q, tmp);

                        printQueue(&Q);
                    }
                }
    }
    // print_D(n, D, pi, size);
    return D;
}

// Run ComputeTDistanceFromIdentity
int main()
{
    int n = 4;

    // int result = rank1(n, pi, pi_inv);
    // printf("Rank: %d\n", result);
    int *pi = (int *)malloc(n * sizeof(int));
    int *pi_inv = (int *)malloc(n * sizeof(int));
    long long size = factorial(n);              // Calculate n!
    int *D = (int *)malloc(size * sizeof(int)); // Allocate memory for D

    D = ComputeTDistanceFromIdentity(n);
    // printf("Distance: %d\n", D[5]);
    print_D(n, D, pi, size);
}

// ------------ End of main functions-------------
//---- End of test using function translocation
// int main()
// {
//     int n = 3;
//     int *pi = (int *)malloc(n * sizeof(int));
//     long long size = factorial(n); // Calculate n!
//     int translocated[MAX_N];

//     initialize_identity_permutation(pi, n);

//     for (int i = 0; i < n; i++)
//     {
//         for (int j = i + 1; j < n; j++)
//         {
//             for (int k = j; k < n; k++)
//             {                                              // Because starting from 0
//                 memcpy(translocated, pi, n * sizeof(int)); // Copy original permutation
//                 translocate(translocated, n, i, j, k);
//                 printf("Translocation (%d, %d, %d) -> ", i, j, k); // Swap interval p[j] to p[k] with p[i] to p[j-1], p[0] to p[j-1] and p[k+1] to p[n-1] remain the same
//                 print_array(translocated, n);
//             }
//         }
//     }

//     free(pi);
//     return 0;
// }
// ---- End of test using function translocation

// ----- Main
// int main()
// {
//     int n = 5;
//     // int perm[] = {1, 2, 3};
//     // generate_adjacent_translocations(perm, n);

//     // int n = 10;

//     // int result = rank1(n, pi, pi_inv);
//     // printf("Rank: %d\n", result);
//     int *pi = (int *)malloc(n * sizeof(int));
//     int *pi_inv = (int *)malloc(n * sizeof(int));
//     long long size = factorial(n);              // Calculate n!
//     int *D = (int *)malloc(size * sizeof(int)); // Allocate memory for D

//     // Initialize all elements of D to infinity
//     for (long long i = 0; i < size; i++)
//     {
//         D[i] = INF;
//     }
//     initialize_identity_permutation(pi, n);
//     // print_array(pi, n);

//     compute_inverse(pi, pi_inv, n);
//     // print_array(pi_inv, n);

//     int pid = rank_safe(n, pi, pi_inv);
//     D[pid] = 0;

//     Queue Q;
//     initializeQueue(&Q, n);

//     enqueue(&Q, pi);
//     printQueue(&Q);

//     int result[MAX_N];
//     while (!isEmpty(&Q))
//     {
//         pop(&Q, result);
//         printf("Popped permutation: ");
//         print_array(result, n);
//         printf("Q: ");
//         printQueue(&Q);
//         int *tmp = (int *)malloc(n * sizeof(int));
//         int *tmp_inv = (int *)malloc(n * sizeof(int));

//         for (int i = 0; i < n; ++i)
//             for (int j = i + 1; j < n; ++j)
//                 for (int k = j; k < n; ++k)
//                 {
//                     /* Build translocated permutation */
//                     int idx = 0;

//                     /* 1. Prefix: [0..i-1] */
//                     for (int x = 0; x < i; ++x)
//                         tmp[idx++] = result[x];

//                     /* 2. Block: [j..k] */
//                     for (int x = j; x <= k; ++x)
//                         tmp[idx++] = result[x];

//                     /* 3. Middle: [i..j-1] */
//                     for (int x = i; x < j; ++x)
//                         tmp[idx++] = result[x];

//                     /* 4. Suffix: [k+1..n-1] */
//                     for (int x = k + 1; x < n; ++x)
//                         tmp[idx++] = result[x];

//                     printf("temp");
//                     print_array(tmp, n);

//                     // Pop rank
//                     compute_inverse(result, tmp_inv, n);
//                     int rank_cur = rank_safe(n, result, tmp_inv);

//                     // Neighbour rank
//                     compute_inverse(tmp, tmp_inv, n);
//                     int rank_tmp = rank_safe(n, tmp, tmp_inv);
//                     printf("Rank temp: %d\n", rank_tmp);

//                     if (D[rank_tmp] > D[rank_cur] + 1)
//                     {
//                         D[rank_tmp] = D[rank_cur] + 1;

//                         enqueue(&Q, tmp);

//                         // printQueue(&Q);
//                     }
//                 }
//     }
//     // print_D(n, D, pi, size);
// }