#ifndef FUNCTION_H
#define FUNCTION_H
#include <stdbool.h>

// Declaration of function

// ---------- Helper function-----------------
void print_array(int arr[], int n);
void compute_inverse(int pi[], int pi_inv[], int n);
void shift_permutation_by_one(int *src, int *dest, int n);
void translocate(const int *src, int *dst, int n, int i, int j, int k);
long long factorial(int n);
void initialize_identity_permutation(int *pid, int n);

//===============Calculating max len=============
int longest_increasing_consecutive_values(int *arr, int n);
int longest_increasing_subsequence(int *arr, int n);
int computeMaxLen(int pi[], int n);    // Using longest_increasing_consecutive_values
int computeMaxLen_v2(int pi[], int n); // Using longest_increasing_subsequence

// ================= Read distances array from files ============
int *load_D_from_file_LehmerAscendingRadix(int n, long long *size_out);
int *load_D_from_file_lex(int n, long long *size_out);
int *load_D_from_file_Lehmer(int n, long long *size_out);

// ----------- Ranking function ----------------
void swap(int *a, int *b);
int rank_lex(int pi[], int n);
void unrank_lex(int n, int r, int pi[]);
int rank1(int n, int pi[], int pi_inv[]);
int rank_safe(int n, const int src[], int *inv_buf);
void unrank1(int n, int r, int pi[]);
int rank2(int n, int pi[], int pi_inv[]);
int rank2_safe(int n, const int src[], int *inv_buf);
void unrank2(int n, int r, int pi[]);

// ============== Computing cycles functions ==========
int *creatingBreakpointGraph(int arr[], int size);
static inline int L(int x);
static inline int R(int x);
static void add_edge(int (*adj)[2], int *deg, int u, int v, const char *label);
static void build_black_edges_from_perm(int (*adj)[2], int *deg, int *pi, int n);
static void build_gray_edges_identity(int (*adj)[2], int *deg, int n);
static void count_cycles_colored(int (*adj)[2], int n);
int count_odd_cycles(int (*adj)[2], int n);
void print_cycles(int (*adj)[2], int n);
void reset_graph(int n, int adj[][2], int *deg);
int computeMostOddCycle(int *perm, int n);

// ============== Find best neighbor pairs =============
// Structure to hold both metrics for a neighbor
typedef struct
{
    int max_len;
    int max_cycles;
    int best_i, best_j, best_k; // Optional: track which translocation was best
} BestNeighborMetrics;

BestNeighborMetrics findBestNeighborMetrics(int *perm, int n);

// ================= Queue functions ===========================
void initQueue();
bool isFull();
bool isEmpty();
void enqueue(int rank);
int dequeue();
long long queueSize();
void printQueue();

// =============== Dealing with memory ========================
bool allocate_memory(int n);
void free_memory();

// =============== Saving files function ======================
void save_D_to_file(const char *filename, int *D, long long size);

//================== Compute distance array ====================
int *ComputeTDistanceFromIdentity_lex(int n);
int *ComputeTDistanceFromIdentity_LehmerAscendingRadix(int n);
int *ComputeTDistanceFromIdentity_Lehmer(int n);

// ================== Computing PAs =========================
// Computing T(n,d) PA - an array A of permutation on [1..n] with dt(A) >= d.
// Using Lehmer Ascending Radix ranking
int distance_between_2_permutations_LehmerAscendingRadix(int n, int *pi, int *sigma, int *D);
long long T_LehmerAscendingRadix(int n, int d, int *D); // Using Lehmer Ascending Radix ranking

// Using  lex ranking
int distance_between_2_permutations_lex(int n, int *pi, int *sigma, int *D);
long long T_lex(int n, int d, int *D); // Using lexigraphical ranking

// Using Lehmer ranking
int distance_between_2_permutations_Lehmer(int n, int *pi, int *sigma, int *D);
long long T_Lehmer(int n, int d, int *D); // Using Lehmer ranking

// ============== Print bad permutation ==================
// Max len based on max consecutive values
void printBadTranslocationMaxConsecutiveValue_Level1(int n, int *distance_array);
void printBadTranslocationMaxConsecutiveValue_Level2(int n, int *distance_array);

// Max odd cycle
void printBadTranslocationOddCycle_Level1(int n, int *distance_array);
void printBadTranslocationOddCycle_Level2(int n, int *distance_array);

// Combine longest_increasing_consecutive_values and max_odd_cycle
void printBadTranslocationFromIdentityCombined_Level1(int n, int *distance_array);
void printBadTranslocationFromIdentityCombined_Level2(int n, int *distance_array);

// ============== Check if bad permutation =====================================
#endif /* FUNCTION_H */
