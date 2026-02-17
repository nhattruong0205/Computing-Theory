#ifndef FUNCTION_H
#define FUNCTION_H
#include <stdbool.h>

// Declaration of function

// ---------- Helper function-----------------
int compare_long_long(const void *a, const void *b);
bool is_rank_in_selected(long long *selected, int size, long long rank);
bool next_permutation(int *arr, int n);

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
int *load_D_from_file(int n, long long *size_out, const char *rank_name);

// ===================Ranking functions=====================
// Rank_name
// === 1. LehmerAscendingRadix ===
// === 2. Lex ===
// === 3. Lehmer ===
// === 4. SJT ===
// === 5 ReverseColexOrder ===
void swap(int *a, int *b);
int rank1(int n, int pi[], int pi_inv[]);
int rank_safe(int n, const int src[], int *inv_buf);
void unrank1(int n, int r, int pi[]);
int rank2(int n, int pi[], int pi_inv[]);
int rank2_safe(int n, const int src[], int *inv_buf);
void unrank2(int n, int r, int pi[]);
int rankSJT(int n, int pi[]);
void unrankSJT(int n, int r, int pi[], int dir[]);
int rankReverseColexOrder(int n, int pi[]);
void unrankRReverseColexOrder(int n, int r, int pi[]);

// ============== Computing cycles functions ==========
int *creatingBreakpointGraph(int arr[], int size);
int L(int x);
int R(int x);
void add_edge(int (*adj)[2], int *deg, int u, int v, const char *label);
void build_black_edges_from_perm(int (*adj)[2], int *deg, int *pi, int n);
void build_gray_edges_identity(int (*adj)[2], int *deg, int n);
void count_cycles_colored(int (*adj)[2], int n);
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
int *ComputeTDistanceFromIdentity(int n, const char *rank_name);

// ================== Computing PAs =========================
// Computing distance between 2 permutation using different rankings (rank_name)
int distance_between_2_permutations(int n, int *pi, int *sigma, int *D, const char *rank_name);

// Computing T(n,d) PA - an array A of permutation on [1..n] with dt(A) >= d.
long long T_n_d(int n, int d, int *D, const char *rank_name); // Using Lehmer Ascending Radix ranking

bool is_rank_in_selected(long long *selected, int size, long long rank);

// Compute rank of list permutation
long long *compute_ranks(int n, int perms[][n], int num_perms);

// A in PA, Adding \pi to A and check if A + \pi is still in PA
bool can_add_to_code_incremental_d2(int n, int *pi, long long *selected, int code_size, const char *rank_name);

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

// Function to compute a (n,2)-PA using greedy construction
// int compute_n2_PA(int n, int perms[][n], int max_size, const char *rank_name);
int compute_n2_PA(int n, long long *selected, int max_size, const char *rank_name);
int compute_n2_PA_Greedy_Neighbor_Deletion(int n);

bool next_permutation(int *arr, int n);
// =========== Verify (n,2)-PA
// bool verify_n2_PA(int n, int perms[][n], int num_perms, const char *rank_name);
bool verify_n2_PA(int n, const char *rank_name);

// ============== Check if bad permutation =====================================
#endif /* FUNCTION_H */
