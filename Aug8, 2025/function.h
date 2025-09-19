#ifndef FUNCTION_H
#define FUNCTION_H

// Declaration of function

// ---------- Helper function-----------------
void print_array(int arr[], int n);
void compute_inverse(int pi[], int pi_inv[], int n);
void shift_permutation_by_one(int *src, int *dest, int n);
void translocate(const int *src, int *dst, int n, int i, int j, int k);
long long factorial(int n);
void initialize_identity_permutation(int *pid, int n);

//===============Calculating values (max_len, max_cycle)=============
int longest_increasing_consecutive_values(int *arr, int n);
int longest_increasing_subsequence(int *arr, int n);
int computeMaxLen(int pi[], int n);    // Using longest_increasing_consecutive_values
int computeMaxLen_v2(int pi[], int n); // Using longest_increasing_subsequence

// ================= Read distances array from files ============
int *load_D_from_file(int n, long long *size_out);

// ----------- Ranking function ----------------
void swap(int *a, int *b);
int rank1(int n, int pi[], int pi_inv[]);
int rank_safe(int n, const int src[], int *inv_buf);
void unrank1(int n, int r, int pi[]);

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

// ============== Find best neighbor pairs =============
// Structure to hold both metrics for a neighbor
typedef struct
{
    int max_len;
    int max_cycles;
    int best_i, best_j, best_k; // Optional: track which translocation was best
} BestNeighborMetrics;

BestNeighborMetrics findBestNeighborMetrics(int *perm, int n);

// ============== Print bad permutation ==================
void printBadTranslocationFromIdentityCombined_Level1(int n, int *distance_array);
void printBadTranslocationFromIdentityCombined_Level2(int n, int *distance_array);

#endif /* FUNCTION_H */
