#include <limits.h>

typedef struct edge {
  int a;
  int b;
  int w;
} edge;

void print_edge(edge *e) {
  if (e->a > e->b) {
    printf("%i %i\n", e->b, e->a);
  }
  else {
    printf("%i %i\n", e->a, e->b);
  }
}

int root(int i, int *T) {
  while (T[i] != i)
     i = T[i];
  return i;
}

int cmp_weight(const void *a, const void *b) {
  return ((edge*)a)->w - ((edge*)b)->w;
}

/** Computing the Minimum Spanning Tree of a graph
 * @param N the number of vertices in the graph
 * @param M the number of edges in the graph
 * @param adj the adjacency matrix
 * @param algoName the name of the algorithm to be executed
 */
void computeMST(
    int N,
    int M,
    int *adj,
    char *algoName)
{
  int procRank, numProcs;
  MPI_Comm_rank(MPI_COMM_WORLD, &procRank);
  MPI_Comm_size(MPI_COMM_WORLD, &numProcs);

  if (strcmp(algoName, "prim-seq") == 0) { // Sequential Prim's algorithm
    if (procRank == 0) {
      if (numProcs != 1) {
        printf("ERROR: Sequential algorithm is ran with %d MPI processes.\n", numProcs);
        MPI_Abort(MPI_COMM_WORLD, 1);
      }
    }
    // BEGIN IMPLEMENTATION HERE

    // initialization
    int  *D   = malloc(N * sizeof(int));
    edge *mst = malloc((N - 1) * sizeof(edge));
    int  *T   = calloc(N, sizeof(int));
    T[0] = -1;

    int count = N - 1;
    int min, mini, a, b;
    edge *e;

    for (int y = 0; y < N; y++) {
        D[y] = adj[y];
    }

    // tree construction
    while (count--) {
      min = INT_MAX;

      // TODO: respect lexicographic priority
      // TODO(maybe): heap
      for (int i = 0; i < N; i++) {
        if (T[i] != -1 && D[i] > 0 && D[i] < min) {
          a = i;
          min = D[i];
        }
      }

      b = T[a];
      e = &(mst[count]);

      e->a = a;
      e->b = b;
      e->w = min;

      print_edge(e);

      T[a] = -1;

      for (int i = 0; i < N; i++) {
        if (T[i] != -1 && (D[i] == 0 || adj[a * N + i] < D[i])) {
          D[i] = adj[a * N + i];
          T[i] = a;
        }
      }
    }

    // garbage
    free(D);
    free(T);
    free(mst);
  } else if (strcmp(algoName, "kruskal-seq") == 0) { // Sequential Kruskal's algorithm
    if (procRank == 0) {
      if (numProcs != 1) {
        printf("ERROR: Sequential algorithm is ran with %d MPI processes.\n", numProcs);
        MPI_Abort(MPI_COMM_WORLD, 1);
      }
    }
    // BEGIN IMPLEMENTATION HERE

    int  *T     = malloc(N * sizeof(int));
    edge *edges = malloc(M * sizeof(edge));
    edge *mst   = malloc((N - 1) * sizeof(edge));
    int count = 0, k = 0, roota, rootb;
    edge *e;

    for (int i = 0; i < N; i++) {
      T[i] = i;
    }

    // initialization
    for (int a = 0; a < N - 1; a++) {
      for (int b = a + 1; b < N; b++) {
        if (adj[a * N + b] > 0) {
          e = &edges[count++];
          e->a = a;
          e->b = b;
          e->w = adj[a * N + b];
        }
      }
    }

    qsort(edges, count, sizeof(edge), cmp_weight);

    count = N - 1;

    while (count--) {
      do {
        e = &edges[k++];
        roota = root(e->a, T);
        rootb = root(e->b, T);
      } while (roota == rootb);

      // TODO: care about lexicographic priority

      T[rootb] = roota;

      print_edge(e);
      mst[count] = *e;
    }

    // garbage
    free(edges);
    free(mst);
    free(T);
  } else if (strcmp(algoName, "prim-par") == 0) { // Parallel Prim's algorithm
    // BEGIN IMPLEMENTATION HERE

  } else if (strcmp(algoName, "kruskal-par") == 0) { // Parallel Kruskal's algorithm
    // BEGIN IMPLEMENTATION HERE

  } else { // Invalid algorithm name
    if (procRank == 0) {
      printf("ERROR: Invalid algorithm name: %s.\n", algoName);
      MPI_Abort(MPI_COMM_WORLD, 1);
    }
  }
}
