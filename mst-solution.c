#include <limits.h>

typedef struct edge {
  int a;
  int b;
  int w;
} edge;

typedef struct p_edge {
  int v;
  int t;
  edge e;
} p_edge;

int root(int i, int *T) {
  while (T[i] != i)
     i = T[i];
  return i;
}

void order_edge(edge *e) {
  if (e->a > e->b) {
    int t = e->a;
    e->a = e->b;
    e->b = t;
  }
}

// pretty much identical, except one is for edges and the other for pedges
int cmp_edges(const void *a, const void *b) {
  edge *u = (edge *)a;
  edge *v = (edge *)b;

  if (u->w == v->w)
    return (u->a < v->a || u->a == v->a && u->b < v->b) ? -1 : 1;
  else
    return u->w - v->w;
}

void min_edge(void *in, void *inout, int *len, MPI_Datatype *dptr) {
  edge *a = (edge *)in;
  edge *b = (edge *)inout;
  int i = *len;
  while(i--) {
    if (a->w < b->w || a->w == b->w && (a->a < b->a || a->a == b->a && a->b < b->b))
      *b = *a;
    a++;
    b++;
  }
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

    p_edge *A, *D = malloc((N - 1) * sizeof(p_edge));
    int addedi, added, count = N - 1;
    edge choice;
    int *V;

    A = D;
    for (int y = 1; y < N; y++, A++) {
        A->e.w = adj[y];
        A->t = A->e.a = 0;
        A->v = A->e.b = y;
    }

    // tree construction
    while (count--) {
      choice.w = INT_MAX;

      // TODO(maybe): heap
      for (int i = 0; i <= count; i++)
        if (D[i].e.w > 0 && cmp_edges(&D[i].e, &choice) < 0) {
          choice = D[i].e;
          added  = D[i].v;
          addedi = i;
        }

      printf("%i %i\n", choice.a, choice.b);
      D[addedi] = D[count];

      V = adj + added * N;
      for (int i = 0; i < count; i++) {
        if (V[D[i].v] > 0 && (D[i].e.w == 0 || V[D[i].v] < D[i].e.w
                                            || V[D[i].v] == D[i].e.w && added < D[i].t)) {
          D[i].t = added;
          D[i].e.w = V[D[i].v];
          if (added < D[i].v) {
            D[i].e.a = added;
            D[i].e.b = D[i].v;
          }
          else {
            D[i].e.a = D[i].v;
            D[i].e.b = added;
          }
        }
      }
    }

    // garbage
    free(D);
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
    int count = 0, k = 0, roota, rootb;
    edge *e;

    for (int i = 0; i < N; i++) {
      T[i] = i;
    }

    for (int b = 1; b < N; b++) {
      for (int a = 0; a < b; a++) {
        if (adj[a * N + b] > 0) {
          e = &edges[count++];
          e->a = a;
          e->b = b;
          e->w = adj[a * N + b];
        }
      }
    }

    qsort(edges, count, sizeof(edge), cmp_edges);

    count = N - 1;

    while (count--) {
      do {
        e = &edges[k++];
        roota = root(e->a, T);
        rootb = root(e->b, T);
      } while (roota == rootb);

      T[rootb] = roota;

      printf("%i %i\n", e->a, e->b);
    }

    free(edges);
    free(T);
  } else if (strcmp(algoName, "prim-par") == 0) { // Parallel Prim's algorithm
    // BEGIN IMPLEMENTATION HERE

    MPI_Op minop;
    MPI_Datatype etype;

    // defining a custom mpi datatype holding edges
    MPI_Type_contiguous(3, MPI_INT, &etype);
    MPI_Type_commit(&etype);

    // defining a custom min reduce operation
    MPI_Op_create(min_edge, 1, &minop);

    int size = ceil((float)N / (float)numProcs);
    int length = size * N;
    int offset = size * procRank;

    if (procRank == numProcs - 1) {
      size = N - offset;
      length = N * N - length;
    }

    p_edge *A, *D = malloc(size * sizeof(int));
    int *V, *T = calloc(N, sizeof(int));
    edge pick, choice;
    int added;
    int local_count = size;

    if (procRank == 0) {
      T[0] = 1;
    }

    for (int y = 0; y < size; y++) {
        D[y].e.w = adj[y * N];
        D[y].t = D[y].e.a = 0;
        D[y].v = y;
        D[y].e.b = y + offset;
    }

    int count = N - 1;

    while (count--) {
      choice.w = INT_MAX;

      if (local_count)
        for (int i = 0; i < size; i++)
          if (!T[i + offset] && D[i].e.w > 0 && cmp_edges(&D[i].e, &choice) < 0) {
            choice = D[i].e;
          }
      }

      MPI_Allreduce(&choice, &pick, 1, etype, minop, MPI_COMM_WORLD);

      if (procRank == 0)
        printf("%i %i\n", pick.a, pick.b);

      if (pick == choice) {
        local_count--;
      }

      added = T[pick.a] ? pick.b : pick.a;
      T[added] = 1;

      for (int i = 0; i < size; i++) {
        if (!T[i + offset] && adj[i * N + added] > 0 && (D[i].e.w == 0 || adj[i * N + added] < D[i].e.w
                           || adj[i * N + added] == D[i].e.w && added < D[i].t)) {
          D[i].t = added;
          D[i].e.w = adj[i * N + added];
          if (added < D[i].e.a) {
            D[i].e.a = added;
            D[i].e.b = D[i].v + offset;
          }
          else {
            D[i].e.a = D[i].v + offset;
            D[i].e.b = added;
          }
        }
      }
    }

    // garbage
    free(T);
    free(D);
  } else if (strcmp(algoName, "kruskal-par") == 0) { // Parallel Kruskal's algorithm
    // BEGIN IMPLEMENTATION HERE

  } else { // Invalid algorithm name
    if (procRank == 0) {
      printf("ERROR: Invalid algorithm name: %s.\n", algoName);
      MPI_Abort(MPI_COMM_WORLD, 1);
    }
  }
}
