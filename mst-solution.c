#include <limits.h>
#include "mst-solution.h"
#include "mst-kruskal.c"

void order_edge(edge *e) {
  if (e->a > e->b) {
    int t = e->a;
    e->a = e->b;
    e->b = t;
  }
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

  MPI_Op minop;
  MPI_Datatype etype;

  // defining a custom mpi datatype holding edges
  MPI_Type_contiguous(3, MPI_INT, &etype);
  MPI_Type_commit(&etype);

  // defining a custom min reduce operation
  MPI_Op_create(min_edge, 1, &minop);

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

    edge *edges = malloc(M * sizeof(edge));
    int count = 0;
    msf tree;
    edge *e;

    for (int b = 1; b < N; b++)
      for (int a = 0; a < b; a++)
        if (adj[a * N + b] > 0) {
          e = &edges[count++];
          e->a = a;
          e->b = b;
          e->w = adj[a * N + b];
        }

    qsort(edges, count, sizeof(edge), cmp_edges);
    kruskal(edges, count, N, &tree);

    // one drawback of reusing methods useful for parallel kruskal
    // is that here we store the tree in memory, when it's not useful
    // since we only have one processor and we would do just as good
    // if we just printed edges as we found them
    for (int i = 0; i < N - 1; i++) {
      e = &(tree.edges[i]);
      printf("%i %i\n", e->a, e->b);
    }

    free(edges);
    free(tree.edges);

  } else if (strcmp(algoName, "prim-par") == 0) { // Parallel Prim's algorithm
    // BEGIN IMPLEMENTATION HERE


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

      for (int i = 0; i < size; i++)
        if (!T[i + offset] && D[i].e.w > 0 && cmp_edges(&D[i].e, &choice) < 0) {
          choice = D[i].e;
        }

      MPI_Allreduce(&choice, &pick, 1, etype, minop, MPI_COMM_WORLD);

      if (procRank == 0)
        printf("%i %i\n", pick.a, pick.b);


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

    int size = ceil((float)N / (float)numProcs);
    int offset = size * procRank;
    int rank = procRank;
    int procs = numProcs;
    int step = 0;
    msf forest;

    if (procRank == numProcs - 1) {
      size = N - offset;
    }



    // FIRST STEP: LOCAL KRUSKAL
    edge *edges = malloc(M * sizeof(edge));
    int count = 0;

    msf tree;
    edge *e;

    for (int a = 0; a < size; a++)
      for (int b = a + offset + 1; b < N; b++)
        if (adj[a * N + b] > 0) {
          e = &edges[count++];
          e->a = a + offset;
          e->b = b;
          e->w = adj[a * N + b];
        }

    qsort(edges, count, sizeof(edge), cmp_edges);
    kruskal(edges, count, N, &forest);



    // SECOND STEP: P2P MERGE
    // reserve memory for received data
    msf incoming;
    incoming.edges = malloc(N * sizeof(edge));

    while (procs > 1) {
      // odd rank, send data to associated processor with even rank
      // (associated even processor always exists)
      if (rank & 1) {
        int target = (rank ^ 1) << step;
        MPI_Send(&forest.size, 1, MPI_INT, target, 0, MPI_COMM_WORLD);
        MPI_Send(forest.edges, forest.size, etype, target, 0, MPI_COMM_WORLD);
        break;
      }

      // if rank is even, check if we should actually expect data from some odd processor
      else if(rank + 1 < procs) {
        int target = (rank | 1) << step;
        MPI_Recv(&incoming.size, 1, MPI_INT, target, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        MPI_Recv(incoming.edges, incoming.size, etype, target, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        merge(&forest, &incoming, N);
      }

      procs = (procs + 1) >> 1;
      rank >>= 1;
      step++;
    }

    // processor 0 now has all the data
    if (procRank == 0)
      for (int i = 0; i < N - 1; i++) {
        e = &forest.edges[i];
        printf("%i %i\n", e->a, e->b);
      }

    free(edges);
    free(forest.edges);

  } else { // Invalid algorithm name
    if (procRank == 0) {
      printf("ERROR: Invalid algorithm name: %s.\n", algoName);
      MPI_Abort(MPI_COMM_WORLD, 1);
    }
  }
}
