#pragma once

#ifdef SEWAS_DISTRIBUTED
#include <mpi.h>
#endif

class ExecutionContext{
public:
  ExecutionContext()
  {
  }

  ~ExecutionContext()
  {
  }

  static inline auto init(int argc, char* argv[])
  {
#if SEWAS_DISTRIBUTED
    /* Start the MPI runtime */
    int provided;
    MPI_Init_thread(&argc, &argv, MPI_THREAD_SERIALIZED, &provided);
    if (MPI_THREAD_SERIALIZED != provided)
      std::cerr << "WARNING: The thread level supported by the MPI implementation is " << provided << "\n";
#endif
  }

  static inline auto rank()
  {
    int r=0;
#if SEWAS_DISTRIBUTED
    MPI_Comm_rank(MPI_COMM_WORLD, &r);
#endif
    return r;
  }

  static inline auto world()
  {
    int w=1;
#if SEWAS_DISTRIBUTED
    MPI_Comm_size(MPI_COMM_WORLD, &w);
#endif
    return w;
  }

  static inline void barrier()
  {
#if SEWAS_DISTRIBUTED
    MPI_Barrier(MPI_COMM_WORLD);
#endif
  }

  static inline void finalize()
  {
#if SEWAS_DISTRIBUTED
  MPI_Finalize();
#endif
  }
};