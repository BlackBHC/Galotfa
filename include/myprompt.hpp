/**
 * @file
 * @brief Define some utility functions for convenient print.
 */
#ifndef MY_UTILS_HEADER
#define MY_UTILS_HEADER
#include <cstdio>

// print with a line break at the end
#define print( ... )                \
    {                               \
        std::printf( __VA_ARGS__ ); \
        std::printf( "\n" );        \
    }

// similar but for fprint
#define fprint( file_ptr, ... )                \
    {                                          \
        std::fprintf( file_ptr, __VA_ARGS__ ); \
        std::fprintf( file_ptr, "\n" );        \
    }
/* // similar but used for mpi case
#define mpi_print( ... )                        \
    {                                           \
        int rank;                               \
        MPI_Comm_rank( MPI_COMM_WORLD, &rank ); \
        if ( rank == 0 )                        \
        {                                       \
            std::printf( __VA_ARGS__ );         \
            std::printf( "\n" );                \
        }                                       \
    }
#define mpi_fprint( file_prt, ... )                \
    {                                              \
        int rank;                                  \
        MPI_Comm_rank( MPI_COMM_WORLD, &rank );    \
        if ( rank == 0 )                           \
        {                                          \
            std::fprintf( file_ptr, __VA_ARGS__ ); \
            std::fprintf( file_prt, "\n" );        \
        }                                          \
    } */
#define test_reach print( "Call from file [%s] line[%d].", __FILE__, __LINE__ )
// #define mpi_test_reach mpi_print( "Call from file [%s] line[%d].", __FILE__, __LINE__ )

// for info, warning, and error.
#define INFO( ... ) print( "[INFO]: " __VA_ARGS__ );
#define WARN( ... )                                                            \
    {                                                                          \
        fprint( stderr, "Call from file [%s] line[%d].", __FILE__, __LINE__ ); \
        fprint( stderr, "[WARN]: " __VA_ARGS__ );                              \
    }
#define ERROR( ... )                                                           \
    {                                                                          \
        fprint( stderr, "Call from file [%s] line[%d].", __FILE__, __LINE__ ); \
        fprint( stderr, "[ERROR]: " __VA_ARGS__ );                             \
    }
/* #define MPI_INFO( ... ) mpi_print( "[INFO]: " __VA_ARGS__ );
#define MPI_WARN( ... )                                                            \
    {                                                                              \
        mpi_fprint( stderr, "Call from file [%s] line[%d].", __FILE__, __LINE__ ); \
        mpi_fprint( stderr, "[WARN]: " __VA_ARGS__ );                              \
    }
#define MPI_ERROR( ... )                                                           \
    {                                                                              \
        mpi_fprint( stderr, "Call from file [%s] line[%d].", __FILE__, __LINE__ ); \
        mpi_fprint( stderr, "[ERROR]: " __VA_ARGS__ );                             \
    } */
#endif
