/**
 * @file test_recenter.cpp
 * @brief Test the recenter APIs
 */
#include "../include/recenter.hpp"
#ifdef DEBUG
#include "../include/myprompt.hpp"
#endif
#include <cassert>
#include <cstring>
#include <mpi.h>
using namespace std;
#define THRESHOLD 1e-6  // the equal threshold of floating numbers

int main( int argc, char* argv[] )
{
    // check whether the double[3] vecctors are equal within the threshold
    auto neq = []( double vec1[ 3 ], double vec2[ 3 ] ) -> bool {
        int eq = 0;
        for ( int i = 0; i < 3; ++i )
            eq += ( int )( abs( vec1[ i ] - vec2[ i ] ) >= THRESHOLD );
        return eq;
    };

    // NOTE: test in rank 4 mpi process program: 10 coordinate in each rank
    int rank = -1, size = 0;
    MPI_Init( &argc, &argv );
    MPI_Comm_rank( MPI_COMM_WORLD, &rank );
    MPI_Comm_size( MPI_COMM_WORLD, &size );
    assert( size = 4 );  // check the mpi size

    // TEST: center of mass calculation
    // TEST: trival case, all==0
    double coordinate[ 120 ] = { 0 };
    double mass[ 40 ]        = { 0 };
    for ( int i = 0; i < 4; ++i )
    {
        if ( rank == i )
        {
            for ( int j = 0; j < 10; ++j )
            {
                for ( int k = 0; k < 3; ++k )
                    coordinate[ 3 * ( i * 10 + j ) + k ] = 0;
                mass[ i * 10 + j ] = 1.25;
            }
        }
    }
    double expected1[ 3 ] = { 0, 0, 0 };
    auto   res1 =
        recenter::center_of_mass( mass + 10 * rank, coordinate + 3 * ( 10 * rank ), 10, 100 );
    MPI_Barrier( MPI_COMM_WORLD );
    assert( !neq( expected1, res1.get() ) );

    // TEST: trival case, all==3.14
    memset( coordinate, 0, sizeof( double ) * 120 );
    for ( int i = 0; i < 4; ++i )
    {
        if ( rank == i )
        {
            for ( int j = 0; j < 10; ++j )
            {
                for ( int k = 0; k < 3; ++k )
                    coordinate[ 3 * ( i * 10 + j ) + k ] = 3.14;
                mass[ i * 10 + j ] = 1.25;
            }
        }
    }

    double expected2[ 3 ] = { 3.14, 3.14, 3.14 };
    auto   res2 =
        recenter::center_of_mass( mass + 10 * rank, coordinate + 3 * ( 10 * rank ), 10, 100 );
    MPI_Barrier( MPI_COMM_WORLD );
    assert( !neq( expected2, res2.get() ) );

    // TEST: equal mass case
    double coordinate2[ 120 ] = {
        0.58801452, 0.69910875, 0.18815196, 0.04380856, 0.20501895, 0.10606287, 0.72724014,
        0.67940052, 0.4738457,  0.44829582, 0.01910695, 0.75259834, 0.60244854, 0.96177758,
        0.66436865, 0.60662962, 0.44915131, 0.22535416, 0.6701743,  0.73576659, 0.25799564,
        0.09554215, 0.96090974, 0.25176729, 0.28216512, 0.76825393, 0.7979234,  0.5440372,
        0.38270763, 0.38165095, 0.28582739, 0.74026815, 0.23898683, 0.4377217,  0.8835387,
        0.28928114, 0.78450686, 0.75895366, 0.41778538, 0.22576877, 0.42009814, 0.06436369,
        0.59643269, 0.83732372, 0.89248639, 0.20052744, 0.50239523, 0.89538184, 0.25592093,
        0.86723234, 0.01648793, 0.55249695, 0.52790539, 0.92335039, 0.24594844, 0.06401838,
        0.9021047,  0.8740398,  0.16366729, 0.99974131, 0.34680397, 0.31287816, 0.84710402,
        0.8802311,  0.67655865, 0.05367515, 0.55921377, 0.69451294, 0.8241973,  0.31142866,
        0.50523054, 0.84900379, 0.29351563, 0.67711955, 0.4209064,  0.68171271, 0.22122799,
        0.5489977,  0.84884672, 0.7365669,  0.49962259, 0.37966499, 0.78752081, 0.16886931,
        0.58635861, 0.43121067, 0.06191019, 0.28945477, 0.7341454,  0.28865545, 0.39039811,
        0.63561732, 0.83114886, 0.319421,   0.15922479, 0.71166422, 0.87270864, 0.59315637,
        0.69471288, 0.17323332, 0.53173259, 0.87095862, 0.84109027, 0.97205554, 0.78225721,
        0.19703051, 0.61062607, 0.47885551, 0.616637,   0.13993324, 0.41123582, 0.77763034,
        0.93972552, 0.10457941, 0.9384822,  0.79738717, 0.33080272, 0.31178575, 0.29015382,
        0.17388959
    };
    /* python code:
    import numpy as np
    np.random.seed(2024)
    mockCoordinate = np.random.rand(120).reshape(40, 3)
    mockCoordinate.reshape(120)
    np.mean(mockCoordinate, axis=0)
    */
    memset( coordinate, 0, sizeof( double ) * 120 );
    for ( int i = 0; i < 4; ++i )
    {
        if ( rank == i )
        {
            for ( int j = 0; j < 10; ++j )
            {
                for ( int k = 0; k < 3; ++k )
                    coordinate[ 3 * ( i * 10 + j ) + k ] = coordinate2[ 3 * ( i * 10 + j ) + k ];
                mass[ i * 10 + j ] = 1;
            }
        }
    }
    double expected3[ 3 ] = { 0.49207988, 0.57682968, 0.49231838 };
    auto   res3 =
        recenter::center_of_mass( mass + 10 * rank, coordinate + 3 * ( 10 * rank ), 10, 100 );
    MPI_Barrier( MPI_COMM_WORLD );
    assert( !neq( expected3, res3.get() ) );

    // TEST: varing mass case
    double coordinate3[ 120 ] = {
        0.58801452, 0.69910875, 0.18815196, 0.04380856, 0.20501895, 0.10606287, 0.72724014,
        0.67940052, 0.4738457,  0.44829582, 0.01910695, 0.75259834, 0.60244854, 0.96177758,
        0.66436865, 0.60662962, 0.44915131, 0.22535416, 0.6701743,  0.73576659, 0.25799564,
        0.09554215, 0.96090974, 0.25176729, 0.28216512, 0.76825393, 0.7979234,  0.5440372,
        0.38270763, 0.38165095, 0.28582739, 0.74026815, 0.23898683, 0.4377217,  0.8835387,
        0.28928114, 0.78450686, 0.75895366, 0.41778538, 0.22576877, 0.42009814, 0.06436369,
        0.59643269, 0.83732372, 0.89248639, 0.20052744, 0.50239523, 0.89538184, 0.25592093,
        0.86723234, 0.01648793, 0.55249695, 0.52790539, 0.92335039, 0.24594844, 0.06401838,
        0.9021047,  0.8740398,  0.16366729, 0.99974131, 0.34680397, 0.31287816, 0.84710402,
        0.8802311,  0.67655865, 0.05367515, 0.55921377, 0.69451294, 0.8241973,  0.31142866,
        0.50523054, 0.84900379, 0.29351563, 0.67711955, 0.4209064,  0.68171271, 0.22122799,
        0.5489977,  0.84884672, 0.7365669,  0.49962259, 0.37966499, 0.78752081, 0.16886931,
        0.58635861, 0.43121067, 0.06191019, 0.28945477, 0.7341454,  0.28865545, 0.39039811,
        0.63561732, 0.83114886, 0.319421,   0.15922479, 0.71166422, 0.87270864, 0.59315637,
        0.69471288, 0.17323332, 0.53173259, 0.87095862, 0.84109027, 0.97205554, 0.78225721,
        0.19703051, 0.61062607, 0.47885551, 0.616637,   0.13993324, 0.41123582, 0.77763034,
        0.93972552, 0.10457941, 0.9384822,  0.79738717, 0.33080272, 0.31178575, 0.29015382,
        0.17388959
    };
    double mass2[ 40 ] = { 0.88114956, 0.82005649, 0.73665888, 0.84696823, 0.60911562, 0.34423301,
                           0.22966899, 0.94462579, 0.29172571, 0.41004459, 0.33766428, 0.9493286,
                           0.23092893, 0.16947204, 0.49089763, 0.90475345, 0.88399969, 0.19491159,
                           0.52021797, 0.17489932, 0.89954588, 0.37819933, 0.70288128, 0.51174545,
                           0.63485106, 0.93399494, 0.21049022, 0.33620401, 0.65946718, 0.41460336,
                           0.25531344, 0.87560225, 0.40179131, 0.56136793, 0.1716831,  0.00877794,
                           0.40878872, 0.10346889, 0.52504683, 0.45986097 };
    /* python code:
    import numpy as np
    np.random.seed(2024)
    mockCoordinate = np.random.rand(120).reshape(40, 3)
    mockCoordinate.reshape(120)
    mass = np.random.rand(40)
    np.sum(mockCoordinate*np.repeat(mass, 3).reshape(40, 3), axis=0) / np.sum(mass)
    */
    memset( coordinate, 0, sizeof( double ) * 120 );
    for ( int i = 0; i < 4; ++i )
    {
        if ( rank == i )
        {
            for ( int j = 0; j < 10; ++j )
            {
                for ( int k = 0; k < 3; ++k )
                    coordinate[ 3 * ( i * 10 + j ) + k ] = coordinate3[ 3 * ( i * 10 + j ) + k ];
                mass[ i * 10 + j ] = mass2[ i * 10 + j ];
            }
        }
    }
    double expected4[ 3 ] = { 0.44997152, 0.55214703, 0.49108783 };
    auto   res4 =
        recenter::center_of_mass( mass + 10 * rank, coordinate + 3 * ( 10 * rank ), 10, 100 );
    MPI_Barrier( MPI_COMM_WORLD );
    assert( !neq( expected4, res4.get() ) );

    MPI_Finalize();
    return 0;
}
