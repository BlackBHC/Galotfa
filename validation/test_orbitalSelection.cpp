/**
 * @file test_orbitalSelection.cpp
 * @brief Test the orbital particle selection.
 */

#include "../include/myprompt.hpp"
#include "../include/selector.hpp"
#include <cmath>
#include <memory>
#include <mpi.h>
using namespace std;
using namespace otf;

int main( int argc, char* argv[] )
{
    // NOTE: test in rank 4 mpi process program: 10 coordinate in each rank
    int rank = -1;
    int size = 0;
    MPI_Init( &argc, &argv );
    MPI_Comm_rank( MPI_COMM_WORLD, &rank );
    MPI_Comm_size( MPI_COMM_WORLD, &size );
    assert( size == 4 );  // check the mpi size

    unsigned int mockTypes[ 40 ];
    for ( auto& mockType : mockTypes )
    {
        mockType = 2;
    }

    unsigned int mockIDs[ 40 ];
    for ( auto i = 0U; i < 40; ++i )
    {
        mockIDs[ i ] = i + 1;
    }

    double mockMass[ 40 ];
    for ( auto i = 0U; i < 40; ++i )
    {
        mockMass[ i ] = ( double )i + 0.1;
    }

    double mockPos[ 120 ];
    for ( auto i = 0U; i < 120; ++i )
    {
        mockPos[ i ] = ( double )i + 0.2;
    }

    double mockVel[ 120 ];
    for ( auto i = 0U; i < 120; ++i )
    {
        mockVel[ i ] = pow( -1, i ) * ( ( double )i + 0.3 );
    }

    unique_ptr< runtime_para > para =
        make_unique< runtime_para >( "../validation/test_orbit_random.toml" );
    orbit_selector orbitSelector( para );
    auto           getData =
        orbitSelector.select( 10, mockIDs + 10 * rank, mockTypes + 10 * rank, mockMass + 10 * rank,
                              mockPos + 3 * 10 * rank, mockVel + 3 * 10 * rank );
    assert( getData->count == 9 );
    MPI_Finalize();
    return 0;
}
