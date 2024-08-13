#include "../include/recenter.hpp"
#include "../include/myprompt.hpp"
#include <algorithm>
#include <cmath>
#include <memory>
#include <mpi.h>
using namespace std;

// Calculate the norm of the coordinate 1x3 vector.
#define NORM( vec3ptr )                                                             \
    sqrt( *( vec3ptr ) * *( vec3ptr ) + *( ( vec3ptr ) + 1 ) * *( ( vec3ptr ) + 1 ) \
          + *( ( vec3ptr ) + 2 ) * *( ( vec3ptr ) + 2 ) )

void recenter::run()
{
    // TODO:
}

/**
 * @brief Calculate the center of mass in specified range.
 *
 * @param mass masses of particles
 * @param coordinates coordinates of particles
 * @param partNum particle number
 * @param rangeSize enclose radius of the chosen range
 * @return uniqure_ptr to the gotton center of mass
 */
auto recenter::center_of_mass( const double* mass, const double* coordinates,
                               const unsigned int& partNum,
                               const double        rangeSize ) -> unique_ptr< double[] >
{
    auto   com( make_unique< double[] >( 3 ) );
    double massSum           = 0;
    double coordMassSum[ 3 ] = { 0, 0, 0 };

    static unsigned int i = 0;
    static unsigned int j = 0;
    for ( i = 0; i < partNum; ++i )
    {
        if ( NORM( coordinates + 3 * i ) < rangeSize )
        {
            massSum += mass[ i ];
            for ( j = 0; j < 3; ++j )
            {
                coordMassSum[ j ] += mass[ i ] * coordinates[ i * 3 + j ];
            }
        }
    }

    MPI_Allreduce( MPI_IN_PLACE, &massSum, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD );
    MPI_Allreduce( MPI_IN_PLACE, coordMassSum, 3, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD );

    if ( massSum != 0 )
    {
        for ( j = 0; j < 3; ++j )
        {
            com[ j ] = coordMassSum[ j ] / massSum;
        }
    }
    else
    {
        int rank = 0;
        MPI_Comm_rank( MPI_COMM_WORLD, &rank );
        MPI_WARN( rank, "Get an Mtot=0 in calculation of CoM, return 0 only." );
    }
    return com;
}

/**
 * @brief Calculate the position of the most bound particle.
 *
 * @param potential potential of particles
 * @param coordinates coordinates of particles
 * @param partNum particle number
 * @return the coordinates of the most bound particle
 */
auto recenter::most_bound_particle( const double* potential, const double* coordinates,
                                    const unsigned int& partNum ) -> std::unique_ptr< double[] >
{
    auto   minPotPosition( make_unique< double[] >( 3 ) );
    int    minLocateId = min_element( potential, potential + partNum ) - potential;
    double min         = potential[ minLocateId ];  // local min
    // global min after reduce
    MPI_Allreduce( MPI_IN_PLACE, &min, 1, MPI_DOUBLE, MPI_MIN, MPI_COMM_WORLD );

    // Get the rank number of the minimal potential
    int minLocateRank = 0;
    MPI_Comm_rank( MPI_COMM_WORLD, &minLocateRank );
    if ( min != potential[ minLocateId ] )  // local min!=global min
    {
        minLocateRank = 0;
    }
    else
    {
        for ( int i = 0; i < 3; ++i )
        {
            minPotPosition[ i ] = coordinates[ 3 * minLocateId + i ];
        }
    }
    MPI_Allreduce( MPI_IN_PLACE, &minLocateRank, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD );
    MPI_Bcast( minPotPosition.get(), 3, MPI_DOUBLE, minLocateRank, MPI_COMM_WORLD );

    return minPotPosition;
}
