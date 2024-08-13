#include "../include/barinfo.hpp"
#include <cmath>
#include <mpi.h>
using namespace std;

/**
 * @brief Calculate the A0 Fourier coefficient
 *
 * @param partNum particle number
 * @param mass masses of partciles
 * @return the A0 value
 */
auto bar_info::A0( const unsigned int partNum, const double* mass ) -> double
{
    double A0sum = 0;
    for ( auto i = 0U; i < partNum; ++i )
    {
        A0sum += mass[ i ];
    }
    MPI_Allreduce( MPI_IN_PLACE, &A0sum, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD );
    return A0sum;
}

/**
 * @brief Calculate the A2 Fourier coefficient
 *
 * @param partNum particle number
 * @param mass masses of partciles
 * @param phi azimuthal angle of the particles
 * @return the A2 value
 */
auto bar_info::A2( const unsigned int partNum, const double* mass, const double* phi ) -> double
{
    double A2sumRe = 0;  // real part
    double A2sumIm = 0;  // imaginary part
    for ( auto i = 0U; i < partNum; ++i )
    {
        A2sumRe += mass[ i ] * cos( 2 * phi[ i ] );
        A2sumIm += mass[ i ] * sin( 2 * phi[ i ] );
    }
    MPI_Allreduce( MPI_IN_PLACE, &A2sumRe, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD );
    MPI_Allreduce( MPI_IN_PLACE, &A2sumIm, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD );
    return sqrt( A2sumRe * A2sumRe + A2sumIm * A2sumIm );
}

/**
 * @brief Calculate the bar strength parameter.
 *
 * @param partNum particle number
 * @param mass masses of partciles
 * @param phi azimuthal angle of the particles
 * @return the value of bar strength
 */
auto bar_info::Sbar( const unsigned int partNum, const double* mass, const double* phi ) -> double
{
    const double A0value = A0( partNum, mass );
    const double A2value = A2( partNum, mass, phi );
    return A2value / A0value;
}

/**
 * @brief Calculate the buckling strength parameter.
 *
 * @param partNum particle number
 * @param mass masses of partciles
 * @param phi azimuthal angle of the particles
 * @param zed z coordinates of the particles
 * @return the value of buckling strength
 */
auto bar_info::Sbuckle( const unsigned int partNum, const double* mass, const double* phi,
                        const double* zed ) -> double
{
    const double A0value     = A0( partNum, mass );
    double       numeratorRe = 0;  // real part
    double       numeratorIm = 0;  // imaginary part
    for ( auto i = 0U; i < partNum; ++i )
    {
        numeratorRe += mass[ i ] * zed[ i ] * cos( 2 * phi[ i ] );
        numeratorIm += mass[ i ] * zed[ i ] * sin( 2 * phi[ i ] );
    }
    MPI_Allreduce( MPI_IN_PLACE, &numeratorRe, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD );
    MPI_Allreduce( MPI_IN_PLACE, &numeratorIm, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD );
    return sqrt( numeratorRe * numeratorRe + numeratorIm * numeratorIm ) / A0value;
}
