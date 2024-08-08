#include "../include/galotfa.hpp"
/**
 * @brief api for n body simulation, without sub-grid physics parameters and redshifts.
 *
 * @param currentTime simulation time, in simulation unit.
 * @param particleNumber number of particles in this mpi.
 * @param particleIDs particle ids.
 * @param particleTypes particle type ids.
 * @param masses particle masses, in simulation units.
 * @param coordiantes coordinates in simulation unit.
 * @param velocities velocities in simulation units.
 * @return
 */
extern "C" auto OnTheFly_Analysis_Nbody( double currentTime, int particleNumber,
                                         const int* particleIDs, const int* particleTypes,
                                         const double* masses, const double* coordiantes,
                                         const double* velocities ) -> int
{
    ( void )currentTime;
    ( void )particleNumber;
    ( void )particleIDs;
    ( void )particleTypes;
    ( void )masses;
    ( void )coordiantes;
    ( void )velocities;
    return 0;
}
