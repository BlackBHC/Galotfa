#include "../include/galotfa.h"
#include "../include/monitor.hpp"

otf::monitor otfServer( "./galotfa.toml" );

/**
 * @brief API for n body simulation, without sub-grid physics parameters and redshifts.
 *
 * @param currentTime simulation time, in simulation unit.
 * @param particleNumber number of particles in this mpi.
 * @param particleIDs particle ids.
 * @param particleTypes particle type ids.
 * @param masses particle masses, in simulation units.
 * @param coordinates coordinates in simulation unit.
 * @param velocities velocities in simulation units.
 * @return
 */
extern "C" void OnTheFly_Analysis_Nbody( const double currentTime, const unsigned particleNumber,
                                         const int* particleIDs, const int* particleTypes,
                                         const double* masses, const double* potentials,
                                         const double* coordinates, const double* velocities )
{
    otfServer.main_analysis_api( currentTime, particleNumber, particleIDs, particleTypes, masses,
                                 potentials, coordinates, velocities );
}
