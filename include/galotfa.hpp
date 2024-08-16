/**
 * @file galotfa.hpp
 * @brief The public C style API of the galotfa library.
 */

#ifndef GALOTFA_H_INCLUDED
#define GALOTFA_H_INCLUDED
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
extern "C" auto OnTheFly_Analysis_Nbody( const double       currentTime,
                                         const unsigned int particleNumber, const int* particleIDs,
                                         const int* particleTypes, const double* masses,
                                         const double* coordiantes,
                                         const double* velocities ) -> int;

#endif
