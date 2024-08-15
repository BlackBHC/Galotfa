/**
 * @file monitor.cpp
 * @brief The organizer of other components to work together.
 */

#include "../include/monitor.hpp"
#include <mpi.h>
#ifdef DEBUG
#include "../include/myprompt.hpp"
#endif
#include <memory>
using namespace std;

namespace otf {

monitor::monitor( const std::string_view& tomlParaFile )
    : mpiRank( -1 ), stepCounter( 0 ), mpiInitialzedByMonitor( false ),
      para( make_unique< runtime_para >( tomlParaFile ) )
{
    // check whether in mpi environment, if not, call MPI_Init
    int initialzed;
    MPI_Initialized( &initialzed );
    if ( not initialzed )
    {
        mpiInitialzedByMonitor = true;
        MPI_Init( NULL, NULL );
    }

    // get the rank id
    MPI_Comm_rank( MPI_COMM_WORLD, &mpiRank );

    // read in the parameters
    MPI_INFO( mpiRank, "Read in parameter from %s", tomlParaFile.data() );
}

monitor::~monitor()
{
    if ( mpiInitialzedByMonitor )
    {
        MPI_Finalize();
    }
}
void monitor::one_analysis_api( const unsigned int* id, const unsigned int* partType,
                                const double* mass, const double* coordinate,
                                const double* velocity )
{
    if ( not para->enableOtf )
        return;

    ( void )id;
    ( void )partType;
    ( void )mass;
    ( void )coordinate;
    ( void )velocity;

    // first: extract the data for orbital log, and the data for each component

    // second: log the orbits

    // third: analyze each component

    // increase the synchronized step counter
    stepCounter++;
}

}  // namespace otf
