/**
 * @file monitor.cpp
 * @brief The organizer of other components to work together.
 */

#ifdef DEBUG
#include "../include/myprompt.hpp"
#endif
#include "../include/h5out.hpp"
#include "../include/monitor.hpp"
#include "../include/para.hpp"
#include "../include/recenter.hpp"
#include "../include/selector.hpp"
#include <memory>
#include <mpi.h>
#include <string_view>
#include <vector>
using namespace std;

/**
 * @brief The global part.
 *
 * @param para reference of the runtime parameter
 */
void print_global_part( otf::runtime_para& para )
{
    INFO( "GLOBAL PARAMETERS:" );
    if ( para.enableOtf )
    {
        INFO( "On the fly analysis is enabled." );
    }
    else
    {
        INFO( "On the fly analysis is forbidden." );
    }
    INFO( "Output to  [%s]/[%s]", para.outputDir.c_str(), para.fileName.c_str() );
    INFO( "Max iteration [%u], epsilon [%g]", para.maxIter, para.epsilon );
}

/**
 * @brief The orbital part.
 *
 * @param para reference of the runtime parameter
 */
void print_orbital_part( otf::runtime_para& para )
{
    INFO( "ORBITAL LOG PARAMETERS:" );
    if ( para.orbit->enable )
    {
        INFO( "Orbital log is enabled." );
    }
    INFO( "Log period: %d", para.orbit->period );
    INFO( "Particle types to be logged:" );
    for ( auto& id : para.orbit->sampleTypes )
    {
        INFO( "%d ", id );
    }
    INFO( "Method of id selection: " );
    switch ( para.orbit->method )
    {
    case otf::orbit::id_selection_method::RANDOM:
        INFO( "Random selection." );
        break;
    case otf::orbit::id_selection_method::TXTFILE:
        INFO( "By a text file of id list." );
        break;
    default:
        ERROR( "Get into an unexpected branch!" );
    }
    INFO( "Random selection fraction: %g.", para.orbit->fraction );
    INFO( "ID list filename : %s", para.orbit->idfile.c_str() );

    if ( para.orbit->recenter.enable )
    {
        INFO( "Enable recenter before log of orbits." );
    }
    switch ( para.orbit->recenter.method )
    {
    case otf::recenter_method::COM:
        INFO( "By conter of mass." );
        break;
    case otf::recenter_method::MBP:
        INFO( "By most bound particle." );
        break;
    default:
        ERROR( "Get into an unexpected branch!" );
    }
    INFO( "Recenter radius: %g", para.orbit->recenter.radius );
    INFO( "Initial guess of the recenter:" );
    for ( auto& coord : para.orbit->recenter.initialGuess )
    {
        INFO( "%g ", coord );
    }
    INFO( "Particle types to be used as recenter anchors:" );
    for ( auto& id : para.orbit->recenter.anchorIds )
    {
        INFO( "%d ", id );
    }
}

/**
 * @brief The component part.
 *
 * @param para reference of the runtime parameter
 */
void print_component_part( otf::runtime_para& para )
{
    INFO( "COMPONENT PARAMETERS" );
    for ( auto& comp : para.comps )
    {
        INFO( "Get an component: %s", comp.first.c_str() );
        INFO( "Analysis period of this component: %d", comp.second->period );
        INFO( "Particle types in this component:" );
        for ( auto& id : comp.second->types )
        {
            INFO( "%d ", id );
        }
        if ( comp.second->recenter.enable )
        {
            INFO( "Recenter of this component is enabled." );
        }
        string method;
        switch ( comp.second->recenter.method )
        {
        case otf::recenter_method::COM:
            method = "center of mass";
            break;
        case otf::recenter_method::MBP:
            method = "most bound particle";
            break;
        default:
            ERROR( "Get into an unexpected branch!" );
        }
        INFO( "Potential method for recenter: %s.", method.c_str() );
        INFO( "Potential enclosed radius for recenter: %g.", comp.second->recenter.radius );
        INFO( "Initial guess of the recenter:" );
        for ( auto& coord : para.orbit->recenter.initialGuess )
        {
            INFO( "%g ", coord );
        }
        if ( comp.second->align.enable )
        {
            INFO( "Alignment of this component is enabled." );
        }
        INFO( "Initial guess of the recenter:" );
        INFO( "Potential enclosed radius for align: %g.", comp.second->align.radius );
        if ( comp.second->image.enable )
        {
            INFO( "Image of this component is enabled." );
        }
        INFO( "Image half length: %g.", comp.second->image.halfLength );
        INFO( "Image bin number: %d.", comp.second->image.binNum );
        if ( comp.second->A2.enable )
        {
            INFO( "A2 of this component is enabled." );
        }
        INFO( "A2 rmin : %g.", comp.second->A2.rmin );
        INFO( "A2 rmax : %g.", comp.second->A2.rmax );

        if ( comp.second->barAngle.enable )
        {
            INFO( "barAngle of this component is enabled." );
        }
        INFO( "barAngle rmin : %g.", comp.second->barAngle.rmin );
        INFO( "barAngle rmax : %g.", comp.second->barAngle.rmax );

        if ( comp.second->buckle.enable )
        {
            INFO( "buckle of this component is enabled." );
        }
        INFO( "buckle rmin : %g.", comp.second->buckle.rmin );
        INFO( "buckle rmax : %g.", comp.second->buckle.rmax );
    }
}

/**
 * @brief Print the basic informations of the runtime parameters.
 *
 * @param para reference of the runtime parameter
 */
void print_para_info( otf::runtime_para& para )
{
    // print global parameters
    print_global_part( para );

    // print orbital log parameters
    print_orbital_part( para );

    // component analysis parameters
    print_component_part( para );
}

namespace otf {

monitor::monitor( const std::string_view& tomlParaFile )
    : mpiRank( -1 ), mpiSize( 0 ), isRootRank( false ), stepCounter( 0 ),
      mpiInitialzedByMonitor( false ), para( runtime_para( tomlParaFile ) ), h5Organizer( nullptr )

{
    // check whether in mpi environment, if not, call MPI_Init
    int initialzed;
    MPI_Initialized( &initialzed );
    if ( initialzed == 0 )
    {
        mpiInitialzedByMonitor = true;
        MPI_Init( nullptr, nullptr );
    }

    // get the rank id
    MPI_Comm_rank( MPI_COMM_WORLD, &mpiRank );
    if ( mpiRank == 0 )
    {
        isRootRank = true;
    }

    // get the mpi size
    MPI_Comm_size( MPI_COMM_WORLD, &mpiSize );

    // read in the parameters
    MPI_INFO( mpiRank, "Read in parameter from %s", tomlParaFile.data() );
    if ( isRootRank )
    {
        print_para_info( para );
        h5Organizer = make_unique< h5_out >( para.outputDir, para.fileName );
    }
}

monitor::~monitor()
{
    if ( mpiInitialzedByMonitor )
    {
        MPI_Finalize();
    }
}

void monitor::one_analysis_api( double time, const unsigned int particleNumber, const int* id,
                                const int* partType, const double* mass, const double* coordinate,
                                const double* velocity )
{
    ( void )time;
    if ( not para.enableOtf )
    {
        return;
    }

    // First: extract the data for orbital log, and the data for each component
    auto orbitalData =
        id_data_process( time, particleNumber, id, partType, mass, coordinate, velocity );
    // if it's the first extraction, create the datasets in the root rank
    if ( isRootRank and stepCounter == 0 )
    {
        // TODO: do something here
    }
    // TODO: selection of each component

    // Second: log the orbits,
    // TODO: log the orbits

    // Third: analyze each component
    // TODO: analyze each component

    // Last: increase the synchronized step counter
    stepCounter++;
}

auto monitor::id_data_process( const double time, const unsigned int particleNumber,
                               const int* particleID, const int* particleType, const double* mass,
                               const double* coordinate,
                               const double* velocity ) const -> vector< orbitPoint >
{
    vector< orbitPoint >             points;
    static const otf::orbit_selector orbitSelector( para );
    auto getData = orbitSelector.select( particleNumber, particleID, particleType, mass, coordinate,
                                         velocity );

    // TODO: MPI collection
    const int                 localNum = getData->count;  // the number of ids in local mpi rank
    const unique_ptr< int[] > numInEachRank( new int[ mpiSize ]() );  // number in each rank
    // collective communication to gather the number of particles in each rank
    MPI_Allgather( &localNum, 1, MPI_INT, numInEachRank.get(), 1, MPI_INT, MPI_COMM_WORLD );

    // get the global total number
    int totalNum = 0;
    for ( int i = 0; i < mpiSize; ++i )
    {
        totalNum += numInEachRank[ i ];
    }

    // the offset of the local mpi rank, used for variable length mpi gathering
    int localOffset = 0;
    for ( int i = 0; i < mpiSize; ++i )
    {
        localOffset += numInEachRank[ i ];
    }

    // the array of offset values
    const unique_ptr< int[] > offsets( new int[ mpiSize ]() );
    MPI_Allgather( &localOffset, 1, MPI_INT, offsets.get(), 1, MPI_INT, MPI_COMM_WORLD );

    // number of data points and offsets of 3D array
    const unique_ptr< int[] > numInEachRank3D( new int[ mpiSize ]() );  // number in each rank
    const unique_ptr< int[] > offsets3D( new int[ mpiSize ]() );
    for ( auto i = 0; i < mpiSize; ++i )
    {
        numInEachRank3D[ i ] = numInEachRank[ i ] * 3;
        offsets3D[ i ]       = offsets[ i ] * 3;
    }

    // global array of datas
    const unique_ptr< int[] >    gIDs( new int[ totalNum ]() );
    const unique_ptr< double[] > gCoordinate( new double[ totalNum * 3 ]() );
    const unique_ptr< double[] > gVelocity( new double[ totalNum * 3 ]() );
    // gather ids
    MPI_Allgatherv( getData->id.data(), localNum, MPI_INT, gIDs.get(), numInEachRank.get(),
                    offsets.get(), MPI_INT, MPI_COMM_WORLD );
    // gather coordinates
    MPI_Allgatherv( getData->coordinate.data(), localNum * 3, MPI_DOUBLE, gCoordinate.get(),
                    numInEachRank3D.get(), offsets3D.get(), MPI_DOUBLE, MPI_COMM_WORLD );
    // gather velocities
    MPI_Allgatherv( getData->velocity.data(), localNum * 3, MPI_DOUBLE, gVelocity.get(),
                    numInEachRank3D.get(), offsets3D.get(), MPI_DOUBLE, MPI_COMM_WORLD );

    ( void )time;
    return points;
}

}  // namespace otf
