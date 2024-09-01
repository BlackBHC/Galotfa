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
using namespace std;

void print_para_info( otf::runtime_para& para )
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

    INFO( "Method of id determination: " );
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

namespace otf {

monitor::monitor( const std::string_view& tomlParaFile )
    : mpiRank( -1 ), isRootRank( false ), stepCounter( 0 ), mpiInitialzedByMonitor( false ),
      para( runtime_para( tomlParaFile ) ), h5Organizer( nullptr )

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

void monitor::one_analysis_api( const unsigned int particleNumber, const int* id,
                                const int* partType, const double* mass, const double* coordinate,
                                const double* velocity )
{
    if ( not para.enableOtf )
    {
        return;
    }

    // First: extract the data for orbital log, and the data for each component
    auto orbitalData = id_data_process( particleNumber, id, partType, mass, coordinate, velocity );
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

auto monitor::id_data_process( unsigned int particleNumber, const int* particleID,
                               const int* particleType, const double* mass,
                               const double* coordinate,
                               const double* velocity ) const -> std::unique_ptr< dataContainer >
{
    static otf::orbit_selector orbitSelector( para );
    auto getData = orbitSelector.select( particleNumber, particleID, particleType, mass, coordinate,
                                         velocity );
    // TODO: MPI collection
    return getData;
}

}  // namespace otf
