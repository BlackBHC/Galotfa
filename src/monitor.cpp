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
#include <H5Tpublic.h>
#include <algorithm>
#include <memory>
#include <mpi.h>
#include <string>
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
    INFO( "(%g, %g, %g)", para.orbit->recenter.initialGuess[ 0 ],
          para.orbit->recenter.initialGuess[ 1 ], para.orbit->recenter.initialGuess[ 2 ] );
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
        INFO( "(%g, %g, %g)", para.orbit->recenter.initialGuess[ 0 ],
              para.orbit->recenter.initialGuess[ 1 ], para.orbit->recenter.initialGuess[ 2 ] );
        if ( comp.second->align.enable )
        {
            INFO( "Alignment of this component is enabled." );
        }
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

        if ( comp.second->sBuckle.enable )
        {
            INFO( "buckle of this component is enabled." );
        }
        INFO( "buckle rmin : %g.", comp.second->sBuckle.rmin );
        INFO( "buckle rmax : %g.", comp.second->sBuckle.rmax );
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
    if ( not para.enableOtf )  // if the on-the-fly analysis is not enabled
    {
        return;
    }

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

/**
 * @brief The main analysis api, which should be called in the main loop of the simulation.
 *
 * @param time time of the simulation
 * @param particleNumber number of particles in the local mpi rank
 * @param particleID ids of particles
 * @param particleType PartTypes of particles
 * @param mass masses of particles
 * @param coordinate coordinates of particles
 * @param velocity velocities of particles
 */
void monitor::main_analysis_api( const double time, const unsigned int particleNumber,
                                 const int* ids, const int* partTypes, const double* masses,
                                 const double* coordinates, const double* velocities )
{
    if ( not para.enableOtf )
    {
        return;
    }

    // First: orbital logs part
    if ( para.orbit->enable )
    {
        orbital_log( time, particleNumber, ids, partTypes, masses, coordinates, velocities );
    }

    // Second: analyze each component
    // TODO: analyze each component
    for ( auto& comp : para.comps )
    {
        component_analysis( time, particleNumber, partTypes, masses, coordinates, velocities,
                            comp.second );
    }

    // Last: increase the synchronized step counter
    stepCounter++;
}

/**
 * @brief The api of orbital log.
 *
 * @param time time of the simulation
 * @param particleNumber number of particles in the local mpi rank
 * @param particleID ids of particles
 * @param particleType PartTypes of particles
 * @param mass masses of particles
 * @param coordinate coordinates of particles
 * @param velocity velocities of particles
 */
void monitor::orbital_log( const double time, const unsigned int particleNumber, const int* ids,
                           const int* partTypes, const double* masses, const double* coordinates,
                           const double* velocities )
{
    if ( stepCounter % para.orbit->period != 0 )  // only log in the chosen steps
    {
        return;
    }

    // First: extract the data for orbital log, and the data for each component
    auto orbitData =
        id_data_process( time, particleNumber, ids, partTypes, masses, coordinates, velocities );
    // if it's the first extraction, create the datasets in the root rank
    if ( isRootRank and stepCounter == 0 )
    {
        for ( auto& data : orbitData )
        {
            // create the datasets of particles' orbit
            const string datasetName = "Particle-" + to_string( data.particleID );
            h5Organizer->create_dataset_in_group( datasetName, "Orbit", { orbitPointDim },
                                                  H5T_NATIVE_DOUBLE );
            // backup the dataset names
            orbitDatasetNames.push_back( datasetName );
        }
    }

    // Second: log the orbits in the root rank
    if ( isRootRank )
    {
        for ( auto i = 0UL; i < orbitDatasetNames.size(); ++i )
        {
#ifdef DEBUG
            auto returnCode = h5Organizer->flush_single_block( "Orbit", orbitDatasetNames[ i ],
                                                               orbitData[ i ].data );
            if ( returnCode != 0 )
            {
                ERROR( "The dataset [Orbit/%s] written faild!", orbitDatasetNames[ i ].c_str() );
                throw;
            }
#else
            h5Organizer->flush_single_block( "Orbit", orbitDatasetNames[ i ].c_str(),
                                             orbitData[ i ].data );
#endif
        }
    }
}

/**
 * @brief The api of data extraction for component analysis.
 *
 * @param time time of the simulation
 * @param particleNumber number of particles in the local mpi rank
 * @param particleType PartTypes of particles
 * @param mass masses of particles
 * @param coordinate coordinates of particles
 * @param velocity velocities of particles
 * @param comp otf::component object, a structure of parameters for a component
 * @return the vector of otf::monitor::compDataContainer objects, which are containers of analysis
 * result
 */
auto monitor::component_data_extract(
    unsigned int particleNumber, const int* partType, const double* mass, const double* coordinate,
    const double* velocity, unique_ptr< otf::component >& comp ) const -> monitor::compDataContainer
{
    compDataContainer compData;
    ( void )particleNumber;
    ( void )partType;
    ( void )mass;
    ( void )coordinate;
    ( void )velocity;
    ( void )comp;
    return compData;
}

/**
 * @brief The api of analysis part for a single component.
 *
 * @param dataContainer reference to the extracted data, in the form of compDataContainer
 * @return the data container of the analysis results
 */
auto monitor::component_data_analyze( monitor::compDataContainer& dataContainer ) const
    -> monitor::compResContainer
{
    compResContainer compRes;
    ( void )dataContainer;
    return compRes;
}

/**
 * @brief The api of analysis part for a single component
 *
 * @param time time of the simulation
 * @param particleNumber number of particles in the local mpi rank
 * @param particleType PartTypes of particles
 * @param mass masses of particles
 * @param coordinate coordinates of particles
 * @param velocity velocities of particles
 * @param comp otf::component object, a structure of parameters for a component
 */
void monitor::component_analysis( double time, unsigned int particleNumber, const int* partType,
                                  const double* mass, const double* coordinate,
                                  const double* velocity, unique_ptr< otf::component >& comp ) const
{
    // NOTE: collect the component data
    auto compDataContainer =
        component_data_extract( particleNumber, partType, mass, coordinate, velocity, comp );

    // NOTE: get the analysis result
    auto compResContainer = component_data_analyze( compDataContainer );

    // NOTE: create the datasets at the first call
    if ( isRootRank and stepCounter == 0 )
    {
        // create the datasets for times
        if ( comp->recenter.enable )
        {
            h5Organizer->create_dataset_in_group( "Time", comp->compName, { 1 },
                                                  H5T_NATIVE_DOUBLE );
        }

        // create the datasets for center positions
        if ( comp->recenter.enable )
        {
            h5Organizer->create_dataset_in_group( "Center", comp->compName, { 3 },
                                                  H5T_NATIVE_DOUBLE );
        }

        // create the datasets for bar infos
        if ( comp->A2.enable )
        {
            h5Organizer->create_dataset_in_group( "A2", comp->compName, { 1 }, H5T_NATIVE_DOUBLE );
        }
        if ( comp->barAngle.enable )
        {
            h5Organizer->create_dataset_in_group( "BarAngle", comp->compName, { 1 },
                                                  H5T_NATIVE_DOUBLE );
        }
        if ( comp->sBuckle.enable )
        {
            h5Organizer->create_dataset_in_group( "Sbuckle", comp->compName, { 1 },
                                                  H5T_NATIVE_DOUBLE );
        }
        // TODO: bar length

        // create the datasets for image
        if ( comp->image.enable )
        {
            // create the datasets for images in x-y, y-z, x-z planes
            h5Organizer->create_dataset_in_group( "ImageXY", comp->compName,
                                                  { comp->image.binNum, comp->image.binNum },
                                                  H5T_NATIVE_DOUBLE );
            h5Organizer->create_dataset_in_group( "ImageXZ", comp->compName,
                                                  { comp->image.binNum, comp->image.binNum },
                                                  H5T_NATIVE_DOUBLE );
            h5Organizer->create_dataset_in_group( "ImageYZ", comp->compName,
                                                  { comp->image.binNum, comp->image.binNum },
                                                  H5T_NATIVE_DOUBLE );
        }
    }

    // NOTE: flush the data
    if ( isRootRank )
    {
        // time of the current iteration
        h5Organizer->flush_single_block( comp->compName, "Time", &time );

        // center positions
        if ( comp->recenter.enable )
        {
            h5Organizer->flush_single_block( comp->compName, "Center", compResContainer.center );
        }

        // bar infos
        if ( comp->A2.enable )
        {
            h5Organizer->flush_single_block( comp->compName, "A2", &compResContainer.A2 );
        }
        if ( comp->barAngle.enable )
        {
            h5Organizer->flush_single_block( comp->compName, "BarAngle",
                                             &compResContainer.barAngle );
        }
        if ( comp->sBuckle.enable )
        {
            h5Organizer->flush_single_block( comp->compName, "Sbuckle", &compResContainer.sBuckle );
        }
        // TODO: bar length

        // images
        if ( comp->image.enable )
        {
            h5Organizer->flush_single_block( comp->compName, "ImageXY",
                                             compResContainer.imageXZ.get() );
            h5Organizer->flush_single_block( comp->compName, "ImageXZ",
                                             compResContainer.imageXZ.get() );
            h5Organizer->flush_single_block( comp->compName, "ImageYZ",
                                             compResContainer.imageXZ.get() );
        }
    }
}

/**
 * @brief Extract the orbital data points and sort them based on the particle IDs, only the root
 * rank will return the effective data.
 *
 * @param time time of the simulation
 * @param particleNumber number of particles in the local mpi rank
 * @param particleID ids of particles
 * @param particleType PartTypes of particles
 * @param mass masses of particles
 * @param coordinate coordinates of particles
 * @param velocity velocities of particles
 * @return the vector of orbitPoint objects
 */
auto monitor::id_data_process( const double time, const unsigned int particleNumber,
                               const int* particleIDs, const int* particleType, const double* mass,
                               const double* coordinate,
                               const double* velocity ) const -> vector< orbitPoint >
{
    vector< orbitPoint >             points;
    static const otf::orbit_selector orbitSelector( para );
    auto getData = orbitSelector.select( particleNumber, particleIDs, particleType, mass,
                                         coordinate, velocity );

    // NOTE: MPI collection
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
    for ( int i = 0; i < mpiRank; ++i )
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
    MPI_Gatherv( getData->id.data(), localNum, MPI_INT, gIDs.get(), numInEachRank.get(),
                 offsets.get(), MPI_INT, 0, MPI_COMM_WORLD );
    // gather coordinates
    MPI_Gatherv( getData->coordinate.data(), localNum * 3, MPI_DOUBLE, gCoordinate.get(),
                 numInEachRank3D.get(), offsets3D.get(), MPI_DOUBLE, 0, MPI_COMM_WORLD );
    // gather velocities
    MPI_Gatherv( getData->velocity.data(), localNum * 3, MPI_DOUBLE, gVelocity.get(),
                 numInEachRank3D.get(), offsets3D.get(), MPI_DOUBLE, 0, MPI_COMM_WORLD );

    if ( not isRootRank )  // if not root rank, directly return
    {
        return points;
    }

    // NOTE: construct the vector of the orbit data points
    for ( auto i = 0; i < totalNum; ++i )
    {
        points.push_back( orbitPoint(
            { gIDs[ i ],
              { time, gCoordinate[ i * 3 + 0 ], gCoordinate[ i * 3 + 1 ], gCoordinate[ i * 3 + 2 ],
                gVelocity[ i * 3 + 0 ], gVelocity[ i * 3 + 1 ], gVelocity[ i * 3 + 2 ] } } ) );
    }

    // NOTE: sort the orbit points based on their id
    static auto cmp = []( orbitPoint p1, orbitPoint p2 ) { return p1.particleID < p2.particleID; };
    std::sort( points.begin(), points.end(), cmp );
    return points;
}

}  // namespace otf
