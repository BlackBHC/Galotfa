/**
 * @file monitor.cpp
 * @brief The organizer of other components to work together.
 */

#include <utility>
#ifdef DEBUG
#include "../include/myprompt.hpp"
#endif
#include "../include/barinfo.hpp"
#include "../include/eigen.hpp"
#include "../include/h5out.hpp"
#include "../include/monitor.hpp"
#include "../include/para.hpp"
#include "../include/recenter.hpp"
#include "../include/selector.hpp"
#include <H5Tpublic.h>
#include <algorithm>
#include <cmath>
#include <memory>
#include <mpi.h>
#include <string>
#include <string_view>
#include <vector>
using namespace std;

// inline function to calculate the determinant of an matrix
inline auto determinant( const double* matrix ) -> double
{
    double det = 0;
    det += matrix[ 0 ] * matrix[ 4 ] * matrix[ 8 ] + matrix[ 1 ] * matrix[ 5 ] * matrix[ 6 ]
           + matrix[ 2 ] * matrix[ 3 ] * matrix[ 7 ];

    det -= matrix[ 2 ] * matrix[ 4 ] * matrix[ 6 ] + matrix[ 1 ] * matrix[ 3 ] * matrix[ 8 ]
           + matrix[ 0 ] * matrix[ 5 ] * matrix[ 7 ];

    return det;
}

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
        if ( comp.second->sBar.enable )
        {
            INFO( "A2 of this component is enabled." );
        }
        INFO( "A2 rmin : %g.", comp.second->sBar.rmin );
        INFO( "A2 rmax : %g.", comp.second->sBar.rmax );

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

monitor::monitor( const string_view& tomlParaFile )
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
 * @brief The main analysis API, which should be called in the main loop of the simulation.
 *
 * @param time time of the simulation
 * @param particleNumber number of particles in the local mpi rank
 * @param particleID ids of particles
 * @param particleType PartTypes of particles
 * @param mass masses of particles
 * @param coordinate coordinates of particles
 * @param velocity velocities of particles
 */
void monitor::main_analysis_api( const double time, const unsigned particleNumber, const int* ids,
                                 const int* partTypes, const double* masses,
                                 const double* potentials, const double* coordinates,
                                 const double* velocities )
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
        component_analysis( time, particleNumber, partTypes, masses, potentials, coordinates,
                            velocities, comp.second );
    }

    // Last: increase the synchronized step counter
    stepCounter++;
}

/**
 * @brief The API of orbital log.
 *
 * @param time time of the simulation
 * @param particleNumber number of particles in the local mpi rank
 * @param particleID ids of particles
 * @param particleType PartTypes of particles
 * @param mass masses of particles
 * @param coordinate coordinates of particles
 * @param velocity velocities of particles
 */
void monitor::orbital_log( const double time, const unsigned particleNumber, const int* ids,
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
 * @brief The API of data extraction for component analysis.
 *
 * @param time time of the simulation
 * @param particleNumber number of particles in the local mpi rank
 * @param particleType PartTypes of particles
 * @param mass masses of particles
 * @param coordinate coordinates of particles
 * @param velocity velocities of particles
 * @param comp otf::component object, a structure of parameters for a component
 * @return an object of otf::monitor::compDataContainer, which is the container of the extracted
 * data
 */
auto monitor::component_data_extract(
    unsigned particleNumber, const int* partType, const double* masses, const double* potentials,
    const double* coordinates, const double* velocities,
    unique_ptr< otf::component >& comp ) -> monitor::compDataContainer
{
    unsigned         count = 0;
    vector< double > extractedMasses;
    vector< double > extractedPotentials;
    vector< double > extractedCoordinates;
    vector< double > extractedVelocities;
    extractedMasses.resize( particleNumber );
    extractedPotentials.resize( particleNumber );
    extractedCoordinates.resize( particleNumber * 3 );
    extractedVelocities.resize( particleNumber * 3 );
    for ( unsigned i = 0; i < particleNumber; ++i )
    {
        // check whether it's a particle with the specified id
        auto find_res = find( comp->types.begin(), comp->types.end(), partType[ i ] );
        // if not found, go to the next loop
        if ( find_res == comp->types.end() )
        {
            continue;
        }

        // get the extracted data
        extractedMasses[ count ]              = masses[ i ];
        extractedPotentials[ count ]          = potentials[ i ];
        extractedCoordinates[ count * 3 + 0 ] = coordinates[ i * 3 + 0 ];
        extractedCoordinates[ count * 3 + 1 ] = coordinates[ i * 3 + 1 ];
        extractedCoordinates[ count * 3 + 2 ] = coordinates[ i * 3 + 2 ];
        extractedVelocities[ count * 3 + 0 ]  = velocities[ i * 3 + 0 ];
        extractedVelocities[ count * 3 + 1 ]  = velocities[ i * 3 + 1 ];
        extractedVelocities[ count * 3 + 2 ]  = velocities[ i * 3 + 2 ];

        // increase the particle count
        ++count;
    }

    // remove the tail garbage values
    extractedMasses.resize( count );
    extractedPotentials.resize( count );
    extractedCoordinates.resize( count * 3 );
    extractedVelocities.resize( count * 3 );

    compDataContainer compData;
    compData.partNum = count;
    unique_ptr< double[] > massPtr( new double[ count ]() );
    unique_ptr< double[] > potPtr( new double[ count ]() );
    unique_ptr< double[] > posPtr( new double[ count * 3 ]() );
    unique_ptr< double[] > velPtr( new double[ count * 3 ]() );

    // get the extracted data
    for ( unsigned i = 0; i < count; ++i )
    {
        massPtr[ i ]        = extractedMasses[ i ];
        potPtr[ i ]         = extractedPotentials[ i ];
        posPtr[ i * 3 + 0 ] = extractedCoordinates[ i * 3 + 0 ];
        posPtr[ i * 3 + 1 ] = extractedCoordinates[ i * 3 + 1 ];
        posPtr[ i * 3 + 2 ] = extractedCoordinates[ i * 3 + 2 ];
        velPtr[ i * 3 + 0 ] = extractedVelocities[ i * 3 + 0 ];
        velPtr[ i * 3 + 1 ] = extractedVelocities[ i * 3 + 1 ];
        velPtr[ i * 3 + 2 ] = extractedVelocities[ i * 3 + 2 ];
    }

    // move the data to the container
    compData.masses      = std::move( massPtr );
    compData.potentials  = std::move( potPtr );
    compData.coordinates = std::move( posPtr );
    compData.velocities  = std::move( velPtr );

    // TODO: test whether this retern can work correctly
    return compData;
}

/**
 * @brief The API of analysis part for a single component.
 *
 * @param dataContainer reference to the extracted data, in the form of compDataContainer
 * @param comp wrapper of parameters for analysis of a single component
 * @return the data container of the analysis results
 */
auto monitor::component_data_analyze( monitor::compDataContainer&        dataContainer,
                                      std::unique_ptr< otf::component >& comp )
    -> monitor::compResContainer
{
    compResContainer compRes;

    // NOTE: recenter the system if necessary
    if ( comp->recenter.enable )  // if not enable, do nothing
    {
        recenter_coordinate( dataContainer, comp );
    }

    // TODO: align the system if necessary
    if ( comp->align.enable )
    {
        align_coordinate( dataContainer, comp );
    }

    // NOTE: calculate the bar info if necessary: Sbar, Sbuckle, bar angle and
    // TODO: bar length
    if ( comp->sBar.enable or comp->barAngle.enable or comp->sBuckle.enable )
    {
        bar_info( dataContainer, comp, compRes );
    }

    // TODO: calculate the image if necessary
    image();

    // TODO: test whether the data in unique_ptr can be return in this form
    return compRes;
}

/**
 * @brief The API to recenter the coordinates in a data container object.
 *
 * @param dataContainer reference to the data container to be recenterred
 * @param comp wrapper of parameters for analysis of a single component
 */
void monitor::recenter_coordinate( monitor::compDataContainer&        dataContainer,
                                   std::unique_ptr< otf::component >& comp )
{
    // get the system center
    auto center = recenter::get_center( comp->recenter.method, dataContainer.partNum,
                                        dataContainer.masses.get(), dataContainer.potentials.get(),
                                        dataContainer.coordinates.get(), comp->recenter.radius );
    // substract the system center
    for ( unsigned i = 0; i < dataContainer.partNum; ++i )
    {
        for ( unsigned j = 0; j < 3; ++j )
        {
            dataContainer.coordinates[ i * 3 + j ] = center[ j ];
        }
    };
}

/**
 * @brief The API to align the coordinates in a data container object.
 *
 * @param dataContainer reference to the data container
 * @param comp wrapper of parameters for analysis of a single component
 */
void monitor::align_coordinate( monitor::compDataContainer&        dataContainer,
                                std::unique_ptr< otf::component >& comp )
{
    // get the intertia tensor
    double inertiaTensor[ 9 ] = { 0, 0, 0, 0, 0, 0, 0, 0, 0 };
    for ( unsigned i = 0; i < dataContainer.partNum; ++i )
    {
        // get the spherical radius of the particle
        static double radius;
        radius = sqrt(
            dataContainer.coordinates[ i * 3 + 0 ] * dataContainer.coordinates[ i * 3 + 0 ]
            + dataContainer.coordinates[ i * 3 + 1 ] * dataContainer.coordinates[ i * 3 + 1 ]
            + dataContainer.coordinates[ i * 3 + 2 ] * dataContainer.coordinates[ i * 3 + 2 ] );

        // check whether the particle locates in the enclosed radius
        if ( comp->align.radius < radius )
        {
            continue;
        }

        // diagonal terms
        inertiaTensor[ 0 * 3 + 0 ] +=
            dataContainer.masses[ i ]
            * ( dataContainer.coordinates[ i * 3 + 1 ] * dataContainer.coordinates[ i * 3 + 1 ]
                + dataContainer.coordinates[ i * 3 + 2 ] * dataContainer.coordinates[ i * 3 + 2 ] );
        inertiaTensor[ 1 * 3 + 1 ] +=
            dataContainer.masses[ i ]
            * ( dataContainer.coordinates[ i * 3 + 0 ] * dataContainer.coordinates[ i * 3 + 0 ]
                + dataContainer.coordinates[ i * 3 + 2 ] * dataContainer.coordinates[ i * 3 + 2 ] );
        inertiaTensor[ 2 * 3 + 2 ] +=
            dataContainer.masses[ i ]
            * ( dataContainer.coordinates[ i * 3 + 0 ] * dataContainer.coordinates[ i * 3 + 0 ]
                + dataContainer.coordinates[ i * 3 + 1 ] * dataContainer.coordinates[ i * 3 + 1 ] );
        // non-diagonal terms
        inertiaTensor[ 0 * 3 + 1 ] += -dataContainer.masses[ i ]
                                      * dataContainer.coordinates[ i * 3 + 0 ]
                                      * dataContainer.coordinates[ i * 3 + 1 ];
        inertiaTensor[ 0 * 3 + 2 ] += -dataContainer.masses[ i ]
                                      * dataContainer.coordinates[ i * 3 + 0 ]
                                      * dataContainer.coordinates[ i * 3 + 2 ];
        inertiaTensor[ 1 * 3 + 0 ] += -dataContainer.masses[ i ]
                                      * dataContainer.coordinates[ i * 3 + 1 ]
                                      * dataContainer.coordinates[ i * 3 + 0 ];
        inertiaTensor[ 1 * 3 + 2 ] += -dataContainer.masses[ i ]
                                      * dataContainer.coordinates[ i * 3 + 1 ]
                                      * dataContainer.coordinates[ i * 3 + 2 ];
        inertiaTensor[ 2 * 3 + 0 ] += -dataContainer.masses[ i ]
                                      * dataContainer.coordinates[ i * 3 + 2 ]
                                      * dataContainer.coordinates[ i * 3 + 0 ];
        inertiaTensor[ 2 * 3 + 1 ] += -dataContainer.masses[ i ]
                                      * dataContainer.coordinates[ i * 3 + 2 ]
                                      * dataContainer.coordinates[ i * 3 + 1 ];
    }
    // reduce the inertiaTensor from all mpi ranks
    MPI_Allreduce( MPI_IN_PLACE, inertiaTensor, 9, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD );

    // get the eigenvalues and eigenvectors
    double eigenValues[ 3 ];
    double eigenVectors[ 9 ];
    eigen::eigens_sym_33( inertiaTensor, eigenValues, eigenVectors );

    // make sure it's a rotation matrix
    if ( determinant( eigenVectors ) < 0 )
    {
        eigenVectors[ 2 ] *= -1;
        eigenVectors[ 5 ] *= -1;
        eigenVectors[ 8 ] *= -1;
    }
    // NOTE: rotation matrix is Transpose(EigenMatrix) x Identity

    // rotate the coordinates and velocities
    static double x = 0;
    static double y = 0;
    static double z = 0;
    for ( unsigned i = 0; i < dataContainer.partNum; ++i )
    {
        // coordinates
        x = dataContainer.coordinates[ i * 3 + 0 ];
        y = dataContainer.coordinates[ i * 3 + 1 ];
        z = dataContainer.coordinates[ i * 3 + 2 ];
        dataContainer.coordinates[ i * 3 + 0 ] =
            eigenVectors[ 0 ] * x + eigenVectors[ 3 ] * y + eigenVectors[ 6 ] * z;
        dataContainer.coordinates[ i * 3 + 1 ] =
            eigenVectors[ 1 ] * x + eigenVectors[ 4 ] * y + eigenVectors[ 7 ] * z;
        dataContainer.coordinates[ i * 3 + 2 ] =
            eigenVectors[ 2 ] * x + eigenVectors[ 5 ] * y + eigenVectors[ 8 ] * z;

        // velocities
        x = dataContainer.velocities[ i * 3 + 0 ];
        y = dataContainer.velocities[ i * 3 + 1 ];
        z = dataContainer.velocities[ i * 3 + 2 ];
        dataContainer.velocities[ i * 3 + 0 ] =
            eigenVectors[ 0 ] * x + eigenVectors[ 3 ] * y + eigenVectors[ 6 ] * z;
        dataContainer.velocities[ i * 3 + 1 ] =
            eigenVectors[ 1 ] * x + eigenVectors[ 4 ] * y + eigenVectors[ 7 ] * z;
        dataContainer.velocities[ i * 3 + 2 ] =
            eigenVectors[ 2 ] * x + eigenVectors[ 5 ] * y + eigenVectors[ 8 ] * z;
    }
    // TODO: test the rotation part
}

void monitor::bar_info( monitor::compDataContainer&        dataContainer,
                        std::unique_ptr< otf::component >& comp, compResContainer& res )
{
    // calculate the azimuthal angles
    unique_ptr< double[] > phis( new double[ dataContainer.partNum ] );
    for ( unsigned i = 0; i < dataContainer.partNum; ++i )
    {
        phis[ i ] =
            atan2( dataContainer.coordinates[ 3 * i + 1 ], dataContainer.coordinates[ 3 * i + 0 ] );
    }

    // BUG: doesn't consider the enclosed radius!!!
    bar_info::A2info info;
    if ( comp->barAngle.enable or comp->sBar.enable )
    {
        info = bar_info::A2( dataContainer.partNum, dataContainer.masses.get(), phis.get() );

        // NOTE: bar angle
        if ( comp->barAngle.enable )
        {
            res.barAngle = info.phase;
        }

        // NOTE: bar strength
        if ( comp->sBar.enable )
        {
            double A0 = bar_info::A0( dataContainer.partNum, dataContainer.masses.get() );
            res.sBar  = info.amplitude / A0;
        }
    };

    if ( comp->sBuckle.enable )
    {
        ;
    }
}

void monitor::image()
{
    ;
}

/**
 * @brief The API of analysis part for a single component
 *
 * @param time time of the simulation
 * @param particleNumber number of particles in the local mpi rank
 * @param particleType PartTypes of particles
 * @param mass masses of particles
 * @param coordinate coordinates of particles
 * @param velocity velocities of particles
 * @param comp otf::component object, a structure of parameters for a component
 */
void monitor::component_analysis( double time, unsigned particleNumber, const int* partTypes,
                                  const double* masses, const double* potentials,
                                  const double* coordinates, const double* velocities,
                                  unique_ptr< otf::component >& comp ) const
{
    // NOTE: collect the component data
    auto compDataContainer = component_data_extract( particleNumber, partTypes, masses, potentials,
                                                     coordinates, velocities, comp );

    // NOTE: get the analysis result
    auto compResContainer = component_data_analyze( compDataContainer, comp );

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
        if ( comp->sBar.enable )
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
        if ( comp->sBar.enable )
        {
            h5Organizer->flush_single_block( comp->compName, "A2", &compResContainer.sBar );
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
auto monitor::id_data_process( const double time, const unsigned particleNumber,
                               const int* particleIDs, const int* particleTypes,
                               const double* masses, const double* coordinates,
                               const double* velocities ) const -> vector< orbitPoint >
{
    vector< orbitPoint >             points;
    static const otf::orbit_selector orbitSelector( para );
    auto getData = orbitSelector.select( particleNumber, particleIDs, particleTypes, masses,
                                         coordinates, velocities );

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
    sort( points.begin(), points.end(), cmp );
    return points;
}

}  // namespace otf
