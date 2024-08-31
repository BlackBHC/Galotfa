#include "../include/selector.hpp"
#include "../include/myprompt.hpp"
#include "../include/para.hpp"
#include <algorithm>
#include <cstdlib>
#include <fstream>
#include <ios>
#include <iterator>
#include <memory>
#include <mpi.h>
#include <random>
#include <set>
#include <string>
#include <sys/unistd.h>
#include <unistd.h>
#include <utility>
#include <vector>
using namespace std;

namespace otf {

/**
 * @brief Select particle ids with a specified fraction, can be used in any mpi rank.
 *
 * @param raw raw id list
 * @param types pointer to the particle types of simulation data
 * @param sampleTypes vector of the targeted particle types to be sampled
 * @param fraction sampling fraction
 * @return a vector<unsigned int> of the selected id list
 */
auto orbit_selector::id_sample( const vector< unsigned int >& raw, const unsigned int* types,
                                const vector< unsigned int >& sampleTypes,
                                double                        fraction ) -> vector< unsigned int >
{
    if ( fraction <= 0 or fraction > 1 )
    {
        ERROR( "Try to select a illegal fraction: [%lf]", fraction );
    }

    auto                   partNum = raw.size();
    vector< unsigned int > filtered( partNum );
    auto                   counter = 0UL;  // counter of effective particles in this mpi rank
    for ( auto i = 0UL; i < partNum; ++i )
    {
        if ( find( sampleTypes.begin(), sampleTypes.end(), types[ i ] ) != sampleTypes.end() )
        {
            filtered[ counter++ ] = raw[ i ];
        }
    }
    filtered.resize( counter );

    if ( fraction == 1 )
    {
        return filtered;
    }

    auto const selectNum = ( size_t )( counter * fraction );
    int        rank      = -1;
    MPI_Comm_rank( MPI_COMM_WORLD, &rank );
    mpi_print( rank, "Get counter: %lu", selectNum );
    vector< unsigned int > res;
    sample( filtered.begin(), filtered.end(), back_inserter( res ), selectNum,
            mt19937{ random_device{}() } );
    return res;
}

/**
 * @brief Read the ids in the given id list file, can be used in any mpi rank to read the ids into
 * it.
 *
 * @param idFilename filename of the id list file (in ASCII txt format)
 * @return std::vector<int> of the ids.
 */
auto orbit_selector::id_read( const string& idFilename ) -> vector< unsigned int >
{
    // check the availability of the file
    if ( access( idFilename.c_str(), F_OK ) != 0 )
    {
        int rank;
        MPI_Comm_rank( MPI_COMM_WORLD, &rank );
        MPI_ERROR( rank, "Particle ID file not found: [%s]", idFilename.c_str() );
        throw "Particle ID File not found!";
    }

    ifstream               fp;
    vector< unsigned int > ids;
    fp.open( idFilename.c_str(), ios::in );
    while ( !fp.eof() )
    {
        string lineStr;
        getline( fp, lineStr );
        if ( !lineStr.empty() )
        {
            try
            {
                ids.push_back( std::stoi( lineStr ) );
            }
            catch ( ... )
            {
                ERROR( "Get an unexpected value: %s", lineStr.c_str() );
                throw "Get an unexpected value";
            }
        }
    }

    // Remove the possible repeated values
    set< int > tmpSet( ids.begin(), ids.end() );
    ids.assign( tmpSet.begin(), tmpSet.end() );

    return ids;
}

orbit_selector::orbit_selector( const runtime_para& para ) : para( para )
{
    ;
}

auto orbit_selector::select( const unsigned int particleNumber, const unsigned int* particleID,
                             const unsigned int* partType, const double* mass,
                             const double* coordinate,
                             const double* velocity ) const -> std::unique_ptr< dataContainer >
{
    if ( not para.orbit->enable )
    {
        return nullptr;
    }

    // get the target id list based on specified parameters
    static vector< unsigned int > targetIDs;
    // BUG: the selected ids may be move to other mpi by Gadget4!
    // reduce them into a static vector!
    static bool firstCall = true;
    if ( firstCall )
    {
        if ( para.orbit->method == otf::orbit::id_selection_method::RANDOM )
        {
            // random sampling
            vector< unsigned int > rawIds( particleNumber );
            for ( auto i = 0U; i < particleNumber; ++i )
            {
                rawIds[ i ] = particleID[ i ];
            }
            targetIDs =
                id_sample( rawIds, partType, para.orbit->sampleTypes, para.orbit->fraction );
        }
        else
        {
            // read from a txt file
            targetIDs = id_read( para.orbit->idfile );
        }
        firstCall = false;
    }

    // count of found particles
    unsigned int counter = 0;

    // temporary variables restoring extracted data
    vector< double >       tmpMass( particleNumber );
    vector< unsigned int > tmpId( particleNumber );
    vector< double >       tmpPos( particleNumber * 3 );
    vector< double >       tmpVel( particleNumber * 3 );

    for ( auto i = 0U; i < particleNumber; ++i )
    {
        // check whether the id is in the target id list
        if ( find( targetIDs.begin(), targetIDs.end(), particleID[ i ] ) != targetIDs.end() )
        {
            tmpMass[ counter ]        = mass[ i ];
            tmpId[ counter ]          = particleID[ i ];
            tmpPos[ counter * 3 + 0 ] = coordinate[ i * 3 + 0 ];
            tmpPos[ counter * 3 + 1 ] = coordinate[ i * 3 + 1 ];
            tmpPos[ counter * 3 + 2 ] = coordinate[ i * 3 + 2 ];
            tmpVel[ counter * 3 + 0 ] = velocity[ i * 3 + 0 ];
            tmpVel[ counter * 3 + 1 ] = velocity[ i * 3 + 1 ];
            tmpVel[ counter * 3 + 2 ] = velocity[ i * 3 + 2 ];
            // increase the count of particles found
            counter++;
        }
    }

    // remove the rubish value at the end
    tmpMass.resize( counter );
    tmpId.resize( counter );
    tmpPos.resize( counter * 3 );
    tmpVel.resize( counter * 3 );

    // get the data container
    unique_ptr< dataContainer > container = make_unique< dataContainer >();
    container->count                      = counter;
    container->mass                       = std::move( tmpMass );
    container->id                         = std::move( tmpId );
    container->coordinate                 = std::move( tmpPos );
    container->velocity                   = std::move( tmpVel );

    return container;
}

}  // namespace otf
