/**
 * @file
 * @brief Test the parameter parser.
 */

#include "../include/myprompt.hpp"
#include "../include/para.hpp"
#include <hdf5.h>
using namespace std;
using namespace otf;

int main()
{
    runtime_para para( "../examples/galotfa.toml" );

    INFO( "GLOBAL PARAMETERS:" );
    if ( para.enableOtf )
    {
        INFO( "On the fly analysis is enabled." );
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
    for ( auto& id : para.orbit->logTypes )
    {
        INFO( "%d ", id );
    }

    INFO( "Method of id determination: " );
    switch ( para.orbit->method )
    {
    case orbit::log_method::RANDOM:
        INFO( "Random selection." );
        break;
    case orbit::log_method::TXTFILE:
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


    return 0;
}
