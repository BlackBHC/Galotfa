#include "../include/para.hpp"
#include "../include/myprompt.hpp"
#include <string>
#include <string_view>
using namespace std;


runtime_para::runtime_para( const std::string_view& tomlParaFile )
    : paraTable( toml::parse_file( tomlParaFile ) )
{
    ;
}

void runtime_para::read_one_by_one()
{
    // auto config = toml::parse_file( "/Users/chenbinhui/codes/galotfa/examples/galotfa.toml" );
    //
    // // get key-value pairs
    // // int
    // int maxIter = config[ "global" ][ "maxiter" ].value_or( 0 );
    // // double
    // double epsilon = config[ "global" ][ "epsilon" ].value_or( 0.0 );
    // // bool
    // bool enableRecenter = config[ "component1" ][ "recenter" ][ "enable" ].value_or( true );
    // // string
    // string_view outDir   = config[ "global" ][ "outdir" ].value_or( "./otfLogs"sv );
    // string_view filename = config[ "global" ][ "filename" ].value_or( "galotfa.hdf5"sv );
    //
    // cout << "Output at [" << outDir << "]: [" << filename << "]" << endl;
    // cout << "Max iteration:" << maxIter << endl;
    // cout << "Equal threshold of floating-point number: " << epsilon << endl;
    // if ( enableRecenter )
    //     cout << "The recenter is enabled." << endl;
    //
    // // array
    // auto   iguess = config[ "component1" ][ "recenter" ][ "iguess" ];
    // double x      = iguess[ 0 ].value< double >().value();
    // double y      = iguess[ 1 ].value< double >().value();
    // double z      = iguess[ 2 ].value< double >().value();
    // std::cout << "iguess: " << x << ", " << y << ", " << z << "\n";
    ;
}

component::component( string_view& compName, toml::table& para ) : compName( compName )
{
    auto compNode = para[ compName ];
    // period
    period = compNode[ "period" ].value< double >().value();

    // recenter parameters
    recenterPara.enable = compNode[ "recenter" ][ "enable" ].value< bool >().value();
    recenterPara.radius = compNode[ "recenter" ][ "radius" ].value< bool >().value();
    auto str            = compNode[ "recenter" ][ "method" ].value< string_view >().value();
    if ( str == "com" )
    {
        recenterPara.method = recenter_method::COM;
    }
    else if ( str == "mbp" )
    {
        recenterPara.method = recenter_method::MBP;
    }
    else
    {
        ERROR( "Get an unknown value for [recenter method] of component [%s]", compName.data() );
        ERROR( "Must be 'com' (for center of mass) or 'mbp' (for most bound particle)." );
        exit( -1 );
    }
    auto iguess = compNode[ "recenter" ][ "iguess" ];
    for ( auto i = 0; i < 3; ++i )
        recenterPara.initialGuess[ i ] = iguess[ i ].value< double >().value();


    // frame
    str = compNode[ "frame" ].value< string_view >().value();
    if ( str == "cyl" )
    {
        frame = coordinate_frame::CYLINDRICAL;
    }
    else if ( str == "sph" )
    {
        frame = coordinate_frame::SPHERICAL;
    }
    else if ( str == "car" )
    {
        frame = coordinate_frame::CARTESIAN;
    }
    else
    {
        ERROR( "Get an unknown value for [coordinate frame] of component [%s]", compName.data() );
        ERROR( "Must be one of 'cyl' (for cylindrical), 'car' (for Cartesian), or 'sph' (for "
               "spherical)." );
        exit( -1 );
    }

    // align
    alignPara.enable = compNode[ "align" ][ "enable" ].value< bool >().value();
    alignPara.radius = compNode[ "align" ][ "radius" ].value< double >().value();

    // image
    imagePara.enable     = compNode[ "image" ][ "enable" ].value< bool >().value();
    imagePara.halfLength = compNode[ "image" ][ "halflength" ].value< double >().value();
    imagePara.binNum     = compNode[ "image" ][ "binnum" ].value< unsigned int >().value();

    // bar info parameters
    A2Para.enable       = compNode[ "A2" ][ "enable" ].value< bool >().value();
    A2Para.rmin         = compNode[ "A2" ][ "rmin" ].value< double >().value();
    A2Para.rmax         = compNode[ "A2" ][ "rmax" ].value< double >().value();
    barAnglePara.enable = compNode[ "barangle" ][ "enable" ].value< bool >().value();
    barAnglePara.rmin   = compNode[ "barangle" ][ "rmin" ].value< double >().value();
    barAnglePara.rmax   = compNode[ "barangle" ][ "rmax" ].value< double >().value();
    bucklePara.enable   = compNode[ "buckle" ][ "enable" ].value< bool >().value();
    bucklePara.rmin     = compNode[ "buckle" ][ "rmin" ].value< double >().value();
    bucklePara.rmax     = compNode[ "buckle" ][ "rmax" ].value< double >().value();
}
