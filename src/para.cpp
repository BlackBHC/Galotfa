#include "../include/para.hpp"
#include "../include/myprompt.hpp"
#include <string_view>
using namespace std;


runtime_para::runtime_para( const std::string_view& tomlParaFile )
    : paraTable( toml::parse_file( tomlParaFile ) )
{
    ;
}

void runtime_para::read_one_by_one()
{
    // auto config = toml::parse_file( "./galotfa.toml" );
    // // get key-value pairs
    // // int
    // int maxIter = config[ "global" ][ "maxiter" ].value_or( 0 );
    // // double
    // double epsilon = config[ "global" ][ "epsilon" ].value_or( 0.0 );
    // // bool
    // // bool enableRecenter = config[ "component1" ][ "recenter.enable" ].value_or( true );
    // bool enableRecenter = *config[ "component1" ][ "recenter" ][ "enable" ].value< bool >();
    // // string
    // string_view outDir   = config[ "global.outdir" ].value_or( "./otfLogs"sv );
    // string_view filename = config[ "global" ][ "filename" ].value_or( "galotfa.hdf5"sv );
    //
    // cout << "Output at [" << outDir << "]: [" << filename << "]" << endl;
    // cout << "Max iteration:" << maxIter << endl;
    // cout << "Equal threshold of floating-point number: " << epsilon << endl;
    // if ( enableRecenter )
    //     cout << "The recenter is enabled." << endl;
    // else
    //     cout << "The recenter is disabled." << endl;
    //
    // if ( filename == "galotfa.hdf5" )
    //     cout << "Using default filename" << endl;
    //
    // // array
    // auto   iguess = config[ "component1" ][ "recenter" ][ "iguess" ];
    // double x      = *iguess[ 0 ].value< unsigned int >();
    // double y      = *iguess[ 1 ].value< double >();
    // double z      = *iguess[ 2 ].value< double >();
    // std::cout << "iguess: " << x << ", " << y << ", " << z << "\n";
    //
    // auto numbers = config[ "orbit" ][ "logtypes" ];
    // if ( toml::array* arr = numbers.as_array() )
    // {
    //     // visitation with for_each() helps deal with heterogeneous data
    //     arr->for_each( []( auto el ) {
    //         if constexpr ( toml::is_number< decltype( el ) > )
    //         {
    //             double tmp = *el;
    //             cout << tmp << endl;
    //         }
    //     } );
    // }
}

component::component( string_view& compName, toml::table& para ) : compName( compName )
{
    auto compNode = para[ compName ];
    // period
    period = *compNode[ "period" ].value< double >();

    // recenter parameters
    recenter.enable = *compNode[ "recenter" ][ "enable" ].value< bool >();
    if ( recenter.enable )
    {
        recenter.radius = *compNode[ "recenter" ][ "radius" ].value< bool >();
        auto str        = *compNode[ "recenter" ][ "method" ].value< string_view >();
        if ( str == "com" )
        {
            recenter.method = recenter_method::COM;
        }
        else if ( str == "mbp" )
        {
            recenter.method = recenter_method::MBP;
        }
        else
        {
            ERROR( "Get an unknown value for [recenter method] of component [%s]: [%s]",
                   compName.data(), str.data() );
            ERROR( "Must be 'com' (for center of mass) or 'mbp' (for most bound particle)." );
            exit( -1 );
        }
        auto iguess = compNode[ "recenter" ][ "iguess" ];
        for ( auto i = 0; i < 3; ++i )
            recenter.initialGuess[ i ] = *iguess[ i ].value< double >();
    }

    // NOTE: the frame parameter is unused at present
    //
    // // frame
    // auto str = *compNode[ "frame" ].value< string_view >();
    // if ( str == "cyl" )
    // {
    //     frame = coordinate_frame::CYLINDRICAL;
    // }
    // else if ( str == "sph" )
    // {
    //     frame = coordinate_frame::SPHERICAL;
    // }
    // else if ( str == "car" )
    // {
    //     frame = coordinate_frame::CARTESIAN;
    // }
    // else
    // {
    //     ERROR( "Get an unknown value for [coordinate frame] of component [%s]: [%s]",
    //            compName.data(), str.data() );
    //     ERROR( "Must be one of 'cyl' (for cylindrical), 'car' (for Cartesian), or 'sph' (for "
    //            "spherical)." );
    //     exit( -1 );
    // }

    // align
    align.enable = *compNode[ "align" ][ "enable" ].value< bool >();
    if ( align.enable )
    {
        align.radius = *compNode[ "align" ][ "radius" ].value< double >();
    }

    // image
    image.enable = *compNode[ "image" ][ "enable" ].value< bool >();
    if ( image.enable )
    {
        image.halfLength = *compNode[ "image" ][ "halflength" ].value< double >();
        image.binNum     = *compNode[ "image" ][ "binnum" ].value< unsigned int >();
    }

    // bar info parameters
    // A2
    A2.enable = *compNode[ "A2" ][ "enable" ].value< bool >();
    if ( A2.enable )
    {
        A2.rmin = *compNode[ "A2" ][ "rmin" ].value< double >();
        A2.rmax = *compNode[ "A2" ][ "rmax" ].value< double >();
    }
    // bar angle
    barAngle.enable = *compNode[ "barangle" ][ "enable" ].value< bool >();
    if ( barAngle.enable )
    {
        barAngle.rmin = *compNode[ "barangle" ][ "rmin" ].value< double >();
        barAngle.rmax = *compNode[ "barangle" ][ "rmax" ].value< double >();
    }
    // buckling strength
    buckle.enable = *compNode[ "buckle" ][ "enable" ].value< bool >();
    if ( buckle.enable )
    {
        buckle.rmin = *compNode[ "buckle" ][ "rmin" ].value< double >();
        buckle.rmax = *compNode[ "buckle" ][ "rmax" ].value< double >();
    }
}

orbit::orbit( toml::table& para )
{
    auto node = para[ "orbit" ];
    // whether enable orbital logs
    enable = *node[ "enable" ].value< bool >();

    if ( not enable )
    {
        return;
    }

    // period
    period = *node[ "period" ].value< double >();

    // particle types to be logged
    auto numbers = node[ "logtypes" ];
    if ( toml::array* arr = numbers.as_array() )
    {
        // visitation with for_each() helps deal with heterogeneous data
        arr->for_each( [ this ]( auto&& el ) {
            if constexpr ( toml::is_number< decltype( el ) > )
            {
                logTypes.push_back( *el );
            }
        } );
    }

    auto str = *node[ "method" ].value< string_view >();
    if ( str == "txtfile" )
    {
        method = log_method::TXTFILE;
    }
    else if ( str == "random" )
    {
        method = log_method::RANDOM;
    }
    else
    {
        ERROR( "Get an unknown value for [orbital log method]: [%s]", str.data() );
        ERROR( "Must be one 'txtfile' (for a text file of id list) or 'random' (for random "
               "selection)." );
        exit( -1 );
    }

    if ( method == log_method::RANDOM )
    {
        // random selection
        random.enable   = *node[ "random" ][ "enable" ].value< bool >();
        random.fraction = *node[ "random" ][ "frac" ].value< double >();
        // recenter.partTypes = node[ "random" ][ "radius" ].value< bool >();
    }
    else
    {
        idfile = *node[ "idfile" ].value< string_view >();
    }

    // recenter parameters
    recenter.enable = *node[ "recenter" ][ "enable" ].value< bool >();
    if ( recenter.enable )
    {
        // particle types to be used as anchors of recenter
        auto numbers = node[ "recenter" ][ "anchorids" ];
        if ( toml::array* arr = numbers.as_array() )
        {
            // visitation with for_each() helps deal with heterogeneous data
            arr->for_each( [ this ]( auto&& el ) {
                if constexpr ( toml::is_number< decltype( el ) > )
                {
                    recenter.anchorIds.push_back( *el );
                }
            } );
        }
        recenter.radius = *node[ "recenter" ][ "radius" ].value< bool >();
        str             = *node[ "recenter" ][ "method" ].value< string_view >();
        if ( str == "com" )
        {
            recenter.method = recenter_method::COM;
        }
        else if ( str == "mbp" )
        {
            recenter.method = recenter_method::MBP;
        }
        else
        {
            ERROR( "Get an unknown value for [recenter method] in orbital logs: [%s]", str.data() );
            ERROR( "Must be 'com' (for center of mass) or 'mbp' (for most bound particle)." );
            exit( -1 );
        }
        auto iguess = node[ "recenter" ][ "iguess" ];
        for ( auto i = 0; i < 3; ++i )
            recenter.initialGuess[ i ] = *iguess[ i ].value< double >();
    }
}
