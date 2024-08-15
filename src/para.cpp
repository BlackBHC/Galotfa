#include "../include/para.hpp"
#include "../include/myprompt.hpp"
#include "../include/recenter.hpp"
#include "../include/toml.hpp"
#include <cassert>
#include <cstdlib>
#include <memory>
#include <string_view>
using namespace std;

namespace otf {

runtime_para::runtime_para( const std::string_view& tomlParaFile )
{
    toml::table paraTable = toml::parse_file( tomlParaFile );

    // check whether enable the on-the-fly analysis
    enableOtf = *paraTable[ "global" ][ "enable" ].value< bool >();
    if ( not enableOtf )
    {
        return;
    }

    // other global parameters
    outputDir = *paraTable[ "global" ][ "outdir" ].value< string_view >();
    fileName  = *paraTable[ "global" ][ "filename" ].value< string_view >();
    constexpr unsigned int defaultMaxIter = 25;
    constexpr double       defaultEpsilon = 1e-8;  // floating-point number equal threshold
    maxIter = paraTable[ "global" ][ "maxiter" ].value_or( defaultMaxIter );
    epsilon = paraTable[ "global" ][ "outdir" ].value_or( defaultEpsilon );
    assert( maxIter > 0 );
    assert( epsilon > 0 );
    // outdir = "./otfLogs"      # path of the log directorys
    // filename = "galotfa.hdf5" # filename of the log file
    // maxiter = 25              # maximal iteration times during analysis
    // epsilon = 1e-12           # the equal threshold for float numbers

    // orbital log parameters
    orbit = make_unique< otf::orbit >( *paraTable[ "orbit" ].as_table() );

    // parameters of each components
    paraTable.for_each( [ this, paraTable ]( auto& key, auto& value ) {
        if constexpr ( toml::is_key< decltype( key ) > and toml::is_table< decltype( value ) > )
        {
            auto key_view = string_view( key );
            if ( key_view.substr( 0, 9 ) == "component" )
            {
                this->comps[ key_view ] =
                    make_unique< otf::component >( key_view, *value.as_table() );
            }
        }
    } );
}

component::component( string_view& compName, toml::table& compNodeTable ) : compName( compName )
{
    // particle types in this component
    auto typeIDs = compNodeTable[ "types" ];
    if ( toml::array* arr = typeIDs.as_array() )
    {
        // visitation with for_each() helps deal with heterogeneous data
        arr->for_each( [ this ]( auto&& element ) {
            if constexpr ( toml::is_number< decltype( element ) > )
            {
                types.push_back( *element );
            }
        } );
    }

    // period
    period = *compNodeTable[ "period" ].value< double >();

    // recenter parameters
    recenter.enable = *compNodeTable[ "recenter" ][ "enable" ].value< bool >();
    if ( recenter.enable )
    {
        recenter.radius = *compNodeTable[ "recenter" ][ "radius" ].value< double >();
        string_view str = *compNodeTable[ "recenter" ][ "method" ].value< string_view >();
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
        auto iguess = compNodeTable[ "recenter" ][ "iguess" ];
        for ( auto i = 0; i < 3; ++i )
        {
            recenter.initialGuess[ i ] = *iguess[ i ].value< double >();
        }
    }

    // NOTE: the frame parameter is unused at present
    //
    // // frame
    // auto str = *compNodeTable[ "frame" ].value< string_view >();
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
    align.enable = *compNodeTable[ "align" ][ "enable" ].value< bool >();
    if ( align.enable )
    {
        align.radius = *compNodeTable[ "align" ][ "radius" ].value< double >();
    }

    // image
    image.enable = *compNodeTable[ "image" ][ "enable" ].value< bool >();
    if ( image.enable )
    {
        image.halfLength = *compNodeTable[ "image" ][ "halflength" ].value< double >();
        image.binNum     = *compNodeTable[ "image" ][ "binnum" ].value< unsigned int >();
    }

    // bar info parameters
    // A2
    A2.enable = *compNodeTable[ "A2" ][ "enable" ].value< bool >();
    if ( A2.enable )
    {
        A2.rmin = *compNodeTable[ "A2" ][ "rmin" ].value< double >();
        A2.rmax = *compNodeTable[ "A2" ][ "rmax" ].value< double >();
    }
    // bar angle
    barAngle.enable = *compNodeTable[ "barangle" ][ "enable" ].value< bool >();
    if ( barAngle.enable )
    {
        barAngle.rmin = *compNodeTable[ "barangle" ][ "rmin" ].value< double >();
        barAngle.rmax = *compNodeTable[ "barangle" ][ "rmax" ].value< double >();
    }
    // buckling strength
    buckle.enable = *compNodeTable[ "buckle" ][ "enable" ].value< bool >();
    if ( buckle.enable )
    {
        buckle.rmin = *compNodeTable[ "buckle" ][ "rmin" ].value< double >();
        buckle.rmax = *compNodeTable[ "buckle" ][ "rmax" ].value< double >();
    }
}

orbit::orbit( toml::table& orbitNode )
{
    // whether enable orbital logs
    enable = *orbitNode[ "enable" ].value< bool >();

    if ( not enable )
    {
        return;
    }

    // period
    period = *orbitNode[ "period" ].value< double >();

    // particle types to be logged
    auto typeIDs = orbitNode[ "logtypes" ];
    if ( toml::array* arr = typeIDs.as_array() )
    {
        // visitation with for_each() helps deal with heterogeneous data
        arr->for_each( [ this ]( auto&& el ) {
            if constexpr ( toml::is_number< decltype( el ) > )
            {
                logTypes.push_back( *el );
            }
        } );
    }

    auto str = *orbitNode[ "method" ].value< string_view >();
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
        random.enable   = *orbitNode[ "random" ][ "enable" ].value< bool >();
        random.fraction = *orbitNode[ "random" ][ "frac" ].value< double >();
        // recenter.partTypes = orbitNode[ "random" ][ "radius" ].value< bool >();
    }
    else
    {
        idfile = *orbitNode[ "idfile" ].value< string_view >();
    }

    // recenter parameters
    recenter.enable = *orbitNode[ "recenter" ][ "enable" ].value< bool >();
    if ( recenter.enable )
    {
        // particle types to be used as anchors of recenter
        auto typeIDs = orbitNode[ "recenter" ][ "anchorids" ];
        if ( toml::array* arr = typeIDs.as_array() )
        {
            // visitation with for_each() helps deal with heterogeneous data
            arr->for_each( [ this ]( auto&& el ) {
                if constexpr ( toml::is_number< decltype( el ) > )
                {
                    recenter.anchorIds.push_back( *el );
                }
            } );
        }
        recenter.radius = *orbitNode[ "recenter" ][ "radius" ].value< bool >();
        str             = *orbitNode[ "recenter" ][ "method" ].value< string_view >();
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
        auto iguess = orbitNode[ "recenter" ][ "iguess" ];
        for ( auto i = 0; i < 3; ++i )
        {
            recenter.initialGuess[ i ] = *iguess[ i ].value< double >();
        }
    }
}

}  // namespace otf
