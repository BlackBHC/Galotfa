#include "../include/para.hpp"

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
