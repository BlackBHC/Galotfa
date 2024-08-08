#define DEBUG
#include "../include/h5out.hpp"
#include "H5Tpublic.h"
#include <hdf5.h>

int main()
{
    int testData[ 30 * 3 ];
    int count = 0;
    for ( int i = 0; i < 30; ++i )
        for ( int j = 0; j < 3; ++j )
            testData[ i * 3 + j ] = count++;


    h5_out organizer( "./test_log/", "galotfa.hdf5" );
    organizer.create_dataset_in_group( "A2", "component1", { 3 }, H5T_NATIVE_INT );
    for ( int i = 0; i < 30; ++i )
    {
        organizer.flush_single_block( "component1", "A2", testData + i * 3 );
    }

    double mockImage[ 4 * 4 ];
    for ( int i = 0; i < 16; ++i )
    {
        mockImage[ i ] = i + 0.314;
    }
    organizer.create_dataset_in_group( "Image", "component1", { 4, 4 }, H5T_NATIVE_DOUBLE );
    for ( int i = 0; i < 4; ++i )
    {
        organizer.flush_single_block( "component1", "Image", mockImage );
    }

    return 0;
}
