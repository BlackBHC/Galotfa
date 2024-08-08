/**
 * @file
 * @brief This file includes a class as a wrapper for statistic functions. At now, mainly the 2D
 * evenly binning statistics for some methods.
 */

#include "../include/statistic.hpp"
#include "../include/myprompt.hpp"
#include <cmath>
#include <cstdlib>
#include <gsl/gsl_statistics_double.h>
#include <memory>
#include <mpi.h>
#include <vector>
using namespace std;

/**
 * @brief 2D binning statistics with chosen method, support count, sum, mean and standard deviation.
 * If need to call the function in MPI mode, it should be called after MPI_Init
 *
 * @param xData: pointing to the first coordinates
 * @param xLowerBound: the lower inclusive limit of the first coordinates to be analyzed
 * @param xUpperBound: the upper exclusive limit of the first coordinates to be analyzed
 * @param xBinNum: binnum of the first coordinate
 * @param yData: pointing to the second coordinates
 * @param yLowerBound: the lower inclusive limit of the second coordinates to be analyzed
 * @param yUpperBound: the upper exclusive limit of the second coordinates to be analyzed
 * @param yBinNum: binnum of the second coordinate
 * @param method: statistic method
 * @param dataNum: number of data points to be analyzed
 * @param data: pointing to target data points
 * @return a unique_ptr pointing to the 1D array of the 2D resutls, in row-major order.
 */
auto statistic::bin2d( double* xData, double xLowerBound, double xUpperBound,
                       const unsigned long& xBinNum, double* yData, double yLowerBound,
                       double yUpperBound, const unsigned long& yBinNum, statistic_method method,
                       const unsigned long& dataNum, double* data ) -> std::unique_ptr< double[] >
{
    int alreadyInitialized;
    MPI_Initialized( &alreadyInitialized );
    if ( !alreadyInitialized )
    {
        MPI_Init( nullptr, nullptr );
    }
    switch ( method )
    {
    case statistic_method::COUNT: {
        return bin2dcount( xData, xLowerBound, xUpperBound, xBinNum, yData, yLowerBound,
                           yUpperBound, yBinNum, dataNum );
    }
    case statistic_method::SUM: {
        return bin2dsum( xData, xLowerBound, xUpperBound, xBinNum, yData, yLowerBound, yUpperBound,
                         yBinNum, dataNum, data );
    }
    case statistic_method::MEAN: {
        return bin2dmean( xData, xLowerBound, xUpperBound, xBinNum, yData, yLowerBound, yUpperBound,
                          yBinNum, dataNum, data );
    }
    case statistic_method::STD: {
        return bin2dstd( xData, xLowerBound, xUpperBound, xBinNum, yData, yLowerBound, yUpperBound,
                         yBinNum, dataNum, data );
    }
    default: {
        ERROR( "Get an unsupported statistic method!" );
        return nullptr;
    }
    }

    if ( !alreadyInitialized )
    {
        MPI_Finalize();
    }
    return nullptr;
}

/**
 * @brief 2D binning statistics for count
 *
 * @param xData: pointing to the first coordinates
 * @param xLowerBound: the lower inclusive limit of the first coordinates to be analyzed
 * @param xUpperBound: the upper exclusive limit of the first coordinates to be analyzed
 * @param xBinNum: binnum of the first coordinate
 * @param yData: pointing to the second coordinates
 * @param yLowerBound: the lower inclusive limit of the second coordinates to be analyzed
 * @param yUpperBound: the upper exclusive limit of the second coordinates to be analyzed
 * @param yBinNum: binnum of the second coordinate
 * @param dataNum: number of data points to be analyzed
 * @return a unique_ptr pointing to the 1D array of the 2D resutls, in row-major order.
 */
auto statistic::bin2dcount( double* xData, double xLowerBound, double xUpperBound,
                            const unsigned long& xBinNum, double* yData, double yLowerBound,
                            double yUpperBound, const unsigned long& yBinNum,
                            const unsigned long& dataNum ) -> std::unique_ptr< double[] >

{
    unsigned long                  idx = 0;
    unsigned long                  idy = 0;
    std::unique_ptr< double[] >    statisticResutls( new double[ xBinNum * yBinNum ]() );
    std::unique_ptr< int[] > const count( new int[ xBinNum * yBinNum ]() );
    std::unique_ptr< int[] > const count_recv( new int[ xBinNum * yBinNum ]() );

    for ( auto i = 0UL; i < dataNum; ++i )
    {
        if ( xData[ i ] < xLowerBound || xData[ i ] >= xUpperBound || yData[ i ] < yLowerBound
             || yData[ i ] >= yUpperBound )
        {
            continue;  // exclude the data points outside the range
        }
        idx = find_index( xLowerBound, xUpperBound, xBinNum, xData[ i ] );
        idy = find_index( yLowerBound, yUpperBound, yBinNum, yData[ i ] );
        ++count[ idx * yBinNum + idy ];
    }

    int rank = -1;
    MPI_Comm_rank( MPI_COMM_WORLD, &rank );
    MPI_Allreduce( MPI_IN_PLACE, count.get(), ( int )xBinNum * yBinNum, MPI_INT, MPI_SUM,
                   MPI_COMM_WORLD );

    for ( auto i = 0U; i < xBinNum * yBinNum; ++i )
    {
        statisticResutls[ i ] = ( double )count[ i ];
    }

    return statisticResutls;
}

/**
 * @brief 2D binning statistics for summation
 *
 * @param xData: pointing to the first coordinates
 * @param xLowerBound: the lower inclusive limit of the first coordinates to be analyzed
 * @param xUpperBound: the upper exclusive limit of the first coordinates to be analyzed
 * @param xBinNum: binnum of the first coordinate
 * @param yData: pointing to the second coordinates
 * @param yLowerBound: the lower inclusive limit of the second coordinates to be analyzed
 * @param yUpperBound: the upper exclusive limit of the second coordinates to be analyzed
 * @param yBinNum: binnum of the second coordinate
 * @param dataNum: number of data points to be analyzed
 * @param data: pointing to target data points
 * @return a unique_ptr pointing to the 1D array of the 2D resutls, in row-major order.
 */
auto statistic::bin2dsum( double* xData, double xLowerBound, double xUpperBound,
                          const unsigned long& xBinNum, double* yData, double yLowerBound,
                          double yUpperBound, const unsigned long& yBinNum,
                          const unsigned long& dataNum,
                          const double*        data ) -> std::unique_ptr< double[] >

{
    unsigned long                     idx = 0;
    unsigned long                     idy = 0;
    std::unique_ptr< double[] >       statisticResutls( new double[ xBinNum * yBinNum ]() );
    std::unique_ptr< double[] > const sum( new double[ xBinNum * yBinNum ]() );

    for ( auto i = 0UL; i < dataNum; ++i )
    {
        if ( xData[ i ] < xLowerBound || xData[ i ] >= xUpperBound || yData[ i ] < yLowerBound
             || yData[ i ] >= yUpperBound )
        {
            continue;  // exclude the data points outside the range
        }
        idx = find_index( xLowerBound, xUpperBound, xBinNum, xData[ i ] );
        idy = find_index( yLowerBound, yUpperBound, yBinNum, yData[ i ] );
        sum[ idx * yBinNum + idy ] += data[ i ];
    }

    for ( auto i = 0U; i < xBinNum * yBinNum; ++i )
    {
        statisticResutls[ i ] = sum[ i ];
    }

    return statisticResutls;
}

/**
 * @brief 2D binning statistics for mean values
 *
 * @param xData: pointing to the first coordinates
 * @param xLowerBound: the lower inclusive limit of the first coordinates to be analyzed
 * @param xUpperBound: the upper exclusive limit of the first coordinates to be analyzed
 * @param xBinNum: binnum of the first coordinate
 * @param yData: pointing to the second coordinates
 * @param yLowerBound: the lower inclusive limit of the second coordinates to be analyzed
 * @param yUpperBound: the upper exclusive limit of the second coordinates to be analyzed
 * @param yBinNum: binnum of the second coordinate
 * @param dataNum: number of data points to be analyzed
 * @param data: pointing to target data points
 * @return a unique_ptr pointing to the 1D array of the 2D resutls, in row-major order.
 */
auto statistic::bin2dmean( double* xData, double xLowerBound, double xUpperBound,
                           const unsigned long& xBinNum, double* yData, double yLowerBound,
                           double yUpperBound, const unsigned long& yBinNum,
                           const unsigned long& dataNum,
                           const double*        data ) -> std::unique_ptr< double[] >

{
    unsigned long                           idx = 0;
    unsigned long                           idy = 0;
    std::unique_ptr< double[] >             statisticResutls( new double[ xBinNum * yBinNum ]() );
    std::unique_ptr< double[] > const       sum( new double[ xBinNum * yBinNum ]() );
    std::unique_ptr< unsigned int[] > const count( new unsigned int[ xBinNum * yBinNum ]() );

    for ( auto i = 0UL; i < dataNum; ++i )
    {
        if ( xData[ i ] < xLowerBound || xData[ i ] >= xUpperBound || yData[ i ] < yLowerBound
             || yData[ i ] >= yUpperBound )
        {
            continue;  // exclude the data points outside the range
        }
        idx = find_index( xLowerBound, xUpperBound, xBinNum, xData[ i ] );
        idy = find_index( yLowerBound, yUpperBound, yBinNum, yData[ i ] );
        ++count[ idx * yBinNum + idy ];
        sum[ idx * yBinNum + idy ] += data[ i ];
    }

    for ( auto i = 0U; i < xBinNum * yBinNum; ++i )
    {
        if ( count[ i ] != 0 )
        {
            statisticResutls[ i ] = sum[ i ] / count[ i ];
        }
        else
        {
            statisticResutls[ i ] = std::nan( "" );
        }
    }

    return statisticResutls;
}

/**
 * @brief 2D binning statistics for standard deviation
 *
 * @param xData: pointing to the first coordinates
 * @param xLowerBound: the lower inclusive limit of the first coordinates to be analyzed
 * @param xUpperBound: the upper exclusive limit of the first coordinates to be analyzed
 * @param xBinNum: binnum of the first coordinate
 * @param yData: pointing to the second coordinates
 * @param yLowerBound: the lower inclusive limit of the second coordinates to be analyzed
 * @param yUpperBound: the upper exclusive limit of the second coordinates to be analyzed
 * @param yBinNum: binnum of the second coordinate
 * @param dataNum: number of data points to be analyzed
 * @param data: pointing to target data points
 * @return a unique_ptr pointing to the 1D array of the 2D resutls, in row-major order.
 */
auto statistic::bin2dstd( double* xData, double xLowerBound, double xUpperBound,
                          const unsigned long& xBinNum, double* yData, double yLowerBound,
                          double yUpperBound, const unsigned long& yBinNum,
                          const unsigned long& dataNum,
                          double*              data ) -> std::unique_ptr< double[] >

{
    unsigned long               idx = 0;
    unsigned long               idy = 0;
    std::unique_ptr< double[] > statisticResutls( new double[ xBinNum * yBinNum ]() );
    std::unique_ptr< std::vector< double >[] > const dataInEachBin(
        new std::vector< double >[ xBinNum * yBinNum ] );

    for ( auto i = 0UL; i < dataNum; ++i )
    {
        if ( xData[ i ] < xLowerBound || xData[ i ] >= xUpperBound || yData[ i ] < yLowerBound
             || yData[ i ] >= yUpperBound )
        {
            continue;  // exclude the data points outside the range
        }
        idx = find_index( xLowerBound, xUpperBound, xBinNum, xData[ i ] );
        idy = find_index( yLowerBound, yUpperBound, yBinNum, yData[ i ] );
        dataInEachBin[ idx * yBinNum + idy ].push_back( data[ i ] );
    }

    for ( auto i = 0U; i < xBinNum * yBinNum; ++i )
    {
        if ( dataInEachBin[ i ].empty() )
        {
            statisticResutls[ i ] = std::nan( "" );
        }
        else
        {
            statisticResutls[ i ] =
                gsl_stats_sd( dataInEachBin[ i ].data(), 1, dataInEachBin[ i ].size() );
        }
    }

    return statisticResutls;
}

/**
 * @brief function that determines the bin a data point should be located, for evenly distributed
 * bins only and there is no boundary check!
 *
 * @param lowerBound: the lower boundary of the data range.
 * @param upperBound: the upper boundary of the data range.
 * @param binNum: the number of bins.
 * @param value: the value of the data point.
 * @return
 */
auto statistic::find_index( double lowerBound, double upperBound, const unsigned long& binNum,
                            double value ) -> unsigned long
{
    return ( value - lowerBound ) / ( upperBound - lowerBound ) * binNum;
}
