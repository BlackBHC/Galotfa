#include "../include/statistic.hpp"
#include "../include/myprompt.hpp"
#include <cmath>
#include <cstdlib>
#include <cstring>
#include <memory>
#include <mpi.h>
using namespace std;

/**
 * @brief 2D binning statistics with chosen method, support count, sum, mean and standard deviation.
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
auto statistic::bin2d( int mpiRank, const double* xData, const double xLowerBound,
                       const double xUpperBound, const unsigned long& xBinNum, const double* yData,
                       const double yLowerBound, const double yUpperBound,
                       const unsigned long& yBinNum, const statistic_method method,
                       const unsigned long& dataNum, const double* data ) -> unique_ptr< double[] >
{
    switch ( method )
    {
    case statistic_method::COUNT: {
        return bin2dcount( mpiRank, xData, xLowerBound, xUpperBound, xBinNum, yData, yLowerBound,
                           yUpperBound, yBinNum, dataNum );
    }
    case statistic_method::SUM: {
        return bin2dsum( mpiRank, xData, xLowerBound, xUpperBound, xBinNum, yData, yLowerBound,
                         yUpperBound, yBinNum, dataNum, data );
    }
    case statistic_method::MEAN: {
        return bin2dmean( mpiRank, xData, xLowerBound, xUpperBound, xBinNum, yData, yLowerBound,
                          yUpperBound, yBinNum, dataNum, data );
    }
    case statistic_method::STD: {
        return bin2dstd( mpiRank, xData, xLowerBound, xUpperBound, xBinNum, yData, yLowerBound,
                         yUpperBound, yBinNum, dataNum, data );
    }
    default: {
        ERROR( "Get an unsupported statistic method!" );
        return nullptr;
    }
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
auto statistic::bin2dcount( int mpiRank, const double* xData, const double xLowerBound,
                            const double xUpperBound, const unsigned long& xBinNum,
                            const double* yData, const double yLowerBound, const double yUpperBound,
                            const unsigned long& yBinNum,
                            const unsigned long& dataNum ) -> unique_ptr< double[] >

{
    unsigned long                idx = 0;
    unsigned long                idy = 0;
    auto                         statisticResutls( make_unique< double[] >( xBinNum * yBinNum ) );
    auto                         count( make_unique< unsigned int[] >( xBinNum * yBinNum ) );
    unique_ptr< unsigned int[] > countRecv = nullptr;

    if ( mpiRank == 0 )
    {
        countRecv = make_unique< unsigned int[] >( xBinNum * yBinNum );
    }

    for ( auto i = 0UL; i < dataNum; ++i )
    {
        if ( xData[ i ] >= xLowerBound and xData[ i ] < xUpperBound and yData[ i ] >= yLowerBound
             and yData[ i ] < yUpperBound )
        {
            idx = find_index( xLowerBound, xUpperBound, xBinNum, xData[ i ] );
            idy = find_index( yLowerBound, yUpperBound, yBinNum, yData[ i ] );
            ++count[ idx * yBinNum + idy ];
        }
    }

    MPI_Reduce( count.get(), countRecv.get(), xBinNum * yBinNum, MPI_UNSIGNED, MPI_SUM, 0,
                MPI_COMM_WORLD );

    if ( mpiRank == 0 )  // effectively update the results in the root process
    {
        for ( auto i = 0U; i < xBinNum * yBinNum; ++i )
        {
            statisticResutls[ i ] = ( double )countRecv[ i ];
        }
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
auto statistic::bin2dsum( int mpiRank, const double* xData, const double xLowerBound,
                          const double xUpperBound, const unsigned long& xBinNum,
                          const double* yData, const double yLowerBound, const double yUpperBound,
                          const unsigned long& yBinNum, const unsigned long& dataNum,
                          const double* data ) -> unique_ptr< double[] >

{
    unsigned long          idx = 0;
    unsigned long          idy = 0;
    auto                   statisticResutls( make_unique< double[] >( xBinNum * yBinNum ) );
    auto                   sum( make_unique< double[] >( xBinNum * yBinNum ) );
    unique_ptr< double[] > sumRecv = nullptr;

    if ( mpiRank == 0 )
    {
        sumRecv = make_unique< double[] >( xBinNum * yBinNum );
    }

    for ( auto i = 0UL; i < dataNum; ++i )
    {
        if ( xData[ i ] >= xLowerBound and xData[ i ] < xUpperBound and yData[ i ] >= yLowerBound
             and yData[ i ] < yUpperBound )
        {
            idx = find_index( xLowerBound, xUpperBound, xBinNum, xData[ i ] );
            idy = find_index( yLowerBound, yUpperBound, yBinNum, yData[ i ] );
            sum[ idx * yBinNum + idy ] += data[ i ];
        }
    }

    MPI_Reduce( sum.get(), sumRecv.get(), xBinNum * yBinNum, MPI_DOUBLE, MPI_SUM, 0,
                MPI_COMM_WORLD );

    if ( mpiRank == 0 )  // effectively update the results in the root process
    {
        for ( auto i = 0U; i < xBinNum * yBinNum; ++i )
        {
            statisticResutls[ i ] = sumRecv[ i ];
        }
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
auto statistic::bin2dmean( int mpiRank, const double* xData, const double xLowerBound,
                           const double xUpperBound, const unsigned long& xBinNum,
                           const double* yData, const double yLowerBound, const double yUpperBound,
                           const unsigned long& yBinNum, const unsigned long& dataNum,
                           const double* data ) -> unique_ptr< double[] >
{
    unsigned long idx = 0;
    unsigned long idy = 0;
    auto          statisticResutls( make_unique< double[] >( xBinNum * yBinNum ) );
    auto          sum( make_unique< double[] >( xBinNum * yBinNum ) );
    auto          count( make_unique< unsigned int[] >( xBinNum * yBinNum ) );
    auto          sumRecv( make_unique< double[] >( xBinNum * yBinNum ) );
    auto          countRecv( make_unique< unsigned int[] >( xBinNum * yBinNum ) );

    if ( mpiRank == 0 )
    {
        sumRecv   = make_unique< double[] >( xBinNum * yBinNum );
        countRecv = make_unique< unsigned int[] >( xBinNum * yBinNum );
    }

    for ( auto i = 0UL; i < dataNum; ++i )
    {
        if ( xData[ i ] >= xLowerBound and xData[ i ] < xUpperBound and yData[ i ] >= yLowerBound
             and yData[ i ] < yUpperBound )
        {
            idx = find_index( xLowerBound, xUpperBound, xBinNum, xData[ i ] );
            idy = find_index( yLowerBound, yUpperBound, yBinNum, yData[ i ] );
            ++count[ idx * yBinNum + idy ];
            sum[ idx * yBinNum + idy ] += data[ i ];
        }
    }

    MPI_Reduce( count.get(), countRecv.get(), xBinNum * yBinNum, MPI_UNSIGNED, MPI_SUM, 0,
                MPI_COMM_WORLD );
    MPI_Reduce( sum.get(), sumRecv.get(), xBinNum * yBinNum, MPI_DOUBLE, MPI_SUM, 0,
                MPI_COMM_WORLD );

    if ( mpiRank == 0 )
    {
        for ( auto i = 0U; i < xBinNum * yBinNum; ++i )
        {
            if ( count[ i ] != 0 )
            {
                statisticResutls[ i ] = sumRecv[ i ] / countRecv[ i ];
            }
            else
            {
                statisticResutls[ i ] = nan( "" );
            }
        }
    }

    return statisticResutls;
}

/**
 * @brief 2D binning statistics for standard deviation, without Bessel correction.
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
auto statistic::bin2dstd( int mpiRank, const double* xData, const double xLowerBound,
                          const double xUpperBound, const unsigned long& xBinNum,
                          const double* yData, const double yLowerBound, const double yUpperBound,
                          const unsigned long& yBinNum, const unsigned long& dataNum,
                          const double* data ) -> unique_ptr< double[] >

{
    unsigned long idx = 0;
    unsigned long idy = 0;
    auto          statisticResutls( make_unique< double[] >( xBinNum * yBinNum ) );
    auto          sum( make_unique< double[] >( xBinNum * yBinNum ) );
    auto          count( make_unique< unsigned int[] >( xBinNum * yBinNum ) );
    auto          sumRecv( make_unique< double[] >( xBinNum * yBinNum ) );
    auto          countRecv( make_unique< unsigned int[] >( xBinNum * yBinNum ) );

    if ( mpiRank == 0 )
    {
        sumRecv   = make_unique< double[] >( xBinNum * yBinNum );
        countRecv = make_unique< unsigned int[] >( xBinNum * yBinNum );
    }

    for ( auto i = 0UL; i < dataNum; ++i )
    {
        if ( xData[ i ] >= xLowerBound and xData[ i ] < xUpperBound and yData[ i ] >= yLowerBound
             and yData[ i ] < yUpperBound )
        {
            idx = find_index( xLowerBound, xUpperBound, xBinNum, xData[ i ] );
            idy = find_index( yLowerBound, yUpperBound, yBinNum, yData[ i ] );
            ++count[ idx * yBinNum + idy ];
            sum[ idx * yBinNum + idy ] += data[ i ];
        }
    }

    MPI_Reduce( count.get(), countRecv.get(), xBinNum * yBinNum, MPI_UNSIGNED, MPI_SUM, 0,
                MPI_COMM_WORLD );
    MPI_Reduce( sum.get(), sumRecv.get(), xBinNum * yBinNum, MPI_DOUBLE, MPI_SUM, 0,
                MPI_COMM_WORLD );

    if ( mpiRank == 0 )
    {
        for ( auto i = 0U; i < xBinNum * yBinNum; ++i )
        {
            if ( count[ i ] != 0 )
            {
                statisticResutls[ i ] = sumRecv[ i ] / countRecv[ i ];
            }
            else
            {
                statisticResutls[ i ] = nan( "" );
            }
        }
    }

    MPI_Bcast( statisticResutls.get(), xBinNum * yBinNum, MPI_DOUBLE, 0, MPI_COMM_WORLD );
    std::memset( sum.get(), 0, sizeof( double ) * xBinNum * yBinNum );  // reset the sum to 0
    for ( auto i = 0UL; i < dataNum; ++i )
    {
        if ( xData[ i ] >= xLowerBound and xData[ i ] < xUpperBound and yData[ i ] >= yLowerBound
             and yData[ i ] < yUpperBound )
        {
            idx = find_index( xLowerBound, xUpperBound, xBinNum, xData[ i ] );
            idy = find_index( yLowerBound, yUpperBound, yBinNum, yData[ i ] );
            sum[ idx * yBinNum + idy ] += ( data[ i ] - statisticResutls[ idx * yBinNum + idy ] )
                                          * ( data[ i ] - statisticResutls[ idx * yBinNum + idy ] );
        }
    }

    MPI_Reduce( sum.get(), sumRecv.get(), xBinNum * yBinNum, MPI_DOUBLE, MPI_SUM, 0,
                MPI_COMM_WORLD );

    if ( mpiRank == 0 )
    {
        for ( auto i = 0U; i < xBinNum * yBinNum; ++i )
        {
            if ( count[ i ] != 0 )
            {
                statisticResutls[ i ] = sqrt( sumRecv[ i ] / countRecv[ i ] );
            }
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

/**
 * @brief Similar to bin2d but for 1D case.
 *
 * @param coord: pointing to the coordinates
 * @param lowerBound: the lower inclusive limit of the first coordinates to be analyzed
 * @param upperBound: the upper exclusive limit of the first coordinates to be analyzed
 * @param binNum: binnum of the first coordinate
 * @param method: statistic method
 * @param dataNum: number of data points to be analyzed
 * @param data: pointing to target data points
 * @return a unique_ptr pointing to the 1D array of resutls.
 */
auto statistic::bin1d( int mpiRank, const double* coord, const double lowerBound,
                       const double upperBound, const unsigned long& binNum,
                       const statistic_method method, const unsigned long& dataNum,
                       const double* data ) -> std::unique_ptr< double[] >
{
    switch ( method )
    {
    case statistic_method::COUNT: {
        return bin1dcount( mpiRank, coord, lowerBound, upperBound, binNum, dataNum );
    }
    case statistic_method::SUM: {
        return bin1dsum( mpiRank, coord, lowerBound, upperBound, binNum, dataNum, data );
    }
    case statistic_method::MEAN: {
        return bin1dmean( mpiRank, coord, lowerBound, upperBound, binNum, dataNum, data );
    }
    case statistic_method::STD: {
        return bin1dstd( mpiRank, coord, lowerBound, upperBound, binNum, dataNum, data );
    }
    default: {
        ERROR( "Get an unsupported statistic method!" );
        return nullptr;
    }
    }

    return nullptr;
}

/**
 * @brief similar to bin2dcount but for 1D case.
 */
auto statistic::bin1dcount( int mpiRank, const double* coord, const double lowerBound,
                            const double upperBound, const unsigned long& binNum,
                            const unsigned long& dataNum ) -> std::unique_ptr< double[] >
{
    unsigned long idx = 0;
    auto          statisticResutls( make_unique< double[] >( binNum ) );
    auto          count( make_unique< unsigned int[] >( binNum ) );
    auto          countRecv( make_unique< unsigned int[] >( binNum ) );

    if ( mpiRank == 0 )
    {
        countRecv = make_unique< unsigned int[] >( binNum );
    }

    for ( auto i = 0UL; i < dataNum; ++i )
    {
        if ( coord[ i ] >= lowerBound and coord[ i ] < upperBound )
        {
            idx = find_index( lowerBound, upperBound, binNum, coord[ i ] );
            ++count[ idx ];
        }
    }

    MPI_Reduce( count.get(), countRecv.get(), binNum, MPI_UNSIGNED, MPI_SUM, 0, MPI_COMM_WORLD );

    if ( mpiRank == 0 )  // effectively update the results in the root process
    {
        for ( auto i = 0U; i < binNum; ++i )
        {
            statisticResutls[ i ] = ( double )countRecv[ i ];
        }
    }
    return statisticResutls;
}

/**
 * @brief similar to bin2dsum but for 1D case.
 */
auto statistic::bin1dsum( int mpiRank, const double* coord, const double lowerBound,
                          const double upperBound, const unsigned long& binNum,
                          const unsigned long& dataNum,
                          const double*        data ) -> std::unique_ptr< double[] >
{
    unsigned long idx = 0;
    auto          statisticResutls( make_unique< double[] >( binNum ) );
    auto          sum( make_unique< double[] >( binNum ) );
    auto          sumRecv( make_unique< double[] >( binNum ) );

    if ( mpiRank == 0 )
    {
        sumRecv = make_unique< double[] >( binNum );
    }

    for ( auto i = 0UL; i < dataNum; ++i )
    {
        if ( coord[ i ] >= lowerBound and coord[ i ] < upperBound )
        {
            idx = find_index( lowerBound, upperBound, binNum, coord[ i ] );
            sum[ idx ] += data[ i ];
        }
    }

    MPI_Reduce( sum.get(), sumRecv.get(), binNum, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD );

    if ( mpiRank == 0 )  // effectively update the results in the root process
    {
        for ( auto i = 0U; i < binNum; ++i )
        {
            statisticResutls[ i ] = ( double )sumRecv[ i ];
        }
    }
    return statisticResutls;
}

/**
 * @brief similar to bin2dmean but for 1D case.
 */
auto statistic::bin1dmean( int mpiRank, const double* coord, const double lowerBound,
                           const double upperBound, const unsigned long& binNum,
                           const unsigned long& dataNum,
                           const double*        data ) -> std::unique_ptr< double[] >
{
    unsigned long idx = 0;
    auto          statisticResutls( make_unique< double[] >( binNum ) );
    auto          sum( make_unique< double[] >( binNum ) );
    auto          count( make_unique< unsigned int[] >( binNum ) );
    auto          sumRecv( make_unique< double[] >( binNum ) );
    auto          countRecv( make_unique< unsigned int[] >( binNum ) );

    if ( mpiRank == 0 )
    {
        sumRecv   = make_unique< double[] >( binNum );
        countRecv = make_unique< unsigned int[] >( binNum );
    }

    for ( auto i = 0UL; i < dataNum; ++i )
    {
        if ( coord[ i ] >= lowerBound and coord[ i ] < upperBound )
        {
            idx = find_index( lowerBound, upperBound, binNum, coord[ i ] );
            ++count[ idx ];
            sum[ idx ] += data[ i ];
        }
    }

    MPI_Reduce( count.get(), countRecv.get(), binNum, MPI_UNSIGNED, MPI_SUM, 0, MPI_COMM_WORLD );
    MPI_Reduce( sum.get(), sumRecv.get(), binNum, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD );

    if ( mpiRank == 0 )  // effectively update the results in the root process
    {
        for ( auto i = 0U; i < binNum; ++i )
        {
            if ( count[ i ] != 0 )
            {
                statisticResutls[ i ] = sumRecv[ i ] / countRecv[ i ];
            }
            else
            {
                statisticResutls[ i ] = nan( "" );
            }
        }
    }
    return statisticResutls;
}

/**
 * @brief similar to bin2dstd but for 1D case.
 */
auto statistic::bin1dstd( int mpiRank, const double* coord, const double lowerBound,
                          const double upperBound, const unsigned long& binNum,
                          const unsigned long& dataNum,
                          const double*        data ) -> std::unique_ptr< double[] >
{
    unsigned long idx = 0;
    auto          statisticResutls( make_unique< double[] >( binNum ) );
    auto          sum( make_unique< double[] >( binNum ) );
    auto          count( make_unique< unsigned int[] >( binNum ) );
    auto          sumRecv( make_unique< double[] >( binNum ) );
    auto          countRecv( make_unique< unsigned int[] >( binNum ) );

    if ( mpiRank == 0 )
    {
        sumRecv   = make_unique< double[] >( binNum );
        countRecv = make_unique< unsigned int[] >( binNum );
    }

    for ( auto i = 0UL; i < dataNum; ++i )
    {
        if ( coord[ i ] >= lowerBound and coord[ i ] < upperBound )
        {
            idx = find_index( lowerBound, upperBound, binNum, coord[ i ] );
            ++count[ idx ];
            sum[ idx ] += data[ i ];
        }
    }

    MPI_Reduce( count.get(), countRecv.get(), binNum, MPI_UNSIGNED, MPI_SUM, 0, MPI_COMM_WORLD );
    MPI_Reduce( sum.get(), sumRecv.get(), binNum, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD );

    if ( mpiRank == 0 )  // effectively update the results in the root process
    {
        for ( auto i = 0U; i < binNum; ++i )
        {
            if ( count[ i ] != 0 )
            {
                statisticResutls[ i ] = sumRecv[ i ] / countRecv[ i ];
            }
            else
            {
                statisticResutls[ i ] = nan( "" );
            }
        }
    }

    MPI_Bcast( statisticResutls.get(), binNum, MPI_DOUBLE, 0, MPI_COMM_WORLD );
    std::memset( sum.get(), 0, sizeof( double ) * binNum );  // reset the sum to 0
    for ( auto i = 0UL; i < dataNum; ++i )
    {
        if ( coord[ i ] >= lowerBound and coord[ i ] < upperBound )
        {
            idx = find_index( lowerBound, upperBound, binNum, coord[ i ] );
            sum[ idx ] +=
                ( data[ i ] - statisticResutls[ idx ] ) * ( data[ i ] - statisticResutls[ idx ] );
        }
    }

    MPI_Reduce( sum.get(), sumRecv.get(), binNum, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD );

    if ( mpiRank == 0 )
    {
        for ( auto i = 0U; i < binNum; ++i )
        {
            if ( count[ i ] != 0 )
            {
                statisticResutls[ i ] = sqrt( sumRecv[ i ] / countRecv[ i ] );
            }
        }
    }

    return statisticResutls;
}
