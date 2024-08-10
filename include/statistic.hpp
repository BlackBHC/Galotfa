/**
 * @file statistic.hpp
 * @brief This file includes a class as a wrapper for statistic functions. At now, mainly the 1D/2D
 * evenly binning statistics for limited methods.
 */

#ifndef STATISTIC_HEADER
#define STATISTIC_HEADER
#include <cstdint>
#include <memory>
enum class statistic_method : std::uint8_t { COUNT = 0, MEAN, STD, SUM };

/**
 * @class statistic
 * @brief Wrapper class of the internal statistic system.
 *
 */
class statistic
{
public:
    static auto bin2d( const int mpiRank, const double* xData, const double xLowerBound,
                       const double xUpperBound, const unsigned long& xBinNum, const double* yData,
                       const double yLowerBound, const double yUpperBound,
                       const unsigned long& yBinNum, const statistic_method method,
                       const unsigned long& dataNum,
                       const double*        data = nullptr ) -> std::unique_ptr< double[] >;
    static auto bin1d( const int mpiRank, const double* coord, const double lowerBound,
                       const double upperBound, const unsigned long& binNum,
                       const statistic_method method, const unsigned long& dataNum,
                       const double* data = nullptr ) -> std::unique_ptr< double[] >;

#ifdef DEBUG
public:
#else
private:
#endif
    // function that determines the bin a data point should be located (evenly distributed bins)
    static auto find_index( const double lowerBound, const double upperBound,
                            const unsigned long& binNum, const double value ) -> unsigned long;
    static auto bin2dcount( const int mpiRank, const double* xData, const double xLowerBound,
                            const double xUpperBound, const unsigned long& xBinNum,
                            const double* yData, const double yLowerBound, const double yUpperBound,
                            const unsigned long& yBinNum,
                            const unsigned long& dataNum ) -> std::unique_ptr< double[] >;
    static auto bin2dsum( const int mpiRank, const double* xData, const double xLowerBound,
                          const double xUpperBound, const unsigned long& xBinNum,
                          const double* yData, const double yLowerBound, const double yUpperBound,
                          const unsigned long& yBinNum, const unsigned long& dataNum,
                          const double* data ) -> std::unique_ptr< double[] >;
    static auto bin2dmean( const int mpiRank, const double* xData, const double xLowerBound,
                           const double xUpperBound, const unsigned long& xBinNum,
                           const double* yData, const double yLowerBound, const double yUpperBound,
                           const unsigned long& yBinNum, const unsigned long& dataNum,
                           const double* data ) -> std::unique_ptr< double[] >;
    static auto bin2dstd( const double* xData, const double xLowerBound, const double xUpperBound,
                          const unsigned long& xBinNum, const double* yData,
                          const double yLowerBound, const double yUpperBound,
                          const unsigned long& yBinNum, const unsigned long& dataNum,
                          const double* data ) -> std::unique_ptr< double[] >;
    static auto bin1dcount( const int mpiRank, const double* coord, const double lowerBound,
                            const double upperBound, const unsigned long& binNum,
                            const unsigned long& dataNum ) -> std::unique_ptr< double[] >;
    static auto bin1dsum( const int mpiRank, const double* coord, const double lowerBound,
                          const double upperBound, const unsigned long& binNum,
                          const unsigned long& dataNum,
                          const double*        data = nullptr ) -> std::unique_ptr< double[] >;
    static auto bin1dmean( const int mpiRank, const double* coord, const double lowerBound,
                           const double upperBound, const unsigned long& binNum,
                           const unsigned long& dataNum,
                           const double*        data = nullptr ) -> std::unique_ptr< double[] >;
    static auto bin1dstd( const double* coord, const double lowerBound, const double upperBound,
                          const unsigned long& binNum, const unsigned long& dataNum,
                          const double* data = nullptr ) -> std::unique_ptr< double[] >;
};
#endif
