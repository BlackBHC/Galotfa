#ifndef STATISTIC_HEADER
#define STATISTIC_HEADER
#include <cstdint>
#include <memory>
enum class statistic_method : std::uint8_t { COUNT = 0, MEAN, STD, SUM };
class statistic
{
public:
    static auto bin2d( int rank, const double* xData, const double xLowerBound,
                       const double xUpperBound, const unsigned long& xBinNum, const double* yData,
                       const double yLowerBound, const double yUpperBound,
                       const unsigned long& yBinNum, statistic_method method,
                       const unsigned long& dataNum,
                       const double*        data = nullptr ) -> std::unique_ptr< double[] >;

#ifdef DEBUG
public:
#else
private:
#endif
    // function that determines the bin a data point should be located (evenly distributed bins)
    static auto find_index( const double lowerBound, const double upperBound,
                            const unsigned long& binNum, const double value ) -> unsigned long;
    static auto bin2dcount( int rank, const double* xData, const double xLowerBound,
                            const double xUpperBound, const unsigned long& xBinNum,
                            const double* yData, const double yLowerBound, const double yUpperBound,
                            const unsigned long& yBinNum,
                            const unsigned long& dataNum ) -> std::unique_ptr< double[] >;
    static auto bin2dsum( int rank, const double* xData, const double xLowerBound,
                          const double xUpperBound, const unsigned long& xBinNum,
                          const double* yData, const double yLowerBound, const double yUpperBound,
                          const unsigned long& yBinNum, const unsigned long& dataNum,
                          const double* data ) -> std::unique_ptr< double[] >;
    static auto bin2dmean( int rank, const double* xData, const double xLowerBound,
                           const double xUpperBound, const unsigned long& xBinNum,
                           const double* yData, const double yLowerBound, const double yUpperBound,
                           const unsigned long& yBinNum, const unsigned long& dataNum,
                           const double* data ) -> std::unique_ptr< double[] >;
    static auto bin2dstd( int rank, const double* xData, const double xLowerBound,
                          const double xUpperBound, const unsigned long& xBinNum,
                          const double* yData, const double yLowerBound, const double yUpperBound,
                          const unsigned long& yBinNum, const unsigned long& dataNum,
                          const double* data ) -> std::unique_ptr< double[] >;
};
#endif
