#ifndef STATISTIC_HEADER
#define STATISTIC_HEADER
#include <cstdint>
#include <memory>
enum class statistic_method : std::uint8_t { COUNT = 0, MEAN, STD, SUM };
class statistic
{
public:
    static auto bin2d( double* xData, double xLowerBound, double xUpperBound,
                       const unsigned long& xBinNum, double* yData, double yLowerBound,
                       double yUpperBound, const unsigned long& yBinNum, statistic_method method,
                       const unsigned long& dataNum,
                       double*              data = nullptr ) -> std::unique_ptr< double[] >;

#ifdef DEBUG
public:
#else
private:
#endif
    // function that determines the bin a data point should be located (evenly distributed bins)
    static auto find_index( double lowerBound, double upperBound, const unsigned long& binNum,
                            double value ) -> unsigned long;
    static auto bin2dcount( double* xData, double xLowerBound, double xUpperBound,
                            const unsigned long& xBinNum, double* yData, double yLowerBound,
                            double yUpperBound, const unsigned long& yBinNum,
                            const unsigned long& dataNum ) -> std::unique_ptr< double[] >;
    static auto bin2dsum( double* xData, double xLowerBound, double xUpperBound,
                          const unsigned long& xBinNum, double* yData, double yLowerBound,
                          double yUpperBound, const unsigned long& yBinNum,
                          const unsigned long& dataNum,
                          const double*        data ) -> std::unique_ptr< double[] >;
    static auto bin2dmean( double* xData, double xLowerBound, double xUpperBound,
                           const unsigned long& xBinNum, double* yData, double yLowerBound,
                           double yUpperBound, const unsigned long& yBinNum,
                           const unsigned long& dataNum,
                           const double*        data ) -> std::unique_ptr< double[] >;
    static auto bin2dstd( double* xData, double xLowerBound, double xUpperBound,
                          const unsigned long& xBinNum, double* yData, double yLowerBound,
                          double yUpperBound, const unsigned long& yBinNum,
                          const unsigned long& dataNum,
                          double*              data ) -> std::unique_ptr< double[] >;
};
#endif
