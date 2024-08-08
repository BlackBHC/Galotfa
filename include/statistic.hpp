#ifndef STATISTIC_HEADER
#define STATISTIC_HEADER
#include <memory>
enum statistic_method { COUNT, MEAN, STD, SUM };
class statistic
{
public:
    static std::unique_ptr< double[] >
    bin2d( double* xData, double xLowerBound, double xUpperBound, const unsigned long& xBinNum,
           double* yData, double yLowerBound, double yUpperBound, const unsigned long& yBinNum,
           statistic_method method, const unsigned long& dataNum, double* data = nullptr );

#ifdef DEBUG
public:
#else
private:
#endif
    // function that determines the bin a data point should be located (evenly distributed bins)
    static unsigned long               find_index( double lowerBound, double upperBound,
                                                   const unsigned long& binNum, double value );
    static std::unique_ptr< double[] > bin2dcount( double* xData, double xLowerBound,
                                                   double xUpperBound, const unsigned long& xBinNum,
                                                   double* yData, double yLowerBound,
                                                   double yUpperBound, const unsigned long& yBinNum,
                                                   const unsigned long& dataNum );
    static std::unique_ptr< double[] > bin2dsum( double* xData, double xLowerBound,
                                                 double xUpperBound, const unsigned long& xBinNum,
                                                 double* yData, double yLowerBound,
                                                 double yUpperBound, const unsigned long& yBinNum,
                                                 const unsigned long& dataNum, const double* data );
    static std::unique_ptr< double[] >
    bin2dmean( double* xData, double xLowerBound, double xUpperBound, const unsigned long& xBinNum,
               double* yData, double yLowerBound, double yUpperBound, const unsigned long& yBinNum,
               const unsigned long& dataNum, const double* data );
    static std::unique_ptr< double[] > bin2dstd( double* xData, double xLowerBound,
                                                 double xUpperBound, const unsigned long& xBinNum,
                                                 double* yData, double yLowerBound,
                                                 double yUpperBound, const unsigned long& yBinNum,
                                                 const unsigned long& dataNum, double* data );
};
#endif
