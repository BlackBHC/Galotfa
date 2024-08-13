/**
 * @file recenter.hpp
 * @brief Some utilities for recenter of N-body system.
 */
#ifndef RECENTER_HEADER
#define RECENTER_HEADER
#include <cstdint>
#include <memory>

enum class recenter_method : std::uint8_t { COM = 0, MBP = 1 };

/**
 * @class recenter
 * @brief Wrapper class of the recenter APIs.
 *
 */
class recenter
{
public:
    static auto run( recenter_method method, const unsigned int& partNum,
                     const double* coordinate = nullptr, const double* mass = nullptr,
                     const double* potential = nullptr,
                     double        rangeSize = -1 ) -> std::unique_ptr< double[] >;

#ifdef DEBUG

#else
private:
#endif
    static auto center_of_mass( const double* mass, const double* coordinates,
                                const unsigned int& partNum,
                                double              rangeSize ) -> std::unique_ptr< double[] >;
    static auto most_bound_particle( const double* potential, const double* coordinates,
                                     const unsigned int& partNum ) -> std::unique_ptr< double[] >;
};
#endif
