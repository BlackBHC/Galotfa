/**
 * @file selector.hpp
 * @brief ID selector and reader.
 */
#ifndef SELECTOR_HEADER
#define SELECTOR_HEADER
#include "../include/para.hpp"
#include <memory>
#include <string>
#include <vector>
namespace otf {

/**
 * @class dataContainer
 * @brief Collection of points to the data for component-analysis and orbital log, aiming for
 * convenience of data selection.
 *
 */
struct dataContainer
{
    unsigned int          count = 0;
    std::vector< double > mass;
    std::vector< double > coordinate;
    std::vector< double > velocity;
};

/**
 * @class orbit_selector
 * @brief Data selection for orbital log.
 *
 */
class orbit_selector
{
public:
    orbit_selector( std::unique_ptr< runtime_para >& para );
    auto select( const unsigned int particleNumber, const unsigned int* particleIDs,
                 const unsigned int* particleTypes, const double* masses, const double* coordiantes,
                 const double* velocities ) -> const std::unique_ptr< dataContainer >;

#ifdef DEBUG

#else
private:
#endif
    std::unique_ptr< runtime_para >& para;
    static auto id_sample( const std::vector< unsigned int >& raw, const unsigned int* types,
                           const std::vector< unsigned int >& sampleTypes,
                           double fraction ) -> std::vector< unsigned int >;
    static auto id_read( const std::string& idFilename ) -> std::vector< unsigned int >;
};

}  // namespace otf
#endif
