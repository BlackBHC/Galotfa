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
public:
    unsigned int                count = 0;
    std::vector< unsigned int > id;
    std::vector< double >       mass;
    std::vector< double >       coordinate;
    std::vector< double >       velocity;
};

/**
 * @class orbit_selector
 * @brief Data selection for orbital log.
 *
 */
class orbit_selector
{
public:
    orbit_selector( const runtime_para& para );
    auto select( unsigned int particleNumber, const unsigned int* particleID,
                 const unsigned int* partType, const double* mass, const double* coordinate,
                 const double* velocity ) const -> std::unique_ptr< dataContainer >;

#ifdef DEBUG

#else
private:
#endif
    const runtime_para& para;
    static auto id_sample( const std::vector< unsigned int >& raw, const unsigned int* types,
                           const std::vector< unsigned int >& sampleTypes,
                           double fraction ) -> std::vector< unsigned int >;
    static auto id_read( const std::string& idFilename ) -> std::vector< unsigned int >;
};

}  // namespace otf
#endif
