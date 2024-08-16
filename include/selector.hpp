/**
 * @file selector.hpp
 * @brief ID selector and reader.
 */
#ifndef SELECTOR_HEADER
#define SELECTOR_HEADER
#include <string>
#include <vector>
namespace otf {

/**
 * @class orbit_selector
 * @brief Data selection for orbital log.
 *
 */
class orbit_selector
{
public:
    void select();

#ifdef DEBUG

#else
private:
#endif
    static auto id_select( const std::vector< unsigned int >& raw,
                           double fraction ) -> std::vector< unsigned int >;
    static auto id_read( const std::string& idFilename ) -> std::vector< unsigned int >;
};

/**
 * @class compDataPtrs
 * @brief Collection of points to the data for component-analysis, aiming for convenience of data
 * selection.
 *
 */
struct compDataPtrs
{
    double* mass       = nullptr;
    double* coordinate = nullptr;
    double* velocity   = nullptr;
};

/**
 * @class orbitDataPtrs
 * @brief Collection of points to the data for orbital log, aiming for convenience of data
 * selection.
 *
 */
struct orbitDataPtrs
{
    double* mass       = nullptr;
    double* coordinate = nullptr;
    double* velocity   = nullptr;
};

}  // namespace otf
#endif
