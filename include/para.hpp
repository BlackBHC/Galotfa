/**
 * @file para.hpp
 * @brief The parameter parser and its utils.
 */

#ifndef PARA_HEADER
#define PARA_HEADER
#include "../include/toml.hpp"
#include "recenter.hpp"
#include <cstdint>
#include <memory>
#include <string_view>
#include <unordered_map>
#include <vector>

namespace otf {

constexpr auto vecDim = 3;

/**
 * @class recenter_para
 * @brief The parameters used for recenter.
 *
 */
struct recenter_para
{
    bool                 enable;
    double               radius;
    double               initialGuess[ vecDim ];
    otf::recenter_method method;
};

/**
 * @class align_para
 * @brief The parameters used for alignment.
 *
 */
struct align_para
{
    bool   enable;
    double radius;
};

/**
 * @class image_para
 * @brief The parameters used for image calculation.
 *
 */
struct image_para
{
    bool         enable;
    double       halfLength;
    unsigned int binNum;
};

/**
 * @class basic_bar_para
 * @brief The parameters used for bar info calculation, A2, Sbar et al.
 *
 */
struct basic_bar_para
{
    bool   enable;
    double rmin;
    double rmax;
};

/**
 * @class orbit_recenter_para
 * @brief The parameters used for recenter in orbital log.
 *
 */
struct orbit_recenter_para : recenter_para
{
    // the anchor type of particles used for recenter
    std::vector< unsigned int > anchorIds;
};

enum class coordinate_frame : std::uint8_t { CYLINDRICAL = 0, SPHERICAL, CARTESIAN };

/**
 * @class component
 * @brief The wrapper of parameter blocks used for each component.
 *
 */
struct component
{
    component( std::string_view& compName, toml::table& compNodeTable );
    std::string                 compName;
    std::vector< unsigned int > types;
    unsigned int                period;
    recenter_para               recenter;
    coordinate_frame            frame;
    align_para                  align;
    image_para                  image;
    basic_bar_para              A2;
    basic_bar_para              barAngle;
    basic_bar_para              buckle;
};

/**
 * @class orbit
 * @brief The wrapper of parameter blocks used for orbital log.
 *
 */
class orbit
{
public:
    orbit( toml::table& orbitNodeTable );
    // method for id log: TXTFILE to use a text file of id list, and RANDOM for random selection
    // according to specified parameters.
    enum class id_selection_method : std::uint8_t { TXTFILE = 0, RANDOM };
    bool                        enable;
    unsigned int                period;
    std::vector< unsigned int > logTypes;
    id_selection_method         method;
    std::string                 idfile   = "not used";
    double                      fraction = -1;
    orbit_recenter_para         recenter;
};

/**
 * @class runtime_para
 * @brief Container of the runtime parameter, designed to be work in each mpi rank.
 *
 */
class runtime_para
{
public:
    runtime_para( const std::string_view& tomlParaFile );
    bool         enableOtf;
    std::string  outputDir;
    std::string  fileName;
    unsigned int maxIter;
    double       epsilon;

    std::unordered_map< std::string, std::unique_ptr< otf::component > > comps;
    std::unique_ptr< otf::orbit >                                        orbit;
};

}  // namespace otf
#endif
