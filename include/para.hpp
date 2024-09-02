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
    double               radius;                  // enclose radius used for coordinate recenter
    double               initialGuess[ vecDim ];  // initial guess of the coordinate center
    otf::recenter_method method;                  // recenter method
};

/**
 * @class align_para
 * @brief The parameters used for alignment.
 *
 */
struct align_para
{
    bool   enable;
    double radius;  // enclosing radius of the inertia tensor calculation
};

/**
 * @class image_para
 * @brief The parameters used for image calculation.
 *
 */
struct image_para
{
    bool         enable;
    double       halfLength;  // half length of the box size to be plotted
    unsigned int binNum;      // binnum of the image
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
    std::string                 compName;  // name of the component
    std::vector< unsigned int > types;     // particle types in this component
    unsigned int                period;    // analysis period
    recenter_para               recenter;  // parameter of coordinate recenter
    coordinate_frame            frame;     // coordinate frame type
    align_para                  align;     // whether align coordinates with the inertia tensor
    image_para                  image;     // parameter of the spatial image part
    basic_bar_para              A2;        // bar strength parameter
    basic_bar_para              barAngle;  // bar angle parameter
    basic_bar_para              sBuckle;   // buckling strength parameter
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

    bool                enable;                 // enable orbital log
    unsigned int        period;                 // log period
    id_selection_method method;                 // id determination method
    std::string         idfile   = "not used";  // if method is txt file, give the file name
    double              fraction = -1;          // if method is random sample, give the fraction
    std::vector< int >  sampleTypes;            // particle types to be sampled
    orbit_recenter_para recenter;               // whether recenter the coordinate of orbits
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
    bool         enableOtf;  // whether enable on-the-fly analysis
    std::string  outputDir;  // output directory of the logs
    std::string  fileName;   // prefix of the log file
    unsigned int maxIter;    // specify the maximal iteration times
    double       epsilon;    // specify the equal threshold of floating-point numbers

    // hash map of parameter for each component
    std::unordered_map< std::string, std::unique_ptr< otf::component > > comps;
    // parameter pointer of orbital logs
    std::unique_ptr< otf::orbit > orbit;
};

}  // namespace otf
#endif
