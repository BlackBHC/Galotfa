/**
 * @file para.hpp
 * @brief The parameter parser and its utils.
 */

#ifndef PARA_HEADER
#define PARA_HEADER
#include "../include/toml.hpp"
#include "recenter.hpp"
#include <cstdint>
#include <string_view>
#include <vector>

/**
 * @class runtime_para
 * @brief Container of the runtime parameter, designed to be work in each mpi rank.
 *
 */
class runtime_para
{
public:
    runtime_para( const std::string_view& tomlParaFile );
    toml::table paraTable;

private:
    void read_one_by_one();
};

/**
 * @class recenter_para
 * @brief The parameters used for recenter.
 *
 */
struct recenter_para
{
    bool   enable;
    double radius;
    double initialGuess[ 3 ];

    // From include/recenter.hpp:
    // enum class recenter_method : std::uint8_t { COM = 0, MBP = 1 };
    recenter_method method;
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
 * @class id_random_selection_para
 * @brief The parameters used for random selection of particle ids.
 *
 */
struct id_random_selection_para
{
    bool                        enable;
    double                      fraction;
    std::vector< unsigned int > partTypes;
};

enum class coordinate_frame : std::uint8_t { CYLINDRICAL = 0, SPHERICAL, CARTESIAN };

/**
 * @class component
 * @brief The wrapper of parameter blocks used for each component.
 *
 */
class component
{
public:
    component( std::string_view& compName, toml::table& para );
    auto operator[]( std::string_view& key );

private:
    std::string_view compName;
    unsigned int     period;
    recenter_para    recenterPara;
    coordinate_frame frame;
    align_para       alignPara;
    image_para       imagePara;
    basic_bar_para   A2Para;
    basic_bar_para   barAnglePara;
    basic_bar_para   bucklePara;
};

/**
 * @class orbit
 * @brief The wrapper of parameter blocks used for orbital log.
 *
 */
class orbit
{
public:
    unsigned int period;
};
#endif
