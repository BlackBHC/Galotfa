/**
 * @file para.hpp
 * @brief The parameter parser.
 */

#ifndef PARA_HEADER
#define PARA_HEADER
#include "../include/toml.hpp"
#ifdef DEBUG  // to avoid repeated definition or macros, must be called after toml.hpp
#include "../include/myprompt.hpp"
#endif
#include <string_view>

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
#endif
