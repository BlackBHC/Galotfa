/**
 * @file selector.hpp
 * @brief ID selector and reader.
 */
#ifndef SELECTOR_HEADER
#define SELECTOR_HEADER
#include <string>
#include <vector>
auto id_selector( const std::vector< int >& raw, const double fraction ) -> std::vector< int >;
auto id_reader( const std::string& idFilename ) -> std::vector< int >;
#endif
