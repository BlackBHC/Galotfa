/**
 * @file selector.hpp
 * @brief ID selector and reader.
 */
#ifndef SELECTOR_HEADER
#define SELECTOR_HEADER
#include <string>
#include <vector>
class id_organizer
{
public:
    static auto select( const std::vector< int >& raw,
                        const double              fraction ) -> std::vector< int >;
    static auto read( const std::string& idFilename ) -> std::vector< int >;
};
#endif
