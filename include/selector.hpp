/**
 * @file selector.hpp
 * @brief ID selector and reader.
 */
#ifndef SELECTOR_HEADER
#define SELECTOR_HEADER
#include <string>
#include <vector>
namespace otf {

class id_organizer
{
public:
    static auto select( const std::vector< unsigned int >& raw,
                        double fraction ) -> std::vector< unsigned int >;
    static auto read( const std::string& idFilename ) -> std::vector< unsigned int >;
};

}  // namespace otf
#endif
