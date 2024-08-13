/**
 * @file barinfo.hpp
 * @brief Wrapper class of bar quantification functions.
 */

#ifndef BARINFO_HEADER
#define BARINFO_HEADER
/**
 * @class bar_info
 * @brief A0, A2, Sbar, Sbuckle, bar length (to be implemented), bar ellipticity (to be implemented)
 *
 */
class bar_info
{
public:
    static auto A0( unsigned int partNum, const double* mass ) -> double;
    static auto A2( unsigned int partNum, const double* mass, const double* phi ) -> double;
    static auto Sbar( unsigned int partNum, const double* mass, const double* phi ) -> double;
    static auto Sbuckle( unsigned int partNum, const double* mass, const double* phi,
                         const double* zed ) -> double;
};
#endif
