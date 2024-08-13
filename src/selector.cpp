#include "../include/selector.hpp"
#include "../include/myprompt.hpp"
#include <algorithm>
#include <cstdlib>
#include <random>
#include <vector>
using namespace std;

/**
 * @brief Select particle ids with a specified fraction.
 *
 * @param raw raw id list
 * @param fraction fraction
 * @return a std::vector<unsigned int> of the selected id list
 */
vector< unsigned int > id_selector( const vector< unsigned int >& raw, const double fraction )
{
    if ( fraction <= 0 or fraction > 1 )
    {
        ERROR( "Try to select a illegal fraction: [%lf]", fraction );
    }
    else if ( fraction == 1 )
    {
        return raw;
    }

    size_t                 selectNum = ( size_t )( raw.size() * fraction );
    vector< unsigned int > res;
    sample( raw.begin(), raw.end(), back_inserter( res ), selectNum, mt19937{ random_device{}() } );
    return res;
}
