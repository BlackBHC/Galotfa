/**
 * @file monitor.hpp
 * @brief The monitor/server of on-the-fly analysis.
 */

#ifndef MONITOR_HEADER
#define MONITOR_HEADER
#include "../include/h5out.hpp"
#include "../include/para.hpp"
#include <memory>
#include <string_view>

namespace otf {

/**
 * @class monitor
 * @brief The main server of the on-the-fly analysis, which organize the other modules to work
 * together
 *
 */
class monitor
{
public:
    monitor( const std::string_view& tomlParaFile );
    ~monitor();
    void one_analysis_api( const unsigned int* id, const unsigned int* partType, const double* mass,
                           const double* coordinate, const double* velocity );

#ifdef DEBUG

#else
private:
#endif
    int                             mpiRank;
    bool                            isRootRank;
    unsigned long                   stepCounter;
    bool                            mpiInitialzedByMonitor;
    std::unique_ptr< runtime_para > para;
    void                            id_data_process();
    void                            component_data_process();
    void                            data_flush();
    void                            bar_info();
    void                            image();
    std::unique_ptr< h5_out >       h5Orgnizer;
};

}  // namespace otf
#endif
