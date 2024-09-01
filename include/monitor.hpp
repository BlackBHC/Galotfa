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
#include <vector>

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
    monitor( const std::string_view& tomlParaFile );  // initialize with the toml file name
    ~monitor();
    void one_analysis_api( double time, unsigned int particleNumber, const int* id,
                           const int* partType, const double* mass, const double* coordinate,
                           const double* velocity );

#ifdef DEBUG

#else
private:
#endif
    int           mpiRank;
    int           mpiSize;
    bool          isRootRank;
    unsigned long stepCounter;             // counter of the synchronized time step
    bool          mpiInitialzedByMonitor;  // whether the MPI_init is called by the monitor object
    runtime_para  para;                    // ptr to the runtime paramter
    static constexpr unsigned int orbitPointDim = 7;
    using orbitPoint                            = struct wrapper
    {
        int    particleID;
        double data[ orbitPointDim ];
    };

    // extract the data used for orbital log
    auto id_data_process( double time, unsigned int particleNumber, const int* particleID,
                          const int* particleType, const double* mass, const double* coordinate,
                          const double* velocity ) const -> std::vector< orbitPoint >;

    // extract the data of a component
    void component_data_process();

    void                      data_flush();  // flush the data to the disk
    void                      bar_info();    // bar info calculation
    void                      image();       // image calculation
    std::unique_ptr< h5_out > h5Organizer;   // the ptr to the HDF5 file organizer
};

}  // namespace otf
#endif
