/**
 * 2D Heat Equation Simulation with VTK file output for ParaView
 *
 * Solves:   dT/dt = alpha * (d2T/dx2 + d2T/dy2)
 *
 * Uses Eigen for the finite-difference computation and writes one
 * VTK ImageData (.vti) file per output step, plus a ParaView Data (.pvd)
 * collection file that lists every step with its simulation time.
 *
 * Field data is encoded as inline base64 binary (format="binary") which is
 * ~4x more compact than ASCII and parsed faster by ParaView.
 *
 * Open heat_output/heat.pvd in ParaView to get the full animated time series.
 *
 * Build:  see CMakeLists.txt
 * Run:    ./heat_equation  [output_dir]
 */

#include <string>
#include <iostream>
#include <cstdlib>

#include "heat_solver.h"


#ifdef _WIN32
#include <direct.h>
#define MKDIR(p) _mkdir(p)
#else
#include <sys/stat.h>
#define MKDIR(p) mkdir(p, 0755)
#endif


int main(int argc, char* argv[])
{
    const std::string out_dir = (argc > 1) ? argv[1] : "heat_output";
    const int simulated_delay_ms = (argc > 2) ? std::max(0, std::atoi(argv[2])) : 5000;

    if (MKDIR(out_dir.c_str()) != 0 && errno != EEXIST)
    {
        std::cerr << "ERROR: could not create directory " << out_dir << "\n";
        return EXIT_FAILURE;
    }

    HeatSolver solver(200, 0.005, 0.4);
    solver.set_simulated_delay_ms(simulated_delay_ms);
    
    if (!solver.run_simulation(out_dir))
    {
        std::cerr << "ERROR: simulation failed\n";
        return EXIT_FAILURE;
    }

}