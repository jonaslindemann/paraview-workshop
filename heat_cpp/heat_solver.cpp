#include "heat_solver.h"

#include <chrono>
#include <iostream>
#include <iomanip>
#include <thread>
#include <vector>

#include "vtk_writer.h"

HeatSolver::HeatSolver(int gridSize, double alpha, double dt)
    : m_nx(gridSize), m_ny(gridSize), m_alpha(alpha)
{
    const double dx = m_lx / (m_nx - 1);
    const double dy = m_ly / (m_ny - 1);

    const double dt_max = (std::min(dx, dy) * std::min(dx, dy)) / (2.0 * m_alpha);
    if (dt >= dt_max) {
        std::cerr << "WARNING: dt=" << dt << " exceeds the stability limit of " << dt_max << "\n";
    }
}

void HeatSolver::set_simulated_delay_ms(int delayMs)
{
    m_simulatedDelayMs = (delayMs < 0) ? 0 : delayMs;
}

void HeatSolver::heat_step(Matrix2D& T, double dx, double dy, double dt)
{
    Matrix2D T_new = T;

    T_new.block(1, 1, m_ny - 2, m_nx - 2) =
        T.block(1, 1, m_ny - 2, m_nx - 2)
        + m_alpha * dt * (
            (T.block(1, 2, m_ny - 2, m_nx - 2)
             - 2.0 * T.block(1, 1, m_ny - 2, m_nx - 2)
             + T.block(1, 0, m_ny - 2, m_nx - 2)) / (dx * dx)
          + (T.block(2, 1, m_ny - 2, m_nx - 2)
             - 2.0 * T.block(1, 1, m_ny - 2, m_nx - 2)
             + T.block(0, 1, m_ny - 2, m_nx - 2)) / (dy * dy)
        );

    T_new.row(0).setZero();
    T_new.row(m_ny - 1).setZero();
    T_new.col(0).setZero();
    T_new.col(m_nx - 1).setZero();

    T = std::move(T_new);
}

void HeatSolver::compute_heat_flux(const Matrix2D& T,
                       double dx, double dy,
                       Matrix2D& qx, Matrix2D& qy)
{
    const Eigen::Index ny = static_cast<Eigen::Index>(m_ny);
    const Eigen::Index nx = static_cast<Eigen::Index>(m_nx);

    qx.setZero();
    qy.setZero();

    qx.block(0, 1, ny, nx - 2) =
        -(T.block(0, 2, ny, nx - 2) - T.block(0, 0, ny, nx - 2)) / (2.0 * dx);

    qy.block(1, 0, ny - 2, nx) =
        -(T.block(2, 0, ny - 2, nx) - T.block(0, 0, ny - 2, nx)) / (2.0 * dy);
}

bool HeatSolver::run_simulation(const std::string& out_dir)
{
    const double dx = m_lx / (m_nx - 1);
    const double dy = m_ly / (m_ny - 1);

    const double dt_max = (std::min(dx, dy) * std::min(dx, dy)) / (2.0 * m_alpha);
    const double dt     = m_cflSafety * dt_max;

    std::cout << "Heat equation simulation  (" << m_nx << "x" << m_ny << " grid)\n"
              << "  alpha=" << m_alpha
              << "  dx=" << dx
              << "  dt=" << dt << " (CFL limit " << dt_max << ", safety=" << m_cflSafety << ")\n"
              << "  " << m_nSteps << " steps, output every " << m_outputFrequency << " steps\n"
              << "  simulated delay=" << m_simulatedDelayMs << " ms per output step\n"
              << "  Writing to: " << out_dir << "/\n\n";

    // -----------------------------------------------------------------------
    // Initial conditions: three Gaussian hot-spots, cold Dirichlet boundaries
    // -----------------------------------------------------------------------
    Matrix2D T(m_ny, m_nx);

    for (int j = 0; j < m_ny; ++j)
    {
        const double y = j * dy;
        for (int i = 0; i < m_nx; ++i)
        {
            const double x = i * dx;
            T(j, i) =       std::exp(-((x-0.5) *(x-0.5)  + (y-0.5) *(y-0.5))  / 0.01)
                    + 0.7 * std::exp(-((x-0.25)*(x-0.25) + (y-0.75)*(y-0.75)) / 0.005)
                    + 0.5 * std::exp(-((x-0.75)*(x-0.75) + (y-0.25)*(y-0.25)) / 0.005);
        }
    }
    T.row(0).setZero();
    T.row(m_ny - 1).setZero();
    T.col(0).setZero();
    T.col(m_nx - 1).setZero();

    Matrix2D qx(m_ny, m_nx), qy(m_ny, m_nx);
    std::vector<std::pair<double, std::string>> pvd_entries;

    // -----------------------------------------------------------------------
    // Time loop
    // -----------------------------------------------------------------------
    for (int step = 0; step <= m_nSteps; ++step)
    {
        const double t = step * dt;

        if (step % m_outputFrequency == 0)
        {
            compute_heat_flux(T, dx, dy, qx, qy);

            std::ostringstream rel_name;
            rel_name << "heat_" << std::setw(4) << std::setfill('0') << step << ".vti";

            write_vti(out_dir, step, T, qx, qy, dx, dy);
            pvd_entries.emplace_back(t, rel_name.str());
            write_pvd(out_dir, pvd_entries);

            std::cout << "  step " << step << "/" << m_nSteps
                      << "  t=" << t
                      << "  T_max=" << T.maxCoeff()
                      << "  -> " << rel_name.str() << "\n";

            if (m_simulatedDelayMs > 0)
                std::this_thread::sleep_for(std::chrono::milliseconds(m_simulatedDelayMs));
        }

        if (step < m_nSteps)
            heat_step(T, dx, dy, dt);
    }

    std::cout << "\nSimulation complete.  Open " << out_dir << "/heat.pvd in ParaView.\n";

    return true;

}