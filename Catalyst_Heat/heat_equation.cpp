/**
 * 2D Heat Equation Simulation with ParaView Catalyst In-Situ Visualization
 *
 * Solves:   dT/dt = alpha * (d2T/dx2 + d2T/dy2)
 *
 * Uses Eigen for the finite-difference computation and ParaView Catalyst v2
 * to push the mesh and field data to a live ParaView session (or to save
 * screenshots) without writing any files from the simulation itself.
 *
 * Catalyst API overview:
 *   catalyst_initialize()  -- start Catalyst, load the pipeline script
 *   catalyst_execute()     -- push data to ParaView at each timestep
 *   catalyst_finalize()    -- shut down Catalyst cleanly
 *
 * Data is described using the Conduit mesh blueprint, a self-describing
 * node-tree format.  The simulation fills a conduit node with:
 *   - mesh topology (uniform 2-D grid)
 *   - temperature scalar field (vertex-centered)
 *   - heat-flux vector field  (vertex-centered)
 * and hands it to catalyst_execute().
 *
 * Build:  see CMakeLists.txt
 * Run:    ./heat_equation  [catalyst_pipeline.py]
 */

#include <catalyst.hpp>         // Catalyst v2 C++ convenience wrapper
#include <conduit_cpp_to_c.hpp> // conduit_cpp::Node -> conduit_node* conversion

#include <Eigen/Dense>

#include <chrono>
#include <cstdlib>
#include <iostream>
#include <thread>

// ---------------------------------------------------------------------------
// Simulation parameters
// ---------------------------------------------------------------------------

static constexpr int    NX          = 200;    // grid points in x
static constexpr int    NY          = 200;    // grid points in y
static constexpr double LX          = 1.0;    // domain width  [m]
static constexpr double LY          = 1.0;    // domain height [m]
static constexpr double ALPHA       = 0.005;  // thermal diffusivity [m²/s]
static constexpr double CFL_SAFETY  = 0.4;    // fraction of stability limit (< 0.5)
static constexpr int    N_STEPS     = 400;    // time steps to run
static constexpr int    OUTPUT_FREQ = 4;      // call Catalyst every N steps
// Pause between Catalyst calls so a live ParaView can keep up (ms)
static constexpr int    PAUSE_MS    = 50;

// Row-major Eigen matrix: T(j, i) is the temperature at grid point (i, j)
// Row-major ensures T.data() is laid out with i (x) varying fastest,
// matching what VTK / conduit expects for a uniform structured grid.
using Matrix2D = Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic,
                               Eigen::RowMajor>;

// ---------------------------------------------------------------------------
// Helpers
// ---------------------------------------------------------------------------

/** Compute one explicit finite-difference step in-place. */
void heat_step(Matrix2D& T, double dx, double dy, double dt)
{
    // Take a copy of the interior so we do not read updated values mid-step
    Matrix2D T_new = T;

    T_new.block(1, 1, NY-2, NX-2) =
        T.block(1, 1, NY-2, NX-2)
        + ALPHA * dt * (
            // d2T/dx2  (column neighbours)
            (T.block(1, 2, NY-2, NX-2)
             - 2.0 * T.block(1, 1, NY-2, NX-2)
             + T.block(1, 0, NY-2, NX-2)) / (dx * dx)
            // d2T/dy2  (row neighbours)
          + (T.block(2, 1, NY-2, NX-2)
             - 2.0 * T.block(1, 1, NY-2, NX-2)
             + T.block(0, 1, NY-2, NX-2)) / (dy * dy)
        );

    // Re-apply Dirichlet (cold) boundary conditions
    T_new.row(0).setZero();
    T_new.row(NY-1).setZero();
    T_new.col(0).setZero();
    T_new.col(NX-1).setZero();

    T = std::move(T_new);
}

/**
 * Compute heat flux q = -grad(T) using central differences.
 * Stores results in pre-allocated row-major matrices.
 */
void compute_heat_flux(const Matrix2D& T,
                       double dx, double dy,
                       Matrix2D& qx, Matrix2D& qy)
{
    qx.setZero();
    qy.setZero();

    // Interior: central difference
    qx.block(0, 1, NY, NX-2) =
        -(T.block(0, 2, NY, NX-2) - T.block(0, 0, NY, NX-2)) / (2.0 * dx);

    qy.block(1, 0, NY-2, NX) =
        -(T.block(2, 0, NY-2, NX) - T.block(0, 0, NY-2, NX)) / (2.0 * dy);
}

// ---------------------------------------------------------------------------
// Catalyst interface
// ---------------------------------------------------------------------------

/**
 * Push the current simulation state to Catalyst.
 *
 * Builds a conduit node following the Conduit Mesh Blueprint for a 2-D
 * uniform (image) grid, then calls catalyst_execute().
 *
 * Node layout:
 *   catalyst/state/timestep          <int>
 *   catalyst/state/time              <double>
 *   catalyst/channels/<name>/type    "mesh"
 *   catalyst/channels/<name>/data/coordsets/coords/type    "uniform"
 *   catalyst/channels/<name>/data/coordsets/coords/dims/i  NX
 *   catalyst/channels/<name>/data/coordsets/coords/dims/j  NY
 *   catalyst/channels/<name>/data/coordsets/coords/origin/x  0.0
 *   catalyst/channels/<name>/data/coordsets/coords/origin/y  0.0
 *   catalyst/channels/<name>/data/coordsets/coords/spacing/dx dx
 *   catalyst/channels/<name>/data/coordsets/coords/spacing/dy dy
 *   catalyst/channels/<name>/data/topologies/mesh/type       "uniform"
 *   catalyst/channels/<name>/data/topologies/mesh/coordset   "coords"
 *   catalyst/channels/<name>/data/fields/temperature/...
 *   catalyst/channels/<name>/data/fields/heat_flux/...
 */
void catalyst_push(int step, double t,
                   const Matrix2D& T,
                   const Matrix2D& qx,
                   const Matrix2D& qy,
                   const Matrix2D& qz,
                   double dx, double dy)
{
    conduit_cpp::Node exec_params;

    // --- Simulation state ---------------------------------------------------
    exec_params["catalyst/state/timestep"].set(step);
    exec_params["catalyst/state/time"].set(t);

    // --- Channel: "heat_grid" ------------------------------------------------
    auto& channel = exec_params["catalyst/channels/heat_grid"];
    channel["type"].set("mesh");

    auto& mesh = channel["data"];

    // Coordinate system: uniform image grid
    mesh["coordsets/coords/type"].set("uniform");
    mesh["coordsets/coords/dims/i"].set(NX);
    mesh["coordsets/coords/dims/j"].set(NY);
    mesh["coordsets/coords/origin/x"].set(0.0);
    mesh["coordsets/coords/origin/y"].set(0.0);
    mesh["coordsets/coords/spacing/dx"].set(dx);
    mesh["coordsets/coords/spacing/dy"].set(dy);

    // Topology references the coordinate system
    mesh["topologies/mesh/type"].set("uniform");
    mesh["topologies/mesh/coordset"].set("coords");

    // --- Field: temperature (scalar, vertex-centered) -----------------------
    // set_external avoids a copy; T must remain valid until catalyst_execute returns
    mesh["fields/temperature/association"].set("vertex");
    mesh["fields/temperature/topology"].set("mesh");
    mesh["fields/temperature/volume_dependent"].set("false");
    mesh["fields/temperature/values"].set_external(
        const_cast<double*>(T.data()), static_cast<conduit_index_t>(NX * NY));

    // --- Field: heat_flux (vector, vertex-centered) -------------------------
    // Conduit stores vector components as separate named children under values
    mesh["fields/heat_flux/association"].set("vertex");
    mesh["fields/heat_flux/topology"].set("mesh");
    mesh["fields/heat_flux/volume_dependent"].set("false");
    mesh["fields/heat_flux/values/x"].set_external(
        const_cast<double*>(qx.data()), static_cast<conduit_index_t>(NX * NY));
    mesh["fields/heat_flux/values/y"].set_external(
        const_cast<double*>(qy.data()), static_cast<conduit_index_t>(NX * NY));
    mesh["fields/heat_flux/values/z"].set_external(
        const_cast<double*>(qz.data()), static_cast<conduit_index_t>(NX * NY));

    // --- Execute Catalyst ---------------------------------------------------
    catalyst_execute(conduit_cpp::c_node(&exec_params));
}

// ---------------------------------------------------------------------------
// Main
// ---------------------------------------------------------------------------

int main(int argc, char* argv[])
{
    const double dx = LX / (NX - 1);
    const double dy = LY / (NY - 1);

    // 2-D CFL stability limit for the explicit scheme:  dt <= dx²/(2α)
    const double dt_max = (std::min(dx, dy) * std::min(dx, dy)) / (2.0 * ALPHA);
    const double dt     = CFL_SAFETY * dt_max;

    std::cout << "Heat equation simulation  (" << NX << "x" << NY << " grid)\n"
              << "  alpha=" << ALPHA
              << "  dx=" << dx
              << "  dt=" << dt << " (CFL limit " << dt_max << ", safety=" << CFL_SAFETY << ")\n"
              << "  " << N_STEPS << " steps, Catalyst every " << OUTPUT_FREQ << " steps\n\n";

    // -----------------------------------------------------------------------
    // Catalyst initialisation
    // -----------------------------------------------------------------------
    {
        conduit_cpp::Node init_params;

        // Path to the pipeline script.  Use argv[1] if provided, else default.
        const char* script = (argc > 1) ? argv[1] : "catalyst_pipeline.py";
        init_params["catalyst/scripts/pipeline/filename"].set(script);

        // Enable Catalyst Live (lets ParaView connect interactively)
        // Comment this out if you only want offline screenshot output.
        init_params["catalyst/catalyst_live/enable"].set(1);
        init_params["catalyst/catalyst_live/port"].set(22222);

        catalyst_initialize(conduit_cpp::c_node(&init_params));
        std::cout << "Catalyst initialised with script: " << script << "\n\n";
    }

    // -----------------------------------------------------------------------
    // Initial conditions
    // -----------------------------------------------------------------------
    Matrix2D T(NY, NX);

    // Build coordinate grids
    for (int j = 0; j < NY; ++j)
    {
        const double y = j * dy;
        for (int i = 0; i < NX; ++i)
        {
            const double x = i * dx;
            // Three Gaussian hot-spots
            T(j, i) =       std::exp(-((x-0.5)*(x-0.5) + (y-0.5)*(y-0.5)) / 0.01)
                    + 0.7 * std::exp(-((x-0.25)*(x-0.25) + (y-0.75)*(y-0.75)) / 0.005)
                    + 0.5 * std::exp(-((x-0.75)*(x-0.75) + (y-0.25)*(y-0.25)) / 0.005);
        }
    }
    // Dirichlet cold boundaries
    T.row(0).setZero();
    T.row(NY-1).setZero();
    T.col(0).setZero();
    T.col(NX-1).setZero();

    // Pre-allocate flux arrays (reused each step)
    Matrix2D qx(NY, NX), qy(NY, NX), qz = Matrix2D::Zero(NY, NX);

    // -----------------------------------------------------------------------
    // Time loop
    // -----------------------------------------------------------------------
    for (int step = 0; step <= N_STEPS; ++step)
    {
        const double t = step * dt;

        if (step % OUTPUT_FREQ == 0)
        {
            compute_heat_flux(T, dx, dy, qx, qy);
            catalyst_push(step, t, T, qx, qy, qz, dx, dy);

            std::cout << "  step " << step << "/" << N_STEPS
                      << "  t=" << t
                      << "  T_max=" << T.maxCoeff() << "\n";

            // Small pause so a live ParaView client can process the update
            std::this_thread::sleep_for(std::chrono::milliseconds(PAUSE_MS));
        }

        if (step < N_STEPS)
            heat_step(T, dx, dy, dt);
    }

    // -----------------------------------------------------------------------
    // Catalyst finalisation
    // -----------------------------------------------------------------------
    {
        conduit_cpp::Node fin_params;
        catalyst_finalize(conduit_cpp::c_node(&fin_params));
    }

    std::cout << "\nSimulation complete.\n";
    return EXIT_SUCCESS;
}
