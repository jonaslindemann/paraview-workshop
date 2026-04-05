#pragma once

#include <string>

#include <Eigen/Dense>

using Matrix2D = Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic,
                               Eigen::RowMajor>;

class HeatSolver {
private:
    int m_nx{200};                 // grid points in x
    int m_ny{200};                 // grid points in y
    double m_lx{1.0};              // domain width  [m]
    double m_ly{1.0};              // domain height [m]
    double m_alpha{0.005};         // thermal diffusivity [m²/s]
    double m_cflSafety{0.4};       // fraction of stability limit (< 0.5)
    int m_nSteps{400};             // time steps to run
    int m_outputFrequency{4};      // write a file every N steps
    int m_simulatedDelayMs{5000};   // artificial delay after each output step

public:
    HeatSolver(int gridSize, double alpha, double dt);
    void set_simulated_delay_ms(int delayMs);
    
    void heat_step(Matrix2D& T, double dx, double dy, double dt);
    void compute_heat_flux(const Matrix2D& T, double dx, double dy, Matrix2D& qx, Matrix2D& qy);    

    bool run_simulation(const std::string& out_dir);
};