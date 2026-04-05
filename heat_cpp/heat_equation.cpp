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

#include <pugixml.hpp>
#include <Eigen/Dense>

#include <cstdint>
#include <cstdlib>
#include <cstring>
#include <iostream>
#include <iomanip>
#include <sstream>
#include <string>
#include <vector>

#ifdef _WIN32
#include <direct.h>
#define MKDIR(p) _mkdir(p)
#else
#include <sys/stat.h>
#define MKDIR(p) mkdir(p, 0755)
#endif

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
static constexpr int    OUTPUT_FREQ = 4;      // write a file every N steps

// Row-major Eigen matrix: T(j, i) is the temperature at grid point (i, j).
// Row-major ensures T.data() is laid out with i (x) varying fastest,
// which matches VTK ImageData's expected memory order.
using Matrix2D = Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic,
                               Eigen::RowMajor>;

// ---------------------------------------------------------------------------
// Simulation helpers
// ---------------------------------------------------------------------------

/** Compute one explicit finite-difference step in-place. */
void heat_step(Matrix2D& T, double dx, double dy, double dt)
{
    Matrix2D T_new = T;

    T_new.block(1, 1, NY-2, NX-2) =
        T.block(1, 1, NY-2, NX-2)
        + ALPHA * dt * (
            (T.block(1, 2, NY-2, NX-2)
             - 2.0 * T.block(1, 1, NY-2, NX-2)
             + T.block(1, 0, NY-2, NX-2)) / (dx * dx)
          + (T.block(2, 1, NY-2, NX-2)
             - 2.0 * T.block(1, 1, NY-2, NX-2)
             + T.block(0, 1, NY-2, NX-2)) / (dy * dy)
        );

    T_new.row(0).setZero();
    T_new.row(NY-1).setZero();
    T_new.col(0).setZero();
    T_new.col(NX-1).setZero();

    T = std::move(T_new);
}

/** Compute heat flux q = -grad(T) using central differences. */
void compute_heat_flux(const Matrix2D& T,
                       double dx, double dy,
                       Matrix2D& qx, Matrix2D& qy)
{
    qx.setZero();
    qy.setZero();

    qx.block(0, 1, NY, NX-2) =
        -(T.block(0, 2, NY, NX-2) - T.block(0, 0, NY, NX-2)) / (2.0 * dx);

    qy.block(1, 0, NY-2, NX) =
        -(T.block(2, 0, NY-2, NX) - T.block(0, 0, NY-2, NX)) / (2.0 * dy);
}

// ---------------------------------------------------------------------------
// Base64 encoding
// ---------------------------------------------------------------------------

/**
 * Encode raw bytes to base64.
 *
 * VTK's "binary" inline format requires base64( uint32_header + raw_data ),
 * where the header is the byte count of the raw data as a little-endian UInt32.
 */
static std::string base64_encode(const uint8_t* data, size_t len)
{
    static constexpr char table[] =
        "ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz0123456789+/";

    std::string out;
    out.reserve(((len + 2) / 3) * 4);

    for (size_t i = 0; i < len; i += 3)
    {
        uint8_t b0 =           data[i];
        uint8_t b1 = i+1<len ? data[i+1] : 0;
        uint8_t b2 = i+2<len ? data[i+2] : 0;

        out += table[b0 >> 2];
        out += table[((b0 & 0x03) << 4) | (b1 >> 4)];
        out += i+1<len ? table[((b1 & 0x0f) << 2) | (b2 >> 6)] : '=';
        out += i+2<len ? table[b2 & 0x3f]                       : '=';
    }
    return out;
}

/** Encode a contiguous array of doubles as VTK inline binary (base64). */
static std::string encode_scalar(const double* data, size_t count)
{
    uint32_t nbytes = static_cast<uint32_t>(count * sizeof(double));

    std::vector<uint8_t> buf(sizeof(nbytes) + nbytes);
    std::memcpy(buf.data(), &nbytes, sizeof(nbytes));
    std::memcpy(buf.data() + sizeof(nbytes), data, nbytes);

    return base64_encode(buf.data(), buf.size());
}

/**
 * Interleave qx, qy into a packed XYZ vector array (Z=0) and encode.
 * VTK requires 3-component vectors even for 2-D data.
 */
static std::string encode_vector(const Matrix2D& qx, const Matrix2D& qy)
{
    const size_t npts = static_cast<size_t>(NX * NY);
    std::vector<double> interleaved(npts * 3);

    for (int j = 0; j < NY; ++j)
        for (int i = 0; i < NX; ++i)
        {
            size_t base = static_cast<size_t>(j * NX + i) * 3;
            interleaved[base + 0] = qx(j, i);
            interleaved[base + 1] = qy(j, i);
            interleaved[base + 2] = 0.0;
        }

    return encode_scalar(interleaved.data(), interleaved.size());
}

// ---------------------------------------------------------------------------
// VTK output
// ---------------------------------------------------------------------------

/**
 * Write one VTK ImageData (.vti) file for the current timestep.
 *
 * Uses pugixml to build the XML tree and inline base64 binary encoding for
 * field data — compact, fast to parse, and no string concatenation.
 */
void write_vti(const std::string& dir, int step,
               const Matrix2D& T,
               const Matrix2D& qx, const Matrix2D& qy,
               double dx, double dy)
{
    std::ostringstream fname;
    fname << dir << "/heat_" << std::setw(4) << std::setfill('0') << step << ".vti";

    // Build extent and spacing strings
    std::ostringstream extent, spacing;
    extent  << "0 " << NX-1 << " 0 " << NY-1 << " 0 0";
    spacing << dx << " " << dy << " 1";

    pugi::xml_document doc;
    doc.append_child(pugi::node_declaration).append_attribute("version") = "1.0";

    auto vtk = doc.append_child("VTKFile");
    vtk.append_attribute("type")        = "ImageData";
    vtk.append_attribute("version")     = "0.1";
    vtk.append_attribute("byte_order")  = "LittleEndian";
    vtk.append_attribute("header_type") = "UInt32";

    auto image = vtk.append_child("ImageData");
    image.append_attribute("WholeExtent") = extent.str().c_str();
    image.append_attribute("Origin")      = "0 0 0";
    image.append_attribute("Spacing")     = spacing.str().c_str();

    auto piece = image.append_child("Piece");
    piece.append_attribute("Extent") = extent.str().c_str();

    auto point_data = piece.append_child("PointData");

    // Temperature (scalar)
    auto temp = point_data.append_child("DataArray");
    temp.append_attribute("type")               = "Float64";
    temp.append_attribute("Name")               = "temperature";
    temp.append_attribute("NumberOfComponents") = "1";
    temp.append_attribute("format")             = "binary";
    temp.append_child(pugi::node_pcdata).set_value(
        encode_scalar(T.data(), NX * NY).c_str());

    // Heat flux (3-component vector)
    auto flux = point_data.append_child("DataArray");
    flux.append_attribute("type")               = "Float64";
    flux.append_attribute("Name")               = "heat_flux";
    flux.append_attribute("NumberOfComponents") = "3";
    flux.append_attribute("format")             = "binary";
    flux.append_child(pugi::node_pcdata).set_value(
        encode_vector(qx, qy).c_str());

    if (!doc.save_file(fname.str().c_str(), "  "))
        std::cerr << "ERROR: could not write " << fname.str() << "\n";
}

/**
 * Write (or overwrite) the ParaView Data collection file (.pvd).
 *
 * A .pvd file maps simulation time to .vti files so ParaView presents them
 * as an animated time series.  Rewritten after each output step so the
 * simulation can be loaded in ParaView while it is still running.
 */
void write_pvd(const std::string& dir,
               const std::vector<std::pair<double, std::string>>& entries)
{
    pugi::xml_document doc;
    doc.append_child(pugi::node_declaration).append_attribute("version") = "1.0";

    auto vtk = doc.append_child("VTKFile");
    vtk.append_attribute("type")       = "Collection";
    vtk.append_attribute("version")    = "0.1";
    vtk.append_attribute("byte_order") = "LittleEndian";

    auto collection = vtk.append_child("Collection");

    for (auto& [t, filename] : entries)
    {
        auto dataset = collection.append_child("DataSet");
        dataset.append_attribute("timestep") = t;
        dataset.append_attribute("group")    = "";
        dataset.append_attribute("part")     = "0";
        dataset.append_attribute("file")     = filename.c_str();
    }

    std::string pvd_path = dir + "/heat.pvd";
    if (!doc.save_file(pvd_path.c_str(), "  "))
        std::cerr << "ERROR: could not write " << pvd_path << "\n";
}

// ---------------------------------------------------------------------------
// Main
// ---------------------------------------------------------------------------

int main(int argc, char* argv[])
{
    const std::string out_dir = (argc > 1) ? argv[1] : "heat_output";

    if (MKDIR(out_dir.c_str()) != 0 && errno != EEXIST)
    {
        std::cerr << "ERROR: could not create directory " << out_dir << "\n";
        return EXIT_FAILURE;
    }

    const double dx = LX / (NX - 1);
    const double dy = LY / (NY - 1);

    const double dt_max = (std::min(dx, dy) * std::min(dx, dy)) / (2.0 * ALPHA);
    const double dt     = CFL_SAFETY * dt_max;

    std::cout << "Heat equation simulation  (" << NX << "x" << NY << " grid)\n"
              << "  alpha=" << ALPHA
              << "  dx=" << dx
              << "  dt=" << dt << " (CFL limit " << dt_max << ", safety=" << CFL_SAFETY << ")\n"
              << "  " << N_STEPS << " steps, output every " << OUTPUT_FREQ << " steps\n"
              << "  Writing to: " << out_dir << "/\n\n";

    // -----------------------------------------------------------------------
    // Initial conditions: three Gaussian hot-spots, cold Dirichlet boundaries
    // -----------------------------------------------------------------------
    Matrix2D T(NY, NX);

    for (int j = 0; j < NY; ++j)
    {
        const double y = j * dy;
        for (int i = 0; i < NX; ++i)
        {
            const double x = i * dx;
            T(j, i) =       std::exp(-((x-0.5) *(x-0.5)  + (y-0.5) *(y-0.5))  / 0.01)
                    + 0.7 * std::exp(-((x-0.25)*(x-0.25) + (y-0.75)*(y-0.75)) / 0.005)
                    + 0.5 * std::exp(-((x-0.75)*(x-0.75) + (y-0.25)*(y-0.25)) / 0.005);
        }
    }
    T.row(0).setZero();
    T.row(NY-1).setZero();
    T.col(0).setZero();
    T.col(NX-1).setZero();

    Matrix2D qx(NY, NX), qy(NY, NX);
    std::vector<std::pair<double, std::string>> pvd_entries;

    // -----------------------------------------------------------------------
    // Time loop
    // -----------------------------------------------------------------------
    for (int step = 0; step <= N_STEPS; ++step)
    {
        const double t = step * dt;

        if (step % OUTPUT_FREQ == 0)
        {
            compute_heat_flux(T, dx, dy, qx, qy);

            std::ostringstream rel_name;
            rel_name << "heat_" << std::setw(4) << std::setfill('0') << step << ".vti";

            write_vti(out_dir, step, T, qx, qy, dx, dy);
            pvd_entries.emplace_back(t, rel_name.str());
            write_pvd(out_dir, pvd_entries);

            std::cout << "  step " << step << "/" << N_STEPS
                      << "  t=" << t
                      << "  T_max=" << T.maxCoeff()
                      << "  -> " << rel_name.str() << "\n";
        }

        if (step < N_STEPS)
            heat_step(T, dx, dy, dt);
    }

    std::cout << "\nSimulation complete.  Open " << out_dir << "/heat.pvd in ParaView.\n";
    return EXIT_SUCCESS;
}
