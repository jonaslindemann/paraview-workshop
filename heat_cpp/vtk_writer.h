#pragma once

#include <string>

#include <Eigen/Dense>

using Matrix2D = Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic,
                               Eigen::RowMajor>;

/**
 * Encode raw bytes to base64.
 *
 * VTK's "binary" inline format requires base64( uint32_header + raw_data ),
 * where the header is the byte count of the raw data as a little-endian UInt32.
 */                               
static std::string base64_encode(const uint8_t* data, size_t len);

/** Encode a contiguous array of doubles as VTK inline binary (base64). */
static std::string encode_scalar(const double* data, size_t count);

/**
 * Interleave qx, qy into a packed XYZ vector array (Z=0) and encode.
 * VTK requires 3-component vectors even for 2-D data.
 */
static std::string encode_vector(const Matrix2D& qx, const Matrix2D& qy);

/**
 * Write one VTK ImageData (.vti) file for the current timestep.
 *
 * Uses pugixml to build the XML tree and inline base64 binary encoding for
 * field data — compact, fast to parse, and no string concatenation.
 */
void write_vti(const std::string& dir, int step,
               const Matrix2D& T,
               const Matrix2D& qx, const Matrix2D& qy,
               double dx, double dy);

/**
 * Write (or overwrite) the ParaView Data collection file (.pvd).
 *
 * A .pvd file maps simulation time to .vti files so ParaView presents them
 * as an animated time series.  Rewritten after each output step so the
 * simulation can be loaded in ParaView while it is still running.
 */
 void write_pvd(const std::string& dir,
               const std::vector<std::pair<double, std::string>>& entries);