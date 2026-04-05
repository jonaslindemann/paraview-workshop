#include "vtk_writer.h"

#include <string>
#include <vector>
#include <cstring>
#include <iostream>
#include <iomanip>
#include <sstream>

#include <pugixml.hpp>

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

static std::string encode_scalar(const double* data, size_t count)
{
    uint32_t nbytes = static_cast<uint32_t>(count * sizeof(double));

    std::vector<uint8_t> buf(sizeof(nbytes) + nbytes);
    std::memcpy(buf.data(), &nbytes, sizeof(nbytes));
    std::memcpy(buf.data() + sizeof(nbytes), data, nbytes);

    return base64_encode(buf.data(), buf.size());
}

static std::string encode_vector(const Matrix2D& qx, const Matrix2D& qy)
{
    const Eigen::Index rows = qx.rows();
    const Eigen::Index cols = qx.cols();

    if (rows != qy.rows() || cols != qy.cols()) {
        return {};
    }

    const size_t npts = static_cast<size_t>(rows * cols);
    std::vector<double> interleaved(npts * 3);

    for (Eigen::Index j = 0; j < rows; ++j)
        for (Eigen::Index i = 0; i < cols; ++i)
        {
            size_t base = static_cast<size_t>(j * cols + i) * 3;
            interleaved[base + 0] = qx(j, i);
            interleaved[base + 1] = qy(j, i);
            interleaved[base + 2] = 0.0;
        }

    return encode_scalar(interleaved.data(), interleaved.size());
}

void write_vti(const std::string& dir, int step,
               const Matrix2D& T,
               const Matrix2D& qx, const Matrix2D& qy,
               double dx, double dy)
{
    const Eigen::Index nx = T.cols();
    const Eigen::Index ny = T.rows();

    std::ostringstream fname;
    fname << dir << "/heat_" << std::setw(4) << std::setfill('0') << step << ".vti";

    // Build extent and spacing strings
    std::ostringstream extent, spacing;
    extent  << "0 " << (nx - 1) << " 0 " << (ny - 1) << " 0 0";
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
        encode_scalar(T.data(), static_cast<size_t>(nx * ny)).c_str());

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
