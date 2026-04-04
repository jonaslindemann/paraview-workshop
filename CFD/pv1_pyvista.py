"""PyVista translation of pv1.py ParaView pipeline.

Pipeline:
  pipe.vtk
    ├── Threshold (U magnitude 0.413–1.425) → Clip → display colored by U magnitude
    └── Glyph (line glyphs oriented/scaled by U) → display
"""
import numpy as np
import pyvista as pv

# ----------------------------------------------------------------
# Load data
# ----------------------------------------------------------------
print("Loading data...")
mesh = pv.read("pipe.vtk")

# U is a cell-data vector field; compute its magnitude for thresholding/coloring
print("Computing U magnitude...")
U = mesh.cell_data["U"]
mesh.cell_data["U_mag"] = np.linalg.norm(U, axis=1)

# ----------------------------------------------------------------
# Threshold: keep cells where U magnitude is within [lower, upper]
# (mirrors ParaView Threshold1 with SelectedComponent='Magnitude')
# ----------------------------------------------------------------
print
lower, upper = 0.4132666245981981, 1.4250573262006831
thresholded = mesh.threshold([lower, upper], scalars="U_mag")

# ----------------------------------------------------------------
# Clip the thresholded mesh with a plane
# (mirrors ParaView Clip2 – default plane normal is X)
# ----------------------------------------------------------------
print("Clipping thresholded mesh...")
clip_origin = [2.8125, 2.8125, 8.28125]
clipped = thresholded.clip(normal="x", origin=clip_origin)

# ----------------------------------------------------------------
# Glyphs on the original mesh
# (mirrors ParaView Glyph1 – Line glyphs oriented/scaled by cell U)
# Cell data must be promoted to point data first for glyph orientation.
# ----------------------------------------------------------------
print("Generating glyphs...")
mesh_pt = mesh.cell_data_to_point_data()
glyphs = mesh_pt.glyph(
    orient="U",
    scale="U",
    factor=0.43062500000000004,
    geom=pv.Line((-0.5, 0, 0), (0.5, 0, 0)),  # unit line, like ParaView 'Line' glyph type
    tolerance=0.05,                             # subsample to keep glyph count manageable
)

# ----------------------------------------------------------------
# Render
# ----------------------------------------------------------------

print("Rendering...")
pl = pv.Plotter()

# Clipped + thresholded volume colored by U magnitude
print("Adding clipped mesh to plotter...")
pl.add_mesh(
    clipped,
    scalars="U_mag",
    clim=[0.0, upper],
    cmap="coolwarm",
    scalar_bar_args={"title": "U Magnitude"},
)

# Velocity glyphs

print("Adding glyphs to plotter...")    
pl.add_mesh(glyphs, color="white")

# Restore the camera from the original ParaView state
pl.camera_position = [
    (14.005785048018312, 13.195445168702804, 26.26675526996097),   # position
    (-2.910355662967486, -2.3182858847157557, -0.712512660960262), # focal point
    (-0.22983066666566032, 0.8989846979299319, -0.37283290834883076), # view up
]

pl.show_axes()
pl.show()
