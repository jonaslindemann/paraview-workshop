"""Read structured grid data and show
the associated vector and scalar fields"""
import pyvista as pv

g = pv.read("pipe.vtk")

# Slice through the middle of the pipe along Z axis
center = g.center
slc = g.slice(normal='z', origin=center)

# Build arrow glyphs on the slice
slc['vectors'] = slc.point_data['U']
slc.point_data.SetActiveVectors('vectors')
arrows = slc.glyph(orient='vectors', scale='vectors', factor=0.5, tolerance=0.05)

# Set up plotter
pl = pv.Plotter()
pl.add_mesh(g, scalars='p', cmap='jet', opacity=0.2)
pl.add_mesh(slc, scalars='p', cmap='jet')
pl.add_mesh(arrows, cmap='hot_r', scalar_bar_args={'title': 'Velocity (U)'})
pl.add_text(__doc__, position='upper_left', font_size=8)
pl.show_axes()
pl.view_isometric()
pl.show()
