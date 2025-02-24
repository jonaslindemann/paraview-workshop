from trame.app import get_server
from trame.ui.vuetify import SinglePageLayout
from trame.widgets import vtk, vuetify

# Import all necessary VTK modules

from vtkmodules.vtkFiltersCore import vtkContourFilter  # Added for contour filter
from vtkmodules.vtkCommonDataModel import vtkStructuredPoints  # Added for structured points
from vtkmodules.vtkRenderingCore import (
    vtkActor,
    vtkPolyDataMapper,
    vtkRenderer,
    vtkRenderWindow,
    vtkRenderWindowInteractor,
    vtkProperty,  # Added for actor properties
)

from vtkmodules.vtkRenderingAnnotation import vtkCubeAxesActor  # For bounding box

from vtk.util import numpy_support

import numpy as np
import pandas as pd

# -----------------------------------------------------------------------------
# IsoSurface visualisation class:
# -----------------------------------------------------------------------------

class IsoSurfaceVis:
    def __init__(self, csv_file):
        """Initialize the isosurface visualization with a CSV file."""
        
        self.csv_file = csv_file
        
        # Initialize VTK pipeline

        self.dataset = self.create_vtk_dataset()
        
        self.contour = vtkContourFilter()  
        self.contour.SetInputData(self.dataset)
        self.contour.SetNumberOfContours(1);
        self.contour.SetValue(0, 50.0)
        
        self.mapper = vtkPolyDataMapper()  
        self.mapper.SetInputConnection(self.contour.GetOutputPort())
        self.actor = vtkActor()  
        self.actor.SetMapper(self.mapper)
        self.actor.GetProperty().SetColor(0.0, 1.0, 0.0)
        
    def create_vtk_dataset(self):
        """Create a VTK structured points dataset from a CSV file."""

        # Read CSV data using Pandas

        df = pd.read_csv(self.csv_file)
        
        # Get unique values to determine dimensions
        x_unique = np.sort(df['X32'].unique())
        y_unique = np.sort(df['Y32'].unique())
        z_unique = np.sort(df['Z32'].unique())
        
        # Create structured points dataset
        dataset = vtkStructuredPoints()  # Use imported class directly
        dataset.SetDimensions(len(x_unique), len(y_unique), len(z_unique))
        
        # Set the origin and spacing based on the data
        x_spacing = x_unique[1] - x_unique[0] if len(x_unique) > 1 else 1
        y_spacing = y_unique[1] - y_unique[0] if len(y_unique) > 1 else 1
        z_spacing = z_unique[1] - z_unique[0] if len(z_unique) > 1 else 1
        
        dataset.SetOrigin(x_unique[0], y_unique[0], z_unique[0])
        dataset.SetSpacing(x_spacing, y_spacing, z_spacing)
        
        # Create and set the scalar data
        scalar_array = df['S'].to_numpy()
        vtk_array = numpy_support.numpy_to_vtk(scalar_array)
        dataset.GetPointData().SetScalars(vtk_array)
        
        return dataset
        
    def update_iso_value(self, iso_value):
        """Update the isosurface value."""

        self.contour.SetNumberOfContours(1)
        self.contour.SetValue(0, iso_value)

# -----------------------------------------------------------------------------
# VTK pipeline
# -----------------------------------------------------------------------------

renderer = vtkRenderer()
renderWindow = vtkRenderWindow()
renderWindow.AddRenderer(renderer)
renderWindowInteractor = vtkRenderWindowInteractor()
renderWindowInteractor.SetRenderWindow(renderWindow)
renderWindowInteractor.GetInteractorStyle().SetCurrentStyleToTrackballCamera()

iso_surface_vis = IsoSurfaceVis("scalar.csv")

renderer.AddActor(iso_surface_vis.actor)

# Add bounding box

cube_axes_actor = vtkCubeAxesActor()
cube_axes_actor.SetCamera(renderer.GetActiveCamera())
cube_axes_actor.SetBounds(iso_surface_vis.dataset.GetBounds())
cube_axes_actor.SetXTitle("X")
cube_axes_actor.SetYTitle("Y")
cube_axes_actor.SetZTitle("Z")
cube_axes_actor.SetFlyModeToOuterEdges()

renderer.AddActor(cube_axes_actor)
renderer.ResetCamera()

# -----------------------------------------------------------------------------
# Trame setup
# -----------------------------------------------------------------------------

server = get_server(client_type="vue2")
state, ctrl = server.state, server.controller

# -----------------------------------------------------------------------------
# Functions
# -----------------------------------------------------------------------------

@state.change("iso_value")
def update_iso_value(iso_value, **kwargs):
    iso_surface_vis.update_iso_value(float(iso_value))
    ctrl.view_update()

def reset_camera():
    ctrl.view_reset_camera()

# -----------------------------------------------------------------------------
# GUI
# -----------------------------------------------------------------------------

with SinglePageLayout(server) as layout:
    layout.title.set_text("Isosurface Visualization")
    
    with layout.content:
        with vuetify.VContainer(
            fluid=True,
            classes="pa-0 fill-height",
        ):
            view = vtk.VtkLocalView(renderWindow)
            ctrl.view_update = view.update
            ctrl.view_reset_camera = view.reset_camera
    
    with layout.toolbar:
        vuetify.VSpacer()
        vuetify.VSlider(
            v_model=("iso_value", 50.0),
            min=10.0,
            max=160.0,
            step=1,
            hide_details=True,
            dense=True,
            style="max-width: 300px",
            label="Iso Value",
        )
        
        with vuetify.VBtn(icon=True, click=reset_camera):
            vuetify.VIcon("mdi-restore")
        
        vuetify.VDivider(vertical=True, classes="mx-2")
        
        vuetify.VSwitch(
            v_model="$vuetify.theme.dark",
            hide_details=True,
            dense=True,
        )
        
        with vuetify.VBtn(icon=True, click=ctrl.view_reset_camera):
            vuetify.VIcon("mdi-crop-free")

# -----------------------------------------------------------------------------
# Main
# -----------------------------------------------------------------------------

if __name__ == "__main__":
    server.start()