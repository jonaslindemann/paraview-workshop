from trame.app import get_server
from trame.ui.vuetify import SinglePageLayout
from trame.widgets import vtk, vuetify

# Import all necessary VTK modules

from vtkmodules.vtkFiltersCore import vtkContourFilter  # Added for contour filter
from vtkmodules.vtkCommonDataModel import vtkStructuredPoints  # Added for structured points
from vtkmodules.vtkFiltersCore import vtkThreshold  # Added for trivial producer
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

class IsoVolumeVis:
    def __init__(self, csv_file):
        """Initialize the isosurface visualization with a CSV file."""
        
        self.csv_file = csv_file
        
        # Initialize VTK pipeline

        self.dataset = self.create_vtk_dataset()
        
        self.threshold = vtkThreshold()  
        self.threshold.SetInputData(self.dataset)
        self.threshold.SetUpperThreshold(50.0)
        self.threshold.SetLowerThreshold(20.0)
        self.threshold.SetThresholdFunction(vtkThreshold.THRESHOLD_BETWEEN)
        self.threshold.Update()
        
        self.mapper = vtkPolyDataMapper()  
        self.mapper.SetInputConnection(self.threshold.GetOutputPort())
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
    
    def update_threshold(self, lower_threshold, upper_threshold):
        """Update the threshold values."""
        
        self.threshold.SetUpperThreshold(upper_threshold)
        self.threshold.SetLowerThreshold(lower_threshold)
        self.threshold.Update()

    def update_threshold_lower(self, lower_threshold):
        """Update the lower threshold value."""
        
        self.threshold.SetLowerThreshold(lower_threshold)
        self.threshold.Update()

    def update_threshold_upper(self, upper_threshold):
        """Update the upper threshold value."""
        
        self.threshold.SetUpperThreshold(upper_threshold)
        self.threshold.Update()



# -----------------------------------------------------------------------------
# VTK pipeline
# -----------------------------------------------------------------------------

renderer = vtkRenderer()
renderWindow = vtkRenderWindow()
renderWindow.AddRenderer(renderer)
renderWindowInteractor = vtkRenderWindowInteractor()
renderWindowInteractor.SetRenderWindow(renderWindow)
renderWindowInteractor.GetInteractorStyle().SetCurrentStyleToTrackballCamera()

iso_volume_vis = IsoVolumeVis("scalar.csv")

renderer.AddActor(iso_volume_vis.actor)

# Add bounding box

cube_axes_actor = vtkCubeAxesActor()
cube_axes_actor.SetCamera(renderer.GetActiveCamera())
cube_axes_actor.SetBounds(iso_volume_vis.dataset.GetBounds())
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

@state.change("threshold_lower")
def update_threshold_lower(threshold_lower, **kwargs):
    iso_volume_vis.update_threshold_lower(threshold_lower)
    ctrl.view_update()

@state.change("threshold_upper")
def update_threshold(threshold_upper, **kwargs):
    iso_volume_vis.update_threshold_upper(threshold_upper)
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
            v_model=("threshold_lower", 20.0),
            min=10.0,
            max=160.0,
            step=1,
            hide_details=True,
            dense=True,
            style="max-width: 300px",
            label="Lower threshold",
        )
        vuetify.VSlider(
            v_model=("threshold_upper", 150.0),
            min=10.0,
            max=160.0,
            step=1,
            hide_details=True,
            dense=True,
            style="max-width: 300px",
            label="Upper threshold",
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