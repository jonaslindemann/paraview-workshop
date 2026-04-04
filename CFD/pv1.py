# state file generated using paraview version 6.0.1
import paraview
paraview.compatibility.major = 6
paraview.compatibility.minor = 0

#### import the simple module from the paraview
from paraview.simple import *
#### disable automatic camera reset on 'Show'
paraview.simple._DisableFirstRenderCameraReset()

# ----------------------------------------------------------------
# setup views used in the visualization
# ----------------------------------------------------------------

# get the material library
materialLibrary1 = GetMaterialLibrary()

# Create a new 'Render View'
renderView1 = CreateView('RenderView')
renderView1.Set(
    ViewSize=[1875, 1218],
    CenterOfRotation=[2.78125, 2.78125, 8.28125],
    CameraPosition=[14.005785048018312, 13.195445168702804, 26.26675526996097],
    CameraFocalPoint=[-2.910355662967486, -2.3182858847157557, -0.712512660960262],
    CameraViewUp=[-0.22983066666566032, 0.8989846979299319, -0.37283290834883076],
    CameraFocalDisk=1.0,
    CameraParallelScale=9.167868055742295,
    OSPRayMaterialLibrary=materialLibrary1,
)

SetActiveView(None)

# ----------------------------------------------------------------
# setup view layouts
# ----------------------------------------------------------------

# create new layout object 'Layout #1'
layout1 = CreateLayout(name='Layout #1')
layout1.AssignView(0, renderView1)
layout1.SetSize(1875, 1218)

# ----------------------------------------------------------------
# restore active view
SetActiveView(renderView1)
# ----------------------------------------------------------------

# ----------------------------------------------------------------
# setup the data processing pipelines
# ----------------------------------------------------------------

# create a new 'Legacy VTK Reader'
pipevtk = LegacyVTKReader(registrationName='pipe.vtk', FileNames=['C:\\Users\\Jonas Lindemann\\Development\\paraview-workshop\\CFD\\pipe.vtk'])

# create a new 'Threshold'
threshold1 = Threshold(registrationName='Threshold1', Input=pipevtk)
threshold1.Set(
    Scalars=['CELLS', 'U'],
    SelectedComponent='Magnitude',
    LowerThreshold=0.4132666245981981,
    UpperThreshold=1.4250573262006831,
)

# create a new 'Clip'
clip1 = Clip(registrationName='Clip1', Input=pipevtk)

# init the 'Plane' selected for 'ClipType'
clip1.ClipType.Origin = [2.78125, 2.78125, 8.28125]

# init the 'Plane' selected for 'HyperTreeGridClipper'
clip1.HyperTreeGridClipper.Origin = [2.78125, 2.78125, 8.28125]

# create a new 'Glyph'
glyph1 = Glyph(registrationName='Glyph1', Input=pipevtk,
    GlyphType='Line')
glyph1.Set(
    OrientationArray=['CELLS', 'U'],
    ScaleArray=['CELLS', 'U'],
    ScaleFactor=0.43062500000000004,
)

# create a new 'Clip'
clip2 = Clip(registrationName='Clip2', Input=threshold1)

# init the 'Plane' selected for 'ClipType'
clip2.ClipType.Origin = [2.8125, 2.8125, 8.28125]

# init the 'Plane' selected for 'HyperTreeGridClipper'
clip2.HyperTreeGridClipper.Origin = [2.8125, 2.8125, 8.28125]

# ----------------------------------------------------------------
# setup the visualization in view 'renderView1'
# ----------------------------------------------------------------

# show data from clip2
clip2Display = Show(clip2, renderView1, 'UnstructuredGridRepresentation')

# get color transfer function/color map for 'U'
uLUT = GetColorTransferFunction('U')
uLUT.Set(
    RGBPoints=GenerateRGBPoints(
        range_min=0.0,
        range_max=1.4250573262006831,
    ),
    ScalarRangeInitialized=1.0,
)

# trace defaults for the display properties.
clip2Display.Set(
    Representation='Surface',
    ColorArrayName=['CELLS', 'U'],
    LookupTable=uLUT,
)

# init the 'Piecewise Function' selected for 'ScaleTransferFunction'
clip2Display.ScaleTransferFunction.Points = [-115102.234375, 0.0, 0.5, 0.0, 1049.041259765625, 1.0, 0.5, 0.0]

# init the 'Piecewise Function' selected for 'OpacityTransferFunction'
clip2Display.OpacityTransferFunction.Points = [-115102.234375, 0.0, 0.5, 0.0, 1049.041259765625, 1.0, 0.5, 0.0]

# show data from glyph1
glyph1Display = Show(glyph1, renderView1, 'GeometryRepresentation')

# trace defaults for the display properties.
glyph1Display.Set(
    Representation='Surface',
    ColorArrayName=['POINTS', ''],
)

# init the 'Piecewise Function' selected for 'ScaleTransferFunction'
glyph1Display.ScaleTransferFunction.Points = [-87912.46875, 0.0, 0.5, 0.0, 6489.4423828125, 1.0, 0.5, 0.0]

# init the 'Piecewise Function' selected for 'OpacityTransferFunction'
glyph1Display.OpacityTransferFunction.Points = [-87912.46875, 0.0, 0.5, 0.0, 6489.4423828125, 1.0, 0.5, 0.0]

# setup the color legend parameters for each legend in this view

# get color legend/bar for uLUT in view renderView1
uLUTColorBar = GetScalarBar(uLUT, renderView1)
uLUTColorBar.Set(
    Title='U',
    ComponentTitle='Magnitude',
)

# set color bar visibility
uLUTColorBar.Visibility = 1

# show color legend
clip2Display.SetScalarBarVisibility(renderView1, True)

# ----------------------------------------------------------------
# setup color maps and opacity maps used in the visualization
# note: the Get..() functions create a new object, if needed
# ----------------------------------------------------------------

# get opacity transfer function/opacity map for 'U'
uPWF = GetOpacityTransferFunction('U')
uPWF.Set(
    Points=[0.0, 0.0, 0.5, 0.0, 1.4250573262006831, 1.0, 0.5, 0.0],
    ScalarRangeInitialized=1,
)

# ----------------------------------------------------------------
# setup animation scene, tracks and keyframes
# note: the Get..() functions create a new object, if needed
# ----------------------------------------------------------------

# get the time-keeper
timeKeeper1 = GetTimeKeeper()

# initialize the timekeeper

# get time animation track
timeAnimationCue1 = GetTimeTrack()

# initialize the animation track

# get animation scene
animationScene1 = GetAnimationScene()

# initialize the animation scene
animationScene1.Set(
    ViewModules=renderView1,
    Cues=timeAnimationCue1,
    AnimationTime=0.0,
)

# initialize the animation scene

# ----------------------------------------------------------------
# restore active source
SetActiveSource(glyph1)
# ----------------------------------------------------------------


##--------------------------------------------
## You may need to add some code at the end of this python script depending on your usage, eg:
#
## Render all views to see them appears
# RenderAllViews()
#
## Interact with the view, usefull when running from pvpython
# Interact()
#
## Save a screenshot of the active view
# SaveScreenshot("path/to/screenshot.png")
#
## Save a screenshot of a layout (multiple splitted view)
# SaveScreenshot("path/to/screenshot.png", GetLayout())
#
## Save all "Extractors" from the pipeline browser
# SaveExtracts()
#
## Save a animation of the current active view
# SaveAnimation()
#
## Please refer to the documentation of paraview.simple
## https://www.paraview.org/paraview-docs/nightly/python/
##--------------------------------------------