# Catalyst v2 pipeline script for the heat equation simulation.
#
# This script is loaded by Catalyst at runtime and defines three hooks:
#   catalyst_initialize() -- called once at startup
#   catalyst_execute(info) -- called every time the simulation calls catalyst_execute()
#   catalyst_finalize()   -- called once at shutdown
#
# The simulation pushes data through a channel named "heat_grid".
# Catalyst makes that channel available as a TrivialProducer source whose
# registration name matches the channel name.
#
# To use this with a live ParaView session:
#   1. Launch ParaView
#   2. Catalyst > Connect  (default port 22222)
#   3. Run the simulation - ParaView shows the live data
#
# To use this in batch/offline mode (screenshot output only):
#   - Comment out the Live block in heat_equation.cpp (catalyst_live/enable)
#   - Screenshots are saved to catalyst_output/

from paraview.simple import *
from paraview import catalyst
import os

# ---------------------------------------------------------------------------
# Options object — required by the Catalyst v2 framework
# ---------------------------------------------------------------------------
options = catalyst.Options()
options.ExtractsOutputDirectory = "catalyst_output"
options.GenerateCinemaSpecification = 0

# ---------------------------------------------------------------------------
# Module-level state
# ---------------------------------------------------------------------------
_pipeline_initialized = False
_producer = None
_render_view = None
_display = None


# ---------------------------------------------------------------------------
# Catalyst hooks
# ---------------------------------------------------------------------------

def catalyst_initialize():
    """Called once when Catalyst starts up."""
    os.makedirs("catalyst_output", exist_ok=True)
    print("[Catalyst] Initialized. Output -> catalyst_output/")


def catalyst_execute(info):
    """
    Called each time the simulation invokes catalyst_execute().

    info.time      -- simulation time (float)
    info.timestep  -- simulation step index (int)
    """
    global _pipeline_initialized, _producer, _render_view, _display

    # ------------------------------------------------------------------
    # First call: build the ParaView pipeline once
    # ------------------------------------------------------------------
    if not _pipeline_initialized:
        _build_pipeline()
        _pipeline_initialized = True

    # ------------------------------------------------------------------
    # Update the source with the latest simulation data
    # ------------------------------------------------------------------
    _producer.UpdatePipeline(info.time)

    # ------------------------------------------------------------------
    # Render and save a screenshot
    # ------------------------------------------------------------------
    Render(_render_view)

    filename = os.path.join(
        "catalyst_output",
        f"heat_{info.timestep:04d}.png"
    )
    SaveScreenshot(filename, _render_view, ImageResolution=[800, 800])

    print(f"[Catalyst] step={info.timestep}  t={info.time:.5f}  -> {filename}")


def catalyst_finalize():
    """Called once when Catalyst shuts down."""
    print("[Catalyst] Finalized.")


# ---------------------------------------------------------------------------
# Pipeline construction (called once on first execute)
# ---------------------------------------------------------------------------

def _build_pipeline():
    global _producer, _render_view, _display

    # "heat_grid" is the channel name used in the C++ catalyst_push() call.
    # Catalyst registers it as a TrivialProducer with that name.
    _producer = GetProducerProxy("heat_grid")

    # --- Render view --------------------------------------------------------
    _render_view = CreateView("RenderView")
    _render_view.ViewSize = [800, 800]
    _render_view.Background = [0.15, 0.15, 0.15]

    # Set a top-down camera for the 2-D field
    _render_view.CameraPosition    = [0.5, 0.5, 3.0]
    _render_view.CameraFocalPoint  = [0.5, 0.5, 0.0]
    _render_view.CameraViewUp      = [0.0, 1.0, 0.0]
    _render_view.CameraParallelProjection = 1
    _render_view.CameraParallelScale = 0.55

    # --- Display: temperature scalar ----------------------------------------
    _display = Show(_producer, _render_view, "UniformGridRepresentation")
    _display.Representation = "Surface"

    ColorBy(_display, ("POINTS", "temperature"))
    _display.RescaleTransferFunctionToDataRange(True)

    # Use a blue-to-red color map
    lut = GetColorTransferFunction("temperature")
    lut.ApplyPreset("Cool to Warm", True)

    # Color legend
    scalar_bar = GetScalarBar(lut, _render_view)
    scalar_bar.Title = "Temperature"
    scalar_bar.ComponentTitle = ""
    scalar_bar.Visibility = 1
    _display.SetScalarBarVisibility(_render_view, True)

    # --- Glyph: heat flux vectors -------------------------------------------
    # Sub-sample with a stride so arrows don't clutter the view
    glyph = Glyph(
        registrationName="FluxGlyph",
        Input=_producer,
        GlyphType="Arrow"
    )
    glyph.OrientationArray  = ["POINTS", "heat_flux"]
    glyph.ScaleArray        = ["POINTS", "heat_flux"]
    glyph.ScaleFactor       = 0.05
    glyph.MaximumNumberOfSamplePoints = 400
    glyph.GlyphMode         = "Every Nth Point"
    glyph.Stride            = 10

    glyph_display = Show(glyph, _render_view, "GeometryRepresentation")
    glyph_display.Representation = "Surface"
    glyph_display.DiffuseColor = [1.0, 1.0, 0.0]   # yellow arrows

    ResetCamera(_render_view)
    print("[Catalyst] Pipeline built: temperature surface + heat-flux glyphs")
