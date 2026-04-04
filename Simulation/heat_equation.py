"""
2D Heat Equation Simulation — File-Based ParaView Visualization Example

Solves the 2D heat equation:
    dT/dt = alpha * (d2T/dx2 + d2T/dy2)

Writes binary VTK image files (.vti) via pyevtk at each timestep.
ParaView detects these as a file series and lets you animate through them.

How to visualize in ParaView:
    1. File > Open > select output/heat_0000.vti
       (ParaView auto-detects the series heat_..vti)
    2. Click Apply in the Properties panel
    3. Color by 'temperature'
    4. To see heat flux vectors: add a Glyph filter, set orientation to 'heat_flux'
    5. Use the Play button or time slider to animate

Dependencies: numpy, pyevtk  (pip install pyevtk)
"""

import numpy as np
import os
import time
from pyevtk.hl import imageToVTK

# ---------------------------------------------------------------------------
# Simulation parameters
# ---------------------------------------------------------------------------

NX, NY       = 80, 80      # grid resolution
LX, LY       = 1.0, 1.0   # physical domain size [m]
ALPHA        = 0.005       # thermal diffusivity [m²/s]
CFL_SAFETY   = 0.4         # fraction of CFL limit used for DT (< 0.5 for stability)
N_STEPS      = 200         # total number of time steps
WRITE_EVERY  = 4           # write a VTK file every N steps
OUTPUT_DIR   = "output"    # output directory (created if missing)
PAUSE        = 0.05        # seconds to sleep between writes (simulate real-time)

# ---------------------------------------------------------------------------
# Grid
# ---------------------------------------------------------------------------

dx = LX / (NX - 1)
dy = LY / (NY - 1)
x  = np.linspace(0, LX, NX)
y  = np.linspace(0, LY, NY)
X, Y = np.meshgrid(x, y)   # shape (NY, NX)

# Automatically choose a stable time step: dt = CFL_SAFETY * dx²/(2α)
# The explicit scheme requires dt <= dx²/(2α) in 2D (factor 2, not 4,
# because both directions contribute).
dt_max = min(dx, dy)**2 / (2.0 * ALPHA)
DT = CFL_SAFETY * dt_max

# ---------------------------------------------------------------------------
# Initial conditions
# ---------------------------------------------------------------------------

# Central hot spot
T = np.exp(-((X - 0.5)**2 + (Y - 0.5)**2) / 0.01)

# Two additional off-center sources
T += 0.7 * np.exp(-((X - 0.25)**2 + (Y - 0.75)**2) / 0.005)
T += 0.5 * np.exp(-((X - 0.75)**2 + (Y - 0.25)**2) / 0.005)

# Fixed cold boundaries (Dirichlet = 0)
T[ 0, :] = 0.0
T[-1, :] = 0.0
T[:,  0] = 0.0
T[:, -1] = 0.0

# ---------------------------------------------------------------------------
# VTK writer
# ---------------------------------------------------------------------------

def write_vtk(path: str, T: np.ndarray) -> None:
    """
    Write a binary VTK image file (.vti) via pyevtk containing:
      - temperature   (scalar, point data)
      - heat_flux     (vector, point data) — negative gradient of T

    'path' should be given WITHOUT extension; pyevtk appends .vti itself.
    """
    # Compute heat flux: q = -grad(T)
    dTdx = np.gradient(T, dx, axis=1)
    dTdy = np.gradient(T, dy, axis=0)

    # pyevtk requires 3D Fortran-contiguous float64 arrays shaped (nx, ny, nz)
    # Our data is (ny, nx) row-major, so transpose to (nx, ny) then add z=1.
    def prep(arr):
        return np.asfortranarray(arr.T[:, :, np.newaxis].astype(np.float64))

    imageToVTK(
        path,
        origin=(0.0, 0.0, 0.0),
        spacing=(dx, dy, 1.0),
        pointData={
            "temperature": prep(T),
            "heat_flux": (prep(-dTdx), prep(-dTdy), prep(np.zeros_like(T))),
        },
    )


# ---------------------------------------------------------------------------
# Time integration — explicit (forward Euler) finite differences
# ---------------------------------------------------------------------------

def heat_step(T: np.ndarray) -> np.ndarray:
    """One explicit finite-difference step of the 2D heat equation."""
    T_new = T.copy()
    T_new[1:-1, 1:-1] = T[1:-1, 1:-1] + ALPHA * DT * (
        (T[1:-1, 2:]  - 2 * T[1:-1, 1:-1] + T[1:-1, :-2]) / dx**2 +
        (T[2:,   1:-1] - 2 * T[1:-1, 1:-1] + T[:-2,  1:-1]) / dy**2
    )
    # Re-apply fixed boundaries
    T_new[ 0, :] = 0.0
    T_new[-1, :] = 0.0
    T_new[:,  0] = 0.0
    T_new[:, -1] = 0.0
    return T_new


# ---------------------------------------------------------------------------
# Main loop
# ---------------------------------------------------------------------------

os.makedirs(OUTPUT_DIR, exist_ok=True)

print(f"Heat equation simulation  ({NX}x{NY} grid)")
print(f"  alpha={ALPHA}, dt={DT:.2e} (CFL limit {dt_max:.2e}, safety={CFL_SAFETY}), {N_STEPS} steps")
print(f"  Writing every {WRITE_EVERY} steps -> {N_STEPS // WRITE_EVERY + 1} VTK files")
print(f"  Output: {os.path.abspath(OUTPUT_DIR)}/\n")

file_index = 0

for step in range(N_STEPS + 1):
    t = step * DT

    if step % WRITE_EVERY == 0:
        filename = os.path.join(OUTPUT_DIR, f"heat_{file_index:04d}")
        write_vtk(filename, T)
        print(
            f"  step {step:4d}/{N_STEPS}  "
            f"t={t:.5f}  "
            f"T_max={T.max():.4f}  "
            f"-> {filename}"
        )
        file_index += 1
        time.sleep(PAUSE)   # let ParaView pick up the file before the next one

    T = heat_step(T)

print(f"\nSimulation complete. {file_index} files in '{OUTPUT_DIR}/'.")
print("\nOpen in ParaView:")
print(f"  File > Open > {OUTPUT_DIR}/heat_0000.vti")
print("  (ParaView will detect the full series automatically)")
