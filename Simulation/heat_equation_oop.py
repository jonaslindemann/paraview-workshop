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

class HeatEquationSimulation:
    CFL_SAFETY = 0.4  # fraction of CFL limit; must be < 0.5 for stability

    def __init__(self, nx=80, ny=80, lx=1.0, ly=1.0, alpha=0.005):
        self.__nx = nx
        self.__ny = ny
        self.__lx = lx
        self.__ly = ly
        self.__alpha = alpha
        self.__dirty = True  # flag to indicate if grid/initial conditions need setup

        self.setup_grid()
        self.setup_initial_conditions()


    def setup_grid(self):
        self.dx = self.lx / (self.nx - 1)
        self.dy = self.ly / (self.ny - 1)
        x = np.linspace(0, self.lx, self.nx)
        y = np.linspace(0, self.ly, self.ny)
        self.X, self.Y = np.meshgrid(x, y)
        # Stable time step: dt = CFL_SAFETY * dx²/(2α)  [2D explicit CFL limit]
        self.__dt = self.CFL_SAFETY * min(self.dx, self.dy)**2 / (2.0 * self.alpha)

    def setup_initial_conditions(self):
        self.T = np.exp(-((self.X - 0.5)**2 + (self.Y - 0.5)**2) / 0.01)
        self.T += 0.7 * np.exp(-((self.X - 0.25)**2 + (self.Y - 0.75)**2) / 0.005)
        self.T += 0.5 * np.exp(-((self.X - 0.75)**2 + (self.Y - 0.25)**2) / 0.005)

        # Apply fixed cold boundaries
        self.T[ 0, :] = 0.0
        self.T[-1, :] = 0.0
        self.T[:,  0] = 0.0
        self.T[:, -1] = 0.0

    def step(self):
        T_new = self.T.copy()
        T_new[1:-1, 1:-1] = self.T[1:-1, 1:-1] + self.alpha * self.dt * (
            (self.T[1:-1, 2:]  - 2 * self.T[1:-1, 1:-1] + self.T[1:-1, :-2]) / self.dx**2 +
            (self.T[2:,   1:-1] - 2 * self.T[1:-1, 1:-1] + self.T[:-2,  1:-1]) / self.dy**2
        )
        # Re-apply fixed boundaries
        T_new[ 0, :] = 0.0
        T_new[-1, :] = 0.0
        T_new[:,  0] = 0.0
        T_new[:, -1] = 0.0
        self.T = T_new

    def write_vtk(self, path):
        dTdx = np.gradient(self.T, self.dx, axis=1)
        dTdy = np.gradient(self.T, self.dy, axis=0)

        def prep(arr):
            return np.asfortranarray(arr.T[:, :, np.newaxis].astype(np.float64))

        imageToVTK(
            path,
            origin=(0.0, 0.0, 0.0),
            spacing=(self.dx, self.dy, 1.0),
            pointData={
                "temperature": prep(self.T),
                "heat_flux": (prep(-dTdx), prep(-dTdy), prep(np.zeros_like(self.T))),
            },
        )

    def run(self, n_steps=200, write_every=4, output_dir="output", pause=0.05):
        os.makedirs(output_dir, exist_ok=True)
        file_index = 0

        if self.__dirty:
            self.setup_grid()
            self.setup_initial_conditions()
            self.__dirty = False

        dt_max = min(self.dx, self.dy)**2 / (2.0 * self.alpha)
        print(f"Heat equation simulation  ({self.nx}x{self.ny} grid)")
        print(f"  alpha={self.alpha}, dt={self.dt:.2e} (CFL limit {dt_max:.2e}, safety={self.CFL_SAFETY}), {n_steps} steps")
        print(f"  Writing every {write_every} steps -> {n_steps // write_every + 1} VTK files")
        print(f"  Output: {os.path.abspath(output_dir)}/\n")

        for step in range(n_steps + 1):
            t = step * self.dt

            if step % write_every == 0:
                filename = os.path.join(output_dir, f"heat_{file_index:04d}")
                self.write_vtk(filename)
                print(
                    f"  step {step:4d}/{n_steps}  "
                    f"t={t:.5f}  "
                    f"T_max={self.T.max():.4f}  "
                    f"-> {filename}"
                )
                file_index += 1
                time.sleep(pause)

            self.step()

        print(f"\nSimulation complete. {file_index} files in '{output_dir}/'.")
        print("\nOpen in ParaView:")
        print(f"  File > Open > {output_dir}/heat_0000.vti")
        print("  (ParaView will detect the full series automatically)")

    @property
    def nx(self):
        return self.__nx
    
    @nx.setter
    def nx(self, value):
        self.__nx = value
        self.__dirty = True
    
    @property
    def ny(self):
        return self.__ny
    
    @ny.setter
    def ny(self, value):
        self.__ny = value
        self.__dirty = True
    
    @property
    def lx(self):
        return self.__lx
    
    @lx.setter
    def lx(self, value):
        self.__lx = value
        self.__dirty = True
    
    @property
    def ly(self):
        return self.__ly    
    
    @ly.setter
    def ly(self, value):
        self.__ly = value
        self.__dirty = True
    
    @property
    def alpha(self):
        return self.__alpha 
    
    @alpha.setter
    def alpha(self, value):
        self.__alpha = value
        self.__dirty = True
    
    @property
    def dt(self):
        return self.__dt


if __name__ == "__main__":
    sim = HeatEquationSimulation()
    sim.nx = 500
    sim.ny = 500
    sim.run()