from pyevtk.hl import pointsToVTK 
import numpy as np  
npoints = 100  
x = np.random.rand(npoints)  
y = np.random.rand(npoints)  
z = np.random.rand(npoints)  
pressure = np.random.rand(npoints)  
temp = np.random.rand(npoints)  
pointsToVTK("./points", x, y, z, data = {"temp" : temp, "pressure" : pressure})