import vtk
import numpy as np
import scipy.io as spi


if __name__ == "__main__":

    arrays = spi.loadmat("Data.mat")
    
    ex = arrays["ex"]
    ey = arrays["ey"]
    ez = arrays["ez"]
    u = arrays["u"]

    print(ex.shape)
    print(u.shape)
