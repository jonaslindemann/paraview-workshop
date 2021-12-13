import pyvtk as vtk
import numpy as np
import scipy.io as spi
import calfem.core as cfc


if __name__ == "__main__":

    arrays = spi.loadmat("Data1.mat")
    print(arrays.keys())
    
    edof = arrays["Edof"]
    ex = arrays["ex"]
    ey = arrays["ey"]
    ez = arrays["ez"]
    u = arrays["u"]

    # --- Creating coords and topo arrays for VTK

    node_hash_coords = {}
    node_hash_numbers = {}
    node_hash_dofs = {}
    el_hash_dofs = []

    for elx, ely, elz, dofs in zip(ex, ey, ez, edof):
        el_dofs = dofs[1:]

        # 0 1 2  3 4 5  6 7 8  9 12 11 

        el_dof1 = el_dofs[:3]
        el_dof2 = el_dofs[3:6]
        el_dof3 = el_dofs[6:9]
        el_dof4 = el_dofs[9:]

        node_hash_coords[hash(tuple(el_dof1))] = [elx[0], ely[0], elz[0]]
        node_hash_coords[hash(tuple(el_dof2))] = [elx[1], ely[1], elz[1]]
        node_hash_coords[hash(tuple(el_dof3))] = [elx[2], ely[2], elz[2]]
        node_hash_coords[hash(tuple(el_dof4))] = [elx[3], ely[3], elz[3]]

        node_hash_numbers[hash(tuple(el_dof1))] = -1
        node_hash_numbers[hash(tuple(el_dof2))] = -1
        node_hash_numbers[hash(tuple(el_dof3))] = -1
        node_hash_numbers[hash(tuple(el_dof4))] = -1

        node_hash_dofs[hash(tuple(el_dof1))] = el_dof1
        node_hash_dofs[hash(tuple(el_dof2))] = el_dof2
        node_hash_dofs[hash(tuple(el_dof3))] = el_dof3
        node_hash_dofs[hash(tuple(el_dof4))] = el_dof4

        el_hash_dofs.append([hash(tuple(el_dof1)), hash(tuple(el_dof2)), hash(tuple(el_dof3)), hash(tuple(el_dof4))])

    coord_count = 0

    coords = []
    node_dofs = []

    for hash in node_hash_numbers.keys():
        node_hash_numbers[hash] = coord_count
        node_dofs.append(node_hash_dofs[hash])
        coord_count +=1

        coords.append(node_hash_coords[hash])

    topo = []

    for el_hashes in el_hash_dofs:
        topo.append([
            node_hash_numbers[el_hashes[0]], 
            node_hash_numbers[el_hashes[1]], 
            node_hash_numbers[el_hashes[2]], 
            node_hash_numbers[el_hashes[3]]
            ]
        )

    # --- Creating vector fields for VTK

    for c in range(200):

        point_data = vtk.PointData()

        vector_field = []

        u_real = u.astype(float)

        for dofs in node_dofs:
            vector_field.append(u_real[dofs-1,c])

        point_data.append(vtk.Vectors(vector_field, name="displacement"))

        structure = vtk.UnstructuredGrid(points=coords, tetra=topo)
        vtk_data = vtk.VtkData(structure, point_data)
        vtk_data.tofile("data%04d.vtk" % c, "ascii")

