import pyvtk as vtk
import numpy as np
import scipy.io as spi
import calfem.core as cfc


def convert_to_node_topo(edof, ex, ey, ez, n_dofs_per_node=3, ignore_first=True):
    """
    Routine to convert dof based topology and element coordinates to node based
    topology required for visualisation with VTK and other visualisation frameworks

    :param array edof: element topology [nel x (n_dofs_per_node)|(n_dofs_per_node+1)*n_nodes ]
    :param array ex: element x coordinates [nel x n_nodes]
    :param array ey: element y coordinates [nel x n_nodes]
    :param array ez: element z coordinates [nel x n_nodes]
    :param array n_dofs_per_node: number of dofs per node. (default = 3)
    :param boolean ignore_first: ignore first column of edof. (default = True)
    :return array coords: Array of node coordinates. [n_nodes x 3]
    :return array topo: Node topology. [nel x n_nodes]
    :return array node_dofs: Dofs for each node. [n_nodes x n_dofs_per_node]
    """

    node_hash_coords = {}
    node_hash_numbers = {}
    node_hash_dofs = {}
    el_hash_dofs = []

    nel, cols = edof.shape

    if ignore_first:
        tot_dofs = cols-1
    else:
        tot_dofs = cols

    n_nodes = int(tot_dofs / n_dofs_per_node)

    print("cols    =", tot_dofs)
    print("nel     =", nel)
    print("n_nodes =", n_nodes)

    for elx, ely, elz, dofs in zip(ex, ey, ez, edof):

        if ignore_first:
            el_dofs = dofs[1:]
        else:
            el_dofs = dofs

        # 0 1 2  3 4 5  6 7 8  9 12 11 

        el_dof = np.zeros((n_nodes, n_dofs_per_node), dtype=int)
        el_hash_topo = []

        for i in range(n_nodes):
            el_dof[i] = el_dofs[ (i*n_dofs_per_node):((i+1)*n_dofs_per_node) ]
            node_hash_coords[hash(tuple(el_dof[i]))] = [elx[i], ely[i], elz[i]]
            node_hash_numbers[hash(tuple(el_dof[i]))] = -1
            node_hash_dofs[hash(tuple(el_dof[i]))] = el_dof[i]
            el_hash_topo.append(hash(tuple(el_dof[i])))

        el_hash_dofs.append(el_hash_topo)

    coord_count = 0

    coords = []
    node_dofs = []

    for node_hash in node_hash_numbers.keys():
        node_hash_numbers[node_hash] = coord_count
        node_dofs.append(node_hash_dofs[node_hash])
        coord_count +=1

        coords.append(node_hash_coords[node_hash])

    topo = []

    for el_hashes in el_hash_dofs:
        topo.append([
            node_hash_numbers[el_hashes[0]], 
            node_hash_numbers[el_hashes[1]], 
            node_hash_numbers[el_hashes[2]], 
            node_hash_numbers[el_hashes[3]]
            ]
        )

    return coords, topo, node_dofs

if __name__ == "__main__":

    arrays = spi.loadmat("Data1.mat")
    print(arrays.keys())
    
    edof = arrays["Edof"]
    ex = arrays["ex"]
    ey = arrays["ey"]
    ez = arrays["ez"]
    u = arrays["u"]

    coords, topo, node_dofs = convert_to_node_topo(edof, ex, ey, ez)

    # --- Creating vector fields for VTK

    for c in range(1):

        point_data = vtk.PointData()

        vector_field = []

        u_real = u.astype(float)

        for dofs in node_dofs:
            vector_field.append(u_real[dofs-1,c])

        point_data.append(vtk.Vectors(vector_field, name="displacement"))

        structure = vtk.UnstructuredGrid(points=coords, tetra=topo)
        vtk_data = vtk.VtkData(structure, point_data)
        vtk_data.tofile("data%04d.vtk" % c, "ascii")

