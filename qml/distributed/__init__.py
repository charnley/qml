"""

Helper for generating config and input for distributed QML jobs

"""

import numpy as np
import os
import sys

import qml

import save_npy

def dump_fchl(jobname,
    A, B,
    two_body_scaling=np.sqrt(8),
    three_body_scaling=1.6,
    two_body_width=0.2,
    three_body_width=np.pi,
    two_body_power=4.0,
    three_body_power=2.0,
    cut_start=1.0,
    cut_distance=5.0,
    fourier_order=1,
    alchemy="periodic-table",
    alchemy_period_width=1.6,
    alchemy_group_width=1.6,
    kernel="gaussian", kernel_args=None,
    properties=None,
    test=False):

    """

    """

    name = "_fchl_"

    atoms_max = A.shape[1]
    neighbors_max = A.shape[3]

    assert B.shape[1] == atoms_max, "ERROR: Check FCHL representation sizes! code = 2"
    assert B.shape[3] == neighbors_max, "ERROR: Check FCHL representation sizes! code = 3"

    nm1 = A.shape[0]
    nm2 = B.shape[0]

    N1 = np.zeros((nm1),dtype=np.int32)
    N2 = np.zeros((nm2),dtype=np.int32)

    for a in range(nm1):
        N1[a] = len(np.where(A[a,:,1,0] > 0.0001)[0])

    for a in range(nm2):
        N2[a] = len(np.where(B[a,:,1,0] > 0.0001)[0])

    # Dump _N
    np.savetxt(jobname + name + "_a_n", N1, fmt='%i')
    np.savetxt(jobname + name + "_b_n", N2, fmt='%i')

    neighbors1 = np.zeros((nm1, atoms_max), dtype=np.int32)
    neighbors2 = np.zeros((nm2, atoms_max), dtype=np.int32)

    for a, representation in enumerate(A):
        ni = N1[a]
        for i, x in enumerate(representation[:ni]):
            neighbors1[a,i] = len(np.where(x[0]< cut_distance)[0])

    for a, representation in enumerate(B):
        ni = N2[a]
        for i, x in enumerate(representation[:ni]):
            neighbors2[a,i] = len(np.where(x[0]< cut_distance)[0])

    # Dump neighbors
    np.savetxt(jobname + name + "_a_neighbors", neighbors1, fmt='%i')
    np.savetxt(jobname + name + "_b_neighbors", neighbors2, fmt='%i')


    doalchemy, pd = qml.fchl.get_alchemy(alchemy, emax=100, r_width=alchemy_group_width, c_width=alchemy_period_width)
    doalchemy = int(doalchemy)

    # 
    kernel_idx, kernel_parameters, n_kernels = qml.fchl.get_kernel_parameters(kernel, kernel_args)

    kernel_parameters_shape = kernel_parameters.shape
    kernel_parameters_shape = np.array(kernel_parameters_shape)

    if len(kernel_parameters_shape) < 2:
        kernel_parameters_shape = list(kernel_parameters_shape) + [1]
        kernel_parameters_shape = np.array(kernel_parameters_shape)


    # Dump pd
    np.savetxt(jobname + name + "pd", pd)

    # other

    # _fchl_doalchemy
    np.savetxt(jobname + name + "doalchemy", [doalchemy], fmt="%i")

    # _fchl_cut_distance
    np.savetxt(jobname + name + "cut_distance", [cut_distance])

    # _fchl_cut_start
    np.savetxt(jobname + name + "cut_start", [cut_start])

    # _fchl_fourier_order
    np.savetxt(jobname + name + "fourier_order", [fourier_order], fmt="%i")

    # _fchl_shape
    # TODO what for?

    # _fchl_three_body_power
    np.savetxt(jobname + name + "three_body_power", [three_body_power])

    # _fchl_three_body_scaling
    np.savetxt(jobname + name + "three_body_scaling", [three_body_scaling])

    # _fchl_three_body_width
    np.savetxt(jobname + name + "three_body_width", [three_body_width])

    # _fchl_two_body_power
    np.savetxt(jobname + name + "two_body_power", [two_body_power])

    # _fchl_two_body_scaling
    np.savetxt(jobname + name + "two_body_scaling", [two_body_scaling])

    # _fchl_two_body_width
    np.savetxt(jobname + name + "two_body_width", [two_body_width])

    # _fchl_sigmas and _fchl_n_kernels
    np.savetxt(jobname + name + "n_kernels", [n_kernels])

    # _fchl_kernel_idx
    np.savetxt(jobname + name + "kernel_idx", [kernel_idx])

    # _fchl_kernel_parameters
    np.savetxt(jobname + name + "kernel_parameters", kernel_parameters)
    np.savetxt(jobname + name + "kernel_parameters_shape", kernel_parameters_shape, fmt="%i")

    if True:

        if properties is None:
            print("error: cannot test without properties")
            quit()

        verbose = True

        # gotta do krr
        fkernels = qml.fchl.fget_kernels_fchl(A, B, verbose, N1, N2, neighbors1,
                neighbors2, nm1, nm2, n_kernels, three_body_width,
                two_body_width, cut_start, cut_distance, fourier_order, pd,
                two_body_scaling, three_body_scaling, doalchemy,
                two_body_power, three_body_power,
                kernel_idx, kernel_parameters)

        fkernel = fkernels[0]
        del fkernels

        np.savetxt(jobname + name + "kernel", fkernel, fmt='%10.5f')

        properties = np.array(properties)

        alpha = qml.math.cho_solve(fkernel, properties)
        alpha = np.array([alpha])
        np.savetxt(jobname + name + "alpha", alpha, fmt='%10.5f')

    return


def make_cm(molecules, max_size=23, **kwargs):

    molecule.generate_coulomb_matrix(size=max_size, sorting="row-norm")

    return molecule.representation


def make_fchl(jobname, molecules, max_size=23, cut_distance=10.0, **kwargs):

    name = "_fchl_"

    representations = []
    atoms = []

    for mol in molecules:

        mol.generate_fchl_representation(max_size=max_size, cut_distance=cut_distance)
        natoms = mol.natoms
        representations.append(mol.representation)
        atoms.append(natoms)

    representations = np.array(representations)
    ni, nj, nk, nl = representations.shape

    # save representations
    save_npy.fsave_fchl_representations(jobname + name + "_a_representations", representations)
    save_npy.fsave_fchl_representations(jobname + name + "_b_representations", representations)

    # sizes
    np.savetxt(jobname + name + "_a_size", [ni], fmt="%i")
    np.savetxt(jobname + name + "_b_size", [ni], fmt="%i")
    np.savetxt(jobname + name + "max_size", [max_size], fmt="%i")
    np.savetxt(jobname + name + "max_neighbors", [nl], fmt="%i")

    # save kernel specific variables
    dump_fchl(jobname, representations, representations, **kwargs)

    return


def dump(jobname, molecules, properties=None, test=False, representation="fchl", **kwargs):


    if properties is not None:
        # save properties
        np.savetxt(jobname + "_properties", properties)

    if representation == "fchl":

        make_fchl(jobname, molecules, properties=properties, **kwargs)





    return

    print(representation)

    repvecs = []

    for molecule in molecules:
        print("convert")
        repvec = make_rep(molecule, **kwargs)
        repvecs.append(repvecs)


    repvecs = np.array(repvecs)

    np.savetxt(args.jobname + "/" + representation + "_representations", repvecs)

    return


def get_energies(filename):
    """ Returns a dictionary with heats of formation for each xyz-file.
    """

    f = open(filename, "r")
    lines = f.readlines()
    f.close()

    energies = dict()

    for line in lines:
        tokens = line.split()

        xyz_name = tokens[0]
        hof = float(tokens[1])

        energies[xyz_name] = hof

    return energies


def main():

    import argparse
    parser = argparse.ArgumentParser()

    parser.add_argument('--xyz', type=str, help="")
    parser.add_argument('--csv', type=str, help="csv of properties")

    parser.add_argument('-r', '--representation', type=str, help="", default="fchl")

    parser.add_argument('--jobname', type=str, help="folder name for the input files")

    parser.add_argument('--test', action="store_true", help="WARNING calculate and dump kernel and alphas as well")
    args = parser.parse_args()


    if os.path.exists(args.jobname):
        print("error: folder already exists", args.jobname)
        # quit()
    else:
        os.mkdir(args.jobname)


    molecules = []
    max_size = 0

    mol_properties = get_energies(args.csv)

    properties = []

    if os.path.isdir(args.xyz):
        # if args.xyz is folder, assume first col is filename

        if args.xyz[-1] != "/":
            args.xyz += "/"

        if args.jobname[-1] != "/":
            args.jobname += "/"

        files = os.listdir(args.xyz)
        files = list(files)
        files.sort()

        for molecule_path in files[:11]:
            molecule_path = molecule_path.strip()

            molecule_property = mol_properties[molecule_path]
            properties.append(molecule_property)

            molobj = qml.Compound(xyz=args.xyz + molecule_path)

            N = len(molobj.nuclear_charges)

            if max_size < N: max_size = N

            molecules.append(molobj)
    else:

        f = open(args.xyz, 'r')
        # TODO
        f.close()


    # save input
    dump(args.jobname, molecules, properties=properties, representation=args.representation)

    # representations = []
    #
    # for molecule in molecules:
    #     molecule.generate_coulomb_matrix(size=max_size, sorting="row-norm")
    #     representations.append(molecule.representation)

    return


if __name__ == "__main__":
    main()
