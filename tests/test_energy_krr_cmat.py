# MIT License
#
# Copyright (c) 2017 Anders Steen Christensen
#
# Permission is hereby granted, free of charge, to any person obtaining a copy
# of this software and associated documentation files (the "Software"), to deal
# in the Software without restriction, including without limitation the rights
# to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
# copies of the Software, and to permit persons to whom the Software is
# furnished to do so, subject to the following conditions:
#
# The above copyright notice and this permission notice shall be included in all
# copies or substantial portions of the Software.
#
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
# OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
# SOFTWARE.

from __future__ import print_function

import sys
import time
import numpy as np
import qml
from qml.kernels import laplacian_kernel
from qml.math import cho_solve


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

if __name__ == "__main__":

    # Parse file containing PBE0/def2-TZVP heats of formation and xyz filenames
    data = get_energies("hof_qm7.txt")

    # Generate a list of qml.Compound() objects
    mols = []

    for xyz_file in sorted(data.keys()):

        # Initialize the qml.Compound() objects
        mol = qml.Compound(xyz="qm7/" + xyz_file)

        # Associate a property (heat of formation) with the object
        mol.properties = data[xyz_file]

        # This is a Molecular Coulomb matrix sorted by row norm
        mol.generate_coulomb_matrix(size=23, sorting="row-norm")

        mols.append(mol)


    # Shuffle molecules
    np.random.seed(666)
    np.random.shuffle(mols)

    # Make training and test sets
    n_test  = 1000
    n_train = 4000

    training = mols[:n_train]
    test  = mols[-n_test:]

    # List of representations
    X  = np.array([mol.coulomb_matrix for mol in training])
    Xs = np.array([mol.coulomb_matrix for mol in test])

    # List of properties
    Y = np.array([mol.properties for mol in training])
    Ys = np.array([mol.properties for mol in test])

    # Set hyper-parameters
    sigma = 10**(4.2)
    llambda = 10**(-10.0)

    # Generate training Kernel
    print("Calculating training kernel ...")
    K = laplacian_kernel(X, X, sigma)

    # Solve alpha
    print("Solving alphas ...")
    K[np.diag_indices_from(K)] += llambda
    alpha = cho_solve(K,Y)

    # Calculate prediction kernel
    print("Calculating prediction kernel ...")
    Ks = laplacian_kernel(X, Xs, sigma)
    Yss = np.dot(Ks.transpose(), alpha)

    # Print final RMSD
    rmsd = np.sqrt(np.mean(np.square(Ys - Yss)))
    print("RMSD = %6.2f kcal/mol" % rmsd)