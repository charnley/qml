#!/usr/bin/env python
import numpy as np
import sys

CGREEN = '\033[32m'
CRED = '\033[91m'
CEND = '\033[0m'

def load(txtfile):
    """ Read fortran scientific notation """
    converters = {0: lambda s: float(s.replace('D', 'E') or float('nan'))}
    arr = np.loadtxt(txtfile, converters=converters)
    return arr

def passed(msg):
    print(CGREEN + msg + CEND)
    return

def warning(msg):
    print(CRED + msg + CEND)
    return

def error(msg):
    print(CRED + msg + CEND)
    return

if len(sys.argv) != 3:
    print ('Usage: %s alphas1 alphas2' % sys.argv[0])
    exit(1)

args = sys.argv[1:]

try:

    txt1 = args[0]
    txt2 = args[0]

    arr1 = load(txt1)
    arr2 = load(txt2)

    if np.allclose(arr1, arr2):
        passed ('PASSED')

    else:
        error ('FAILED, Value error')

except ValueError:
    error ("FAILED, Shape error")

