#!/usr/bin/env python
import numpy as np
import sys


CGREEN = '\033[32m'
CRED = '\033[91m'
CEND = '\033[0m'

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

try:
    if np.allclose(*map(np.loadtxt, sys.argv[1:])):
        passed ('PASSED')
    else:
        error ('Failed value error')
except ValueError:
    error ("Failed, Shape error")

