# This code is taken from the file `sxs_bbh_example.ipynb` in the git
# repository sxs-collarboration/catalog_tools/Examples (on 10/12/19)
# License
#
# MIT License

# Copyright (c) 2019 Simulating eXtreme Spacetimes

# Permission is hereby granted, free of charge, to any person obtaining a copy
# of this software and associated documentation files (the "Software"), to deal
# in the Software without restriction, including without limitation the rights
# to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
# copies of the Software, and to permit persons to whom the Software is
# furnished to do so, subject to the following conditions:

# The above copyright notice and this permission notice shall be included in all
# copies or substantial portions of the Software.

# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
# OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
# SOFTWARE.

# The following lines can be uncommented if any of the necessary libraries are
# not installed
# !pip install sxs
# !pip install h5py matplotlib scipy numpy tqdm
# !pip install tqdm

# Unlike the rest of python files in this project, which should perform in
# both python 2 and python 3, this executable may not work in python 2,
# On clusters which have python 2 installed, one option is to locally install
# Anaconda with python 3, and then use this to call this script

# Set which simulations will have their data downloaded

Simulations  = ['SXS:BBH:1975']

# For downloading data
import sxs
from sxs import zenodo as zen

# For interacting with the data
import h5py
import numpy as np
from matplotlib import pyplot as plt
import json

zen.download.matching("rhOverM_Asymptotic_GeometricUnits_CoM.h5", \
                      sxs_ids = Simulations, \
                      dry_run = False, \
                      highest_lev_only = False)
