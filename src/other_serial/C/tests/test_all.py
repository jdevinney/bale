#/******************************************************************
#
#
#  Copyright(C) 2020, Institute for Defense Analyses
#  4850 Mark Center Drive, Alexandria, VA; 703-845-2500
# 
#
#  All rights reserved.
#  
#   This file is a part of Bale.  For license information see the
#   LICENSE file in the top level directory of the distribution.
#  
# 
# *****************************************************************/ 
import subprocess
import os

def test_all(path, implementation_mask):

  apps = ["histo", "ig", "topo", "randperm", "transpose_matrix", "permute_matrix", "triangles", "sssp", "unionfind"]
  for app in apps:

    assert(os.path.exists(os.path.join(path,app)))
    runs = []
    runs.append('--help ')
    if app == 'histo' or app == "ig":
      runs.append("-T 10 -N 13084 ")
      runs.append("-T 1 -N 12 ")
      runs.append("-T 5 -N 24242 ")
      runs.append("-N 1289 ")
    if app == 'topo':
      runs.append("-N 1082 -F -e 0.1 ")
      runs.append("-N 2341 -F -e 0.03")
      runs.append("-N 663 -F -e 0.3 ")
      runs.append("-N 1331 -G -e 0.08 ")
      runs.append("-N 133 -G -e 0.3 ")
      runs.append("-N 3332 -G -e 0.03 ")
    if app == 'transpose' or app == 'permute_matrix':
      runs.append("-N 1303 -e 0.1 ")
      runs.append("-N 233 -e 0.1 ")
      runs.append("-N 21 -e 0.3 ")
      runs.append("-N 2100 -e 0.03 ")
    if app == 'randperm':
      runs.append("-N 13093 ")
      runs.append("-N 133 ")
      runs.append("-N 2382 ")
      runs.append("-N 12 ")
    if app == 'triangles':
      runs.append("-N 3109 -e 0.05 ")
      runs.append("-N 601 -e 0.5 ")
      runs.append("-N 4409 -e 0.02 ")
      runs.append("-N 3109 -G -e 0.05 ")
      runs.append("-N 601 -G -e 0.5 ")
      runs.append("-N 4409 -G -e 0.02 ")
      runs.append('-K "0:3x4x5" ')
      runs.append('-K "1:3x4x5" ')
      runs.append('-K "2:3x4x5" ')
      runs.append('-K "0:13x17" ')
    if app == 'sssp':
      runs.append("-N 719 -F -e 0.05 ")
      runs.append("-N 601 -F -e 0.5 ")
      runs.append("-N 1213 -F -e 0.02 ")
      runs.append("-N 1109 -G -e 0.05 ")
      runs.append("-N 601 -G -e 0.5 ")
      runs.append("-N 1409 -G -e 0.02 ")
    if app == 'unionfind':
      runs.append("-N 719 -F -e 0.05 ")
      runs.append("-N 719 -G -e 0.05 ")
      runs.append("-N 601 -F -e 0.02 ")
      runs.append("-N 601 -G -e 0.02 ")
      runs.append("-N 1213 -F -e 0.01 ")
      runs.append("-N 1213 -G -e 0.01 ")
      
    for run in runs:
      cmd = "{0} {1} -M {2}".format(os.path.join(path, app), run, implementation_mask)
      print(cmd)
      cp = subprocess.run(cmd, shell=True)
      assert(cp.returncode == 0)
