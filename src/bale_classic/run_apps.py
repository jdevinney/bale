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
# script to run a suite of bale apps for experimentation or testing

import sys
if sys.version_info[0] < 3 or sys.version_info[1] < 7:
  print("This script requires at least Python version 3.7")
  sys.exit(1)

import subprocess
import os
import argparse
from shutil import which

def run_command(cmd):
  #ret = subprocess.run(cmd,shell=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
  ret = subprocess.run(cmd,shell=True, capture_output=True)
  return(ret)


def determine_launcher(launcher):
  if launcher is not None:
    return(""+launcher+" -n {0} {1} {2}")
  ret = run_command('srun --help')
  if ret.returncode == 0: return('srun -n {0} {1} {2}')
  ret = run_command('aprun --help')
  if ret.returncode== 0: return('aprun -n {0} {1} {2}')
  ret = run_command('oshrun --help')
  if ret.returncode== 0: return('oshrun -n {0} {1} {2}')
  ret = run_command('upcrun --help')
  if ret.returncode== 0: return('upcrun -n {0} {1} {2}')
  return("{2} -n {0} {1} --")

all_apps = []
all_apps.append("histo")
all_apps.append("ig")
all_apps.append("randperm")
all_apps.append("transpose_matrix")
all_apps.append("permute_matrix")
all_apps.append("topo")
all_apps.append("triangles")
all_apps.append("sssp")
all_apps.append("sparse_matrix_io")


def append_to_file(master_file, file_to_add, first):
    fin = open(file_to_add, 'r')
    newdata = fin.read()
    fin.close()
    fout = open(master_file, 'a')
    if not first:
        fout.write(",")
    fout.write(newdata)
    fout.close()        


def run_apps(app_dict, pes_list, cores_per_node, launcher_opts, option_str, impl_mask, json_file):
  
  # add user options
  base_cmd = "{0} ".format(option_str)
  
  # add implementation mask to the base_cmd
  if impl_mask is not None:
    base_cmd = "{0} -M {1}".format(base_cmd, impl_mask)
  if cores_per_node is not None:
    base_cmd = "{0} -c {1}".format(base_cmd, cores_per_node)
    
  # if the user wants the output to go to json, add that option
  # and write the initial opening brace in the json file
  if json_file is not None:
    # write initial [ to master json file
    fout = open(json_file,'w')
    fout.write("[\n")
    fout.close()
    tmp_json_file = "tmp_{0}".format(json_file)
    # add -j option to base_cmd
    base_cmd = "{0} -j {1}".format(base_cmd, tmp_json_file)
    
  for app in app_dict:
    first = True
    for pes in pes_list:
      cmd = launcher.format(pes, launcher_opts, app_dict[app]) +" "+base_cmd
      print(cmd)
      ret = run_command(cmd)
      assert(ret.returncode == 0)
      lines = ret.stderr.decode('utf-8')
      with open('run_all.log','a') as wp:
        wp.write(lines)
      if json_file is not None:
        append_to_file(json_file, tmp_json_file, first)
        first=False

  if json_file is not None:
    # write closing brace and close json file
    fout = open(json_file,"a")
    fout.write("]")
    fout.close()



############################################################
if __name__ == '__main__':
  
    
  parser = argparse.ArgumentParser(description="""Script to run any or all of bale apps.""")

  parser.add_argument('-a','--app',  action="append", dest='app_list',  help="Specify an app to run. "
                      " Adding multiple -a options will create a list of apps to run. "
                      " Apps must be in {0}".format(all_apps), default=[])
  parser.add_argument('-c','--cores_per_node', action="store", dest="cores_per_node",
                      help="Specify the number of cores per node.", default=None)
  parser.add_argument('-j','--json_file', action="store", dest='json_file',
                      help="Pass -j option to all apps. The final json file will have a record"
                      " for each run and will be written to filename which is the argument of this option",
                      default=None)
  parser.add_argument('-L','--launcher', action="store", dest="launcher",
                      help="Specify the launcher to use (for example: srun). "
                      "We do our best to automatically detect the launcher, but "
                      "there are cases where we will fail (for instance if you "
                      "have a UPC build but you have oshrun in your PATH).",default=None)
  parser.add_argument(   '--launcher_opts', action="store", dest='launcher_opts',
                         help="Pass these options onto the launcher, these could be srun options for example."
                         " The -n (num_tasks) option to launchers is handled separately"
                         " (using the --pes_list or --nodes_list options)."
                      " Do not specify the -n option for the launcher here.", default="")
  parser.add_argument('-M','--impl_mask', action="store", dest='impl_mask',
                      help="Pass \"-M IMPL_MASK\" to all apps", default=None)
  #parser.add_argument(     '--pe_range', action="store", dest='pe_range',
  #                         help="Specify the number of PEs to run on as a range <start>,<end>,<stride>. "
  #                         "The job will run each app with x PEs where x is in range(pe_range) PEs",
  #                         default=None)
  parser.add_argument(    '--pes_list', action="store", dest='pes_list',
                          help="Specify a comma separated list of PE sizes to run on.", default=None)
  parser.add_argument(    '--nodes_list', action="store", dest="nodes_list",
                          help="Specify a comma separated list of number of nodes to run on. "
                          "Overrides pes_list. Requires -c option to be specified.", default=None)
  parser.add_argument('-o',"--option_str", action="store", dest='option_str',
                      help="Specify a string to pass to all apps. Must be valid for all apps!", default="")
  parser.add_argument('-P','--path',      action="store", dest='path',
                      help="Specify path to binaries. If no path is specified "
                      " this script checks ./ and also the PATH environment variable.", default="")

  
  args = parser.parse_args()

  launcher = determine_launcher(args.launcher)
  
  app_dict = {}
  if len(args.app_list) == 0:
    for app in all_apps:
      app_dict[app] = None
  else:
    for app in args.app_list:
      if app not in all_apps:
        print("I don't know app {0}".format(app))
        exit(1)
      else:
        app_dict[app] = None


  # verify that all app binaries are in our given PATH
  if args.path != "":
    for app in app_dict:
      if not os.path.exists(os.path.join(args.path, app)):
        print("Can't find {0} in path {1}.".format(app, args.path))
        exit(1)
      else:
        app_dict[app] = os.path.join(args.path,app)

  else:
    # the user didn't supply a path, make sure the apps are in the PATH
    for app in app_dict:
      if os.path.exists(app):
        app_dict[app] = "./{0}".format(app)

      elif which(app) is not None:
        app_dict[app] = which(app)
        
      else:
        print("Can't find {0} in user $PATH.".format(app))
        exit(1)
  print()
  for app in app_dict:
    print("Using path {0} for {1}".format(app_dict[app], app))
  print()

  # figure out the sizes of the jobs to be run
  pes_list = []
  if args.nodes_list is not None:
    if args.cores_per_node is None:
      print("You must use the -c option with the --nodes_list option")
      exit(1)
    pes_list = [int(i)*int(args.cores_per_node) for i in args.nodes_list.split(',')]
  elif args.pes_list is not None:
    pes_list = [int(i) for i in args.pes_list.split(',')]
    if args.cores_per_node:
      print("WARNING: The cores_per_node option is not meant to be used with the pes_list option because "\
            "it assumes you are running on full nodes. If you are doing this, used the --nodes_list option.")
  else:
    print("You must specify a list of job sizes with --nodes_list or --pes_list")
    exit(1)

  # remove the run_all.log
  if os.path.exists('run_all.log'):
    os.remove('run_all.log')

  # call run_apps!
  run_apps(app_dict=app_dict,
           pes_list=pes_list,
           cores_per_node=args.cores_per_node,
           launcher_opts=args.launcher_opts,
           option_str=args.option_str,
           impl_mask=args.impl_mask,
           json_file=args.json_file)
