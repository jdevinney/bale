def pytest_addoption(parser):
    parser.addoption("-P", "--path", action="store", default="./",
                     help="Specify the path to bale binaries")
    parser.addoption("-L", "--launcher_cmd", action="store", default="",
                     help="Specify the job launcher on your system (i.e. srun, oshrun, etc)")
    parser.addoption("--launcher_opts", action="store", default="",
                     help="Options to give to the launcher other than -n.")
    parser.addoption("--node_range", action="store", default="",
                     help="A range given with <start>,<end>,<stride> for the number of PEs to run on.")
    parser.addoption("-M", "--implementation_mask", action="store", default="31",
                     help="A bit mask of implementations to run. AGP = 1, "
                     "exstack=2, exstack2=4, conveyors=8,alternates=16")

def pytest_generate_tests(metafunc):
    option_value = metafunc.config.option.path
    if 'path' in metafunc.fixturenames and option_value is not None:
        metafunc.parametrize("path",[option_value])

    option_value = metafunc.config.option.launcher_cmd
    if 'launcher_cmd' in metafunc.fixturenames and option_value is not None:
        metafunc.parametrize("launcher_cmd",[option_value])
        
    option_value = metafunc.config.option.launcher_opts
    if 'launcher_opts' in metafunc.fixturenames and option_value is not None:
        metafunc.parametrize("launcher_opts",[option_value])

    option_value = metafunc.config.option.node_range
    if 'node_range' in metafunc.fixturenames and option_value is not None:
        metafunc.parametrize("node_range",[option_value])

    option_value = metafunc.config.option.implementation_mask
    if 'implementation_mask' in metafunc.fixturenames and option_value is not None:
        metafunc.parametrize("implementation_mask",[option_value])
