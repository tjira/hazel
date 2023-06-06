from src.integral import *
import os, psi4

# get the plugin directory
plugdir = os.path.split(os.path.abspath(__file__))[0]

# load the plugin
psi4.core.plugin_load(plugdir + '/' + os.path.split(plugdir)[1] + ".so")

# close the output
psi4.core.close_outfile()

# clear the output
open(psi4.core.get_output_file(), "w").close()

# reopen the output
psi4.core.reopen_outfile()
