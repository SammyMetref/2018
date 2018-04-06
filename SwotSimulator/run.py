import sys
import os
import shutil

import logging
logger = logging.getLogger()
handler = logging.StreamHandler()
handler.setLevel(logging.INFO)
logger.addHandler(handler)
logger.setLevel(logging.INFO)


if  (len(sys.argv) < 2):
    file_param = os.getcwd() + os.sep + 'example' + os.sep + 'params_OSMOSIS_nonoise_20120901_20130901_1d.py'
    print("no params file specified, default is " + file_param)
else:
    file_param = str(sys.argv[1])
if os.path.isfile(file_param):
    #    basedir=os.path.dirname(swotsimulator.__file__)
    shutil.copyfile(file_param, 'params.py')
    # os.path.join(basedir,'params.py'))
else:
    print("Error: No such file: '%s'" % file_param)
    sys.exit()

import swotsimulator.run_simulator as run_simulator
run_simulator.run_simulator(file_param)
