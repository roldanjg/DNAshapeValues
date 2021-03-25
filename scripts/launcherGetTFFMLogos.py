import os
import subprocess
import time
import sys
import TFFM.tffm_module as tffm_module
from TFFM.constants import TFFM_KIND

#Python 2.7

rootdir = '/home/bgp01/tfbsshape/data'

for subdir, dirs, files in os.walk(rootdir):
    for file in files:
        if 'tffm_detailed.xml' == file:
            if not '2000' in subdir:
                if not '201' in subdir:
                    if not 'old' in subdir:
                        tffm_detailed = tffm_module.tffm_from_xml(os.path.join(subdir, file), TFFM_KIND.DETAILED)
                        picturefile = os.path.join(subdir[27:].replace('/', '_'))

                        with open(
                                os.path.join(
                                    "/home/bgp01/tfbsshape/pictures", picturefile + "_tffm_detailed_summary_logo.svg"
                                ), "w"
                        ) as pictureFile:
                            tffm_detailed.print_summary_logo(pictureFile)

                        with open(
                                os.path.join(
                                    "/home/bgp01/tfbsshape/pictures", picturefile + "_tffm_detailed_dense_logo.svg"
                                ), "w"
                        ) as pictureFile:
                            tffm_detailed.print_dense_logo(pictureFile)
