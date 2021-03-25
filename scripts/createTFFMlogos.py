import sys
import TFFM.tffm_module as tffm_module
from TFFM.constants import TFFM_KIND

tffm_detailed = tffm_module.tffm_from_xml("tffm_detailed.xml",
        TFFM_KIND.DETAILED)

with open("/home/bgp01/tfbsshape/pictures/tffm_detailed_summary_logo.svg", "w") as pictureFile:
    tffm_detailed.print_summary_logo(pictureFile)
with open("/home/bgp01/tfbsshape/pictures/"+tffm_detailed_dense_logo.svg, "w") as pictureFile:
    tffm_detailed.print_dense_logo(pictureFile)
