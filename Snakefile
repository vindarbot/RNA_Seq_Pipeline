import os
import sys
import re
import glob

CONDITIONS = os.listdir("Experience/")

SAMPLES = [ os.path.basename(x) for x in glob.glob("Experience/*/*") ] 

print(SAMPLES)

