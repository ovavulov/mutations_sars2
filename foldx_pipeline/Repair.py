import os
import sys
from scripts.functions import *
from datetime import datetime

init = datetime.now()
projectPath = os.getcwd()

# pipeline parameters reading
command = sys.argv
cfg = command[command.index("-f")+1]
params = getParams(cfg)

# шаг 1 - complex structure correction
makeRepair(**params)
os.chdir(projectPath)

fin = datetime.now()
print(f"\n\nTotal time spend: {fin - init}")




