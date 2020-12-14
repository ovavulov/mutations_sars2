import os
import sys
from scripts.functions import *
from datetime import datetime
import numpy as np

ts1 = datetime.now()
projectPath = os.getcwd()

# читаем параметры пайплайна
command = sys.argv
cfg = command[command.index("-f")+1]
params = getParams(cfg)
assert params["repair"] == "done"
np.random.seed(int(params["rs"]))

# шаг 2 - сканируем замены остатков
makePositionScan(**params)
os.chdir(projectPath)

ts2 = datetime.now()
print(f"\n\nPosition Scan time spend: {ts2 - ts1}\n")
print("*****************************************************************")

# шаг 3 - сканируем замены остатков
buildModel(obj="rbd", **params)
os.chdir(projectPath)

ts3 = datetime.now()
print(f"\n\nBuild Model for RBD time spend: {ts3 - ts2}\n")
print("*****************************************************************")

# шаг 4 - сканируем замены остатков
buildModel(obj="complex", **params)
os.chdir(projectPath)

ts4 = datetime.now()
print(f"\n\nBuild Model for complex time spend: {ts4 - ts3}\n")
print("*****************************************************************")

# шаг 5 - сканируем замены остатков
analyseComplex(**params)
os.chdir(projectPath)

ts5 = datetime.now()
print(f"\n\nBuild Model for complex time spend: {ts5 - ts4}\n")
print("*****************************************************************")

print(f"\n\nTotal time spend: {ts5 - ts1}\n")




