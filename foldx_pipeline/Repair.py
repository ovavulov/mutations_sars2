import os
import sys
from scripts.functions import *
from datetime import datetime

init = datetime.now()
projectPath = os.getcwd()

# читаем параметры пайплайна
command = sys.argv
cfg = command[command.index("-f")+1]
params = getParams(cfg)

# шаг 1 - корректируем структуру комплекса
makeRepair(**params)
os.chdir(projectPath)

fin = datetime.now()
print(f"\n\nTotal time spend: {fin - init}")




