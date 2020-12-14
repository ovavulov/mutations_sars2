import os
import sys
import glob
import pandas as pd

command = sys.argv
path = ""
n = 1
if "-from" in command:
    path = command[command.index("-from")+1]
if "-n" in command:
    n = int(command[command.index("-n")+1])

def doit(n):
    file_names = glob.glob(f'{path}/*.pdbqt')
    everything = []
    failures = []
    print('Found', len(file_names), 'pdbqt files')
    for file_name in file_names:
        file = open(file_name)
        lines = file.readlines()
        file.close()
        try:
            line = lines[1]
            result = float(line.split(':')[1].split()[0])
            everything.append([result, file_name])
        except:
            failures.append(file_name)
    everything = sorted(everything, key = lambda x: x[0])
    part = everything[:n]
    report = pd.DataFrame(index=["ligand", "affinity"])
    i = 0
    for aff, p in part:
        report[i] = [p, aff]
        i += 1
        if "-to" in command:
            to_path = command[command.index("-to")+1]
            if not os.path.exists(to_path):
                os.mkdir(to_path)
            log = ".".join(p.split(".")[:-1]+["log"])
            os.system(f"cp {p} {to_path}")
            os.system(f"cp {log} {to_path}")
    report = report.T
    print(report)
    if "-to" in command:
        report.to_csv(os.path.join(to_path, "report.csv"), index=False)
    print()
    if len(failures) > 0:
        print('WARNING:', len(failures), 'pdbqt files could not be processed')
    if "-to" in command:
        report = pd.DataFrame(index=["ligand", "affinity"])
        n = int(command[command.index("-n")+1])

doit(n)