import os, sys

command = sys.argv
suffix = ""
nthreads = 1
if "-db" in command:
    suffix = "_" + command[command.index("-db")+1]
if "-threads" in command:
    nthreads = int(command[command.index("-threads")+1])

fs = os.listdir("fda_library" + suffix)
zincs = []
for f_ in fs:
    if f_.startswith("ZINC000") and f_.endswith(".pdbqt"):
        f_ = f_.replace(".pdbqt","")
        zincs.append(f_)
fns = "screening" + suffix

if not os.path.exists(fns):
    os.mkdir(fns)

f = open("config_template.txt")
config = f.readlines()
config = "".join(config)
config = config.replace("FOLDER","fda_library"+suffix)

fo_run = open(f"start_screening{suffix}.sh",'w')
fo_run.write("date\n")
fo_run.write("cd "+fns+"\n")
for i in range(len(zincs)):
    zinc_id = zincs[i]
    config_ = config.replace("LIGAND",zinc_id)
    fo = open(fns+"/config_"+zinc_id+".txt",'w')
    fo.write(config_)
    fo_run.write(f"../vina --config config_{zinc_id}.txt &>/dev/null &\n")
    if (i + 1) % nthreads == 0 or i == len(zincs) - 1:
        fo_run.write("wait\n")
        proc = int((i+1)/len(zincs)*100)
        fo_run.write(f"""printf "{"#"*proc}{"_"*(100-proc)}({proc}%%)\\r"\n""")
        
fo_run.write("date\n")
os.system(f"chmod 777 start_screening{suffix}.sh")
