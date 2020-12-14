def getParams(cfg):
    import os
    mainParams = ["complex", "rbd", "posRange", "mutTake", "repair", "posscan", "bmRBD", "bmComplex", "analyseComplex"]
    with open(cfg) as f:
        params = {x[0].strip(): x[1].strip() for x in [l.strip().split("=") for l in f.readlines()]}
    assert sum([p not in params.keys() for p in mainParams]) == 0
    assert params["complex"].endswith(".pdb")
    assert params["rbd"].endswith(".pdb")
    dataList = os.listdir("./data")
    assert params["complex"] in dataList
    assert params["rbd"] in dataList
    low, upp = params["posRange"].split(":")
    try:
        low = int(low)
        upp = int(upp)
    except ValueError as err:
        print(err)
        assert 1 == 0
    mutTake = params["mutTake"]
    try:
        mutTake = int(mutTake)
    except ValueError as e1:
        try:
            mutTake = float(mutTake)
        except ValueError as e2:
            assert mutTake in ["input", "all"]
    assert params["repair"] in ["todo", "done"]
    assert params["posscan"] in ["todo", "done", "skip"]
    assert params["bmRBD"] in ["todo", "done", "skip"]
    assert params["bmComplex"] in ["todo", "done", "skip"]
    assert params["analyseComplex"] in ["todo", "done", "skip"]
    return params


def makeRepair(**params):
    import os
    structure = params["complex"]
    structureRep = params["complex"].split(".")[0]+"_Repair.pdb"
    repairPath = "./stages/1_RepairPDB"
    if params["repair"] == "todo":
        print(f"\n\nRepairing of {structure} structure is in progress...\n\n")
        os.system(f"cp ./data/{structure} {repairPath}")
        os.chdir(repairPath)
        command = f"../../fxbin/foldx --command=RepairPDB --pdb={structure}"
        os.system(command)
        os.system(f"cp ./{structureRep} ../../data")
    else:
        print(f"\n\nRepairing of {structure} structure was already done\n\n")


def getDist(seq1, seq2):
    assert len(seq1) == len(seq2)
    result = 0
    for i in range(len(seq1)):
        result += seq1[i] != seq2[i]
    return result


def makePositionScan(**params):
    import os
    import pickle
    import numpy as np
    import pandas as pd
    from subprocess import Popen, PIPE
    from itertools import combinations, product
    nthreads = int(params["nthreads"])
    structure = params["rbd"]
    if params["posscan"] == "done":
        print(f"\n\nPosition scanning in {structure} structure was already done\n\n")
    elif params["posscan"] == "skip":
        print(f"\n\nPosition scanning in {structure} structure is skipped\n\n")
    else:
        print(f"\n\nPosition scanning in {structure} structure is in progress...\n")
        positions = params["positions"].split(",")
        finalPositions = []
        # справочные словари по обозначениям и соответствию кодов и АК
        codons, one2three, three2one = pickle.load(open("./scripts/gencode.pkl", "rb"))
        # списки АК по флагам Position Scan
        psflags = pickle.load(open("./scripts/psflags.pkl", "rb"))
        # нуклеотидная и аминокислотная последовательности S-белка
        geneS, protS = pickle.load(open("./scripts/Sseq.pkl", "rb"))
        for position in positions:
            chain = position[1]
            index = position[2:-1]
            assert protS[int(index) - 1] == position[0]
            wtRes1 = position[0]
            wtCodon = geneS[int(index) * 3 - 3:int(index) * 3]
            flag = position[-1]
            assert flag in psflags.keys() or flag in one2three.keys()
            try:
                mutResList = psflags[flag]
            except KeyError as err:
                mutResList = [one2three[flag]]
            for mutRes in mutResList:
                mutCodons = codons[mutRes]
                for mutCodon in mutCodons:
                    dist = getDist(wtCodon, mutCodon)
                    if dist == 1:  # отбираем только SNP
                        finalPositions.append(wtRes1 + chain + index + three2one[mutRes])
                        break
        print(f"{len(finalPositions)} point mutations was selected:\n")
        print(",".join(finalPositions))
        equalN = len(finalPositions)//nthreads
        binnedPositions = [[]]*nthreads
        for i in range(nthreads):
            mut = finalPositions[:equalN]
            binnedPositions[i] = mut
            finalPositions = finalPositions[equalN:]
        i = 0
        while finalPositions:
            mut = finalPositions[0]
            binnedPositions[i].append(mut)
            finalPositions = finalPositions[1:]
            i += 1
        binnedPositions = [",".join(posList) for posList in binnedPositions]
        projectPath = os.getcwd()
        psPath = os.path.join(projectPath, "stages", "2_PositionScan")
        threadPaths = []
        for i in range(nthreads):
            threadPath = f"{psPath}/thread_{i+1}"
            threadPaths.append(threadPath)
            os.system(f"mkdir {threadPath}")
            os.system(f"cp ./data/{structure} {threadPath}")
        cdCmdList = [f"cd {threadPath}" for threadPath in threadPaths]
        fxCmdList = [
            f"../../../fxbin/foldx --command=PositionScan --pdb={structure} --positions={ps}" for ps in binnedPositions
        ]
        cmdList = [" && ".join([cdCmdList[i], fxCmdList[i]]) for i in range(nthreads)]
        procList = [Popen(cmd, stdout=PIPE, stderr=PIPE, shell=True) for cmd in cmdList]
        for i in range(nthreads):
            process = procList[i]
            stderr = process.stderr.read()
            if stderr:
                print(stderr)
            process.wait()
        os.chdir(psPath)
        output = None
        for i in range(nthreads):
            if output is None:
                output = pd.read_csv(f"""./thread_{i+1}/PS_{structure.split(".")[0]}_scanning_output.txt""", sep="\t", header=None)
            else:
                output = pd.concat(
                    [output, pd.read_csv(f"""./thread_{i+1}/PS_{structure.split(".")[0]}_scanning_output.txt""", sep="\t", header=None)]
                )
        output.columns = ["mutation", "diff"]
        if output.iloc[-1, 0][:3] not in three2one.keys():
            output = output.iloc[:-1, :]
        output["command"] = output["mutation"].apply(lambda x: three2one[x[:3]] + x[3:])
        output["position"] = output["mutation"].apply(lambda x: x[4:-1])
        df = output[output["diff"] < 0]  # отбираем только замены, не уменьшающие стабильность RBD
        df = df[df["command"].apply(lambda x: x[0] != x[-1])]
        positions = np.unique(df.position)
        posCombinations = []
        low, upp = params["posRange"].split(":")
        low = int(low)
        upp = int(upp)
        upp = min(len(positions), upp)
        for i in range(low, upp + 1):
            for x in combinations(positions, i):
                posCombinations.append(list(x))
        mutCombinations = []
        for posCombination in posCombinations:
            mutVariants = []
            for pos in posCombination:
                mutVariants.append(list(df[df.position == pos]['command'].values))
            mutVariants = product(*mutVariants)
            for mutVariant in mutVariants:
                mutCombinations.append(list(mutVariant))
        print(f"\nTotal number of mutations sets: {len(mutCombinations)}\n")
        if params["mutTake"] == "all":
            commands = '\n'.join([';'.join(x) + ';' for x in mutCombinations])
            with open("../3_BuildModel_RBD/individual_list.txt", "w") as f:
                f.write(commands)
            print("Admissible mutation sets:\n")
            os.system("cat ../3_BuildModel_RBD/individual_list.txt")
        else:
            if params["mutTake"] == "input":
                mutTake = input("Input number or fraction of mutations to analyse: ")
            else:
                mutTake = params["mutTake"]
            try:
                mutTake = min(int(mutTake), len(mutCombinations))
                finMutCombinations = np.random.choice([';'.join(x) + ';' for x in mutCombinations], size=mutTake)
                commands = '\n'.join(finMutCombinations)
                with open("../3_BuildModel_RBD/individual_list.txt", "w") as f:
                    f.write(commands)
                print("Admissible mutations sets:\n")
                os.system("cat ../3_BuildModel_RBD/individual_list.txt")
            except ValueError as e1:
                try:
                    mutTake = float(mutTake)
                    assert 0 < mutTake < 1
                    finMutCombinations = np.random.choice(
                        [';'.join(x) + ';' for x in mutCombinations], size=1+int(len(mutCombinations)*mutTake)
                    )
                    commands = '\n'.join(finMutCombinations)
                    with open("../3_BuildModel_RBD/individual_list.txt", "w") as f:
                        f.write(commands)
                    print("Admissible mutations sets:\n")
                    os.system("cat ../3_BuildModel_RBD/individual_list.txt")
                except ValueError as e2:
                    print("Bad input")


def buildModel(obj, **params):
    assert obj in ['rbd', 'complex']
    import os
    import numpy as np
    import pandas as pd
    from subprocess import Popen, PIPE
    if obj == "rbd":
        structure = params["rbd"]
    else:
        structure = params["complex"].split(".")[0]+"_Repair.pdb"
    projectPath = os.getcwd()
    reportPath = os.path.join(projectPath, "reports")
    if obj == "rbd":
        bmPath = os.path.join(projectPath, "stages", "3_BuildModel_RBD")
    else:
        bmPath = os.path.join(projectPath, "stages", "4_BuildModel_complex")
    nthreads = int(params["nthreads"])
    flag = params["bmRBD"] if obj == "rbd" else params["bmComplex"]
    if flag == "todo":
        print(f"\n\nModel building for mutated {structure} structure is in progress...\n")
        with open(f"{bmPath}/individual_list.txt", "r") as f:
            mutations = f.readlines()
        mutations[-1] = mutations[-1] + "\n"
        mutations = list(np.random.choice(mutations, size=len(mutations), replace=False))
        equalN = len(mutations)//nthreads
        binnedMutations = [[]]*nthreads
        for i in range(nthreads):
            mut = mutations[:equalN]
            binnedMutations[i] = mut
            mutations = mutations[equalN:]
        i = 0
        while mutations:
            mut = mutations[0]
            binnedMutations[i].append(mut)
            mutations = mutations[1:]
            i += 1
        threadPaths = []
        thread2task = {}
        for i in range(nthreads):
            threadPath = f"{bmPath}/thread_{i+1}"
            threadPaths.append(threadPath)
            thread2task[threadPath] = []
            os.system(f"mkdir {threadPath}")
            for j in range(len(binnedMutations[i])):
                taskPath = os.path.join(threadPath, f"task_{j+1}")
                thread2task[threadPath].append(taskPath)
                os.system(f"mkdir {taskPath}")
                os.system(f"cp ./data/{structure} {taskPath}")
                with open(f"{taskPath}/individual_list.txt", "w") as f:
                    f.write(binnedMutations[i][j].strip())
        cdCmdList = []
        for threadPath in threadPaths:
            for taskPath in thread2task[threadPath]:
                cdCmdList.append(f"cd {taskPath}")
        fxCmd = f"../../../../fxbin/foldx --command=BuildModel --pdb={structure} --mutant-file=individual_list.txt"
        cmdList = [" && ".join([cdCmdList[i], fxCmd]) for i in range(len(cdCmdList))]
        while cmdList:
            todoList = cmdList[:nthreads]
            procList = [Popen(cmd, stdout=PIPE, stderr=PIPE, shell=True) for cmd in todoList]
            for i in range(len(todoList)):
                process = procList[i]
                stderr = process.stderr.read()
                if stderr:
                    print(stderr)
                process.wait()
            cmdList = cmdList[nthreads:]
        for i in range(nthreads):
            threadPath = threadPaths[i]
            for j in range(len(thread2task[threadPath])):
                taskPath = thread2task[threadPath][j]
                os.chdir(taskPath)
                tag = structure.split(".")[0]
                with open(f'Dif_{tag}.fxout', 'r') as f:
                    lines = f.readlines()
                with open("tmp.txt", "w") as f:
                    f.write(''.join(lines[8:]))
                diff = pd.read_csv("tmp.txt", sep="\t")
                os.system("rm tmp.txt")
                with open('individual_list.txt', 'r') as f:
                    mutations = [x.strip() for x in f.readlines()]
                mutations = [m for m in mutations if m != ""]
                diff.insert(loc=0, column="Mut", value=mutations)
                diff.insert(loc=1, column="Thread", value=[i+1]*len(mutations))
                diff.insert(loc=2, column="Task", value=[j + 1] * len(mutations))
                if obj == "rbd":
                    getNum = lambda x: int(x.split("_")[-1].split(".")[0]) - 1
                    mutIdxs = diff[diff["total energy"] < 0]['Pdb'].apply(getNum).values
                    approved = [mutations[i] for i in mutIdxs]
                    if approved:
                        with open("../../../4_BuildModel_complex/individual_list.txt", "a") as f:
                            f.write('\n'.join(approved)+"\n")
                else:
                    diff.to_csv(os.path.join(reportPath, f"bm_complex_report_{i+1}_{j+1}.csv"))
        if obj == "rbd":
            with open("../../../4_BuildModel_complex/individual_list.txt", "r") as f:
                selected = f.readlines()
            print("Number of stable mutations sets:", len(selected))
            print("\n")
        else:
            diff = None
            for i in range(nthreads):
                threadPath = threadPaths[i]
                for j in range(len(thread2task[threadPath])):
                    if diff is None:
                        diff = pd.read_csv(os.path.join(reportPath, f"bm_complex_report_{i+1}_{j+1}.csv"))
                    else:
                        diff = pd.concat(
                            [diff, pd.read_csv(os.path.join(reportPath, f"bm_complex_report_{i+1}_{j+1}.csv"))]
                        )
            os.chdir(reportPath)
            os.system("rm ./*")
            diff.to_csv(os.path.join(reportPath, "bm_complex_report.csv"), index=False, index_label=range(len(diff)))
    elif flag == "done":
        print(f"\n\nModel building for mutated {structure} structure was already done\n\n")
    else:
        print(f"\n\nModel building for mutated {structure} structure is skipped\n\n")


def analyseComplex(**params):
    import os
    import pandas as pd
    from subprocess import Popen, PIPE
    structure = params["complex"].split(".")[0] + "_Repair.pdb"
    projectPath = os.getcwd()
    acPath = os.path.join(projectPath, "stages", "5_AnalyseComplex")
    nthreads = int(params["nthreads"])
    if params["analyseComplex"] == "todo":
        print(f"\n\nComplex analysis for mutated {structure} structure is in progress...\n\n")
        mainTag = structure.split(".")[0]
        os.system(f"""touch {os.path.join(acPath, "individual_list.txt")}""")
        counter = 1
        baseStructureMoved = False
        for i in range(nthreads):
            bmThreadPath = os.path.join(projectPath, "stages", "4_BuildModel_complex", f"thread_{i+1}")
            taskPaths = [os.path.join(bmThreadPath, f) for f in os.listdir(bmThreadPath) if f.startswith("task")]
            for taskPath in taskPaths:
                if not baseStructureMoved:
                    os.system(f"cp {os.path.join(taskPath, structure)} {acPath}")
                    baseStructureMoved = True
                sourceName = structure.split(".")[0]+f"_1.pdb"
                targetName = structure.split(".")[0]+f"_{counter}.pdb"
                os.system(f"""cp {os.path.join(taskPath, sourceName)} {os.path.join(acPath, targetName)}""")
                with open(os.path.join(taskPath, "individual_list.txt"), "r") as f:
                    mutations = [m for m in f.readlines() if m != "\n"]
                with open(os.path.join(acPath, "individual_list.txt"), "a") as f:
                    f.writelines("".join(mutations)+"\n")
                counter += 1
        os.chdir(acPath)
        tags = [f for f in os.listdir() if f.startswith(mainTag)]
        while tags:
            if len(tags) > nthreads:
                cur_tags = tags[:nthreads]
            else:
                cur_tags = tags
            cmdList = [f"../../fxbin/foldx --command=AnalyseComplex --pdb={tag} --analyseComplexChains=A,B" for tag in cur_tags]
            procList = [Popen(cmd, stdout=PIPE, stderr=PIPE, shell=True) for cmd in cmdList]
            for i in range(len(cur_tags)):
                process = procList[i]
                stderr = process.stderr.read()
                if stderr:
                    print(stderr)
                process.wait()
            tags = tags[nthreads:]

        tags = sorted([x[8:-9] for x in os.listdir() if x.startswith("Summary")])
        interaction = pd.DataFrame(
            index=[
                'Mutations', 'Pdb', 'Group1', 'Group2', 'IntraclashesGroup1', 'IntraclashesGroup2',
                'Interaction Energy', 'StabilityGroup1', 'StabilityGroup2']
        )
        with open('individual_list.txt', 'r') as f:
            mutations = [x.strip() for x in f.readlines()]
        for tag in tags:
            try:
                i = int(tag.split("_")[-1])
            except ValueError as err:
                i = 0
            with open(f'Summary_{tag}_AC.fxout', 'r') as f:
                lines = f.readlines()
            with open("tmp.txt", "w") as f:
                f.write(''.join(lines[8:]))
            summary = pd.read_csv("tmp.txt", sep="\t")
            os.system("rm tmp.txt")
            if i == 0:
                interaction[i] = ["Wild-Type"] + list(summary.T[0].values)
            else:
                interaction[i] = [mutations[i - 1]] + list(summary.T[0].values)
        interaction = interaction.T
        interaction.columns = ["Mutations"] + list(summary.columns)
        interaction.sort_values(by="Interaction Energy", ascending=True, inplace=True)
        interaction.reset_index(drop=True, inplace=True)
        interaction.to_csv("../../reports/screening_report.csv")
        print("\n\nMUTATION SCREENING REPORT\n")
        print(interaction)
        print("\n\n")

    elif params["analyseComplex"] == "done":
        print(f"\n\nComplex analysis for mutated {structure} structure was already done\n\n")
    else:
        print(f"\n\nComplex analysis for mutated {structure} structure is skipped\n\n")


