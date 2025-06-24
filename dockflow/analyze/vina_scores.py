
def scores_from_pdbqt(pdbqt_path: str = "vina_out.pdbqt"):
    scores=[]
    with open(pdbqt_path, "r") as f:
        for line in f:
            if line.startswith("REMARK VINA RESULT"):
                scores.append(float(line.split()[3]))
    return scores