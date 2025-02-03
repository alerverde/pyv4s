import pandas as pd

def read_coordinates(dir):
    return pd.read_csv(dir, skiprows=2, sep=" ", header=None, names=["a","x","y","z"])

def read_LJ(dir, units="kJ"):
    df= pd.read_csv(dir, comment="#", sep=" ", skip_blank_lines=True, header=None, names=["a","e","s"])
    if(units == "kcal" or units == "Cal"):
        df["e"]*= 4.184
    return df

def parse_atoms_and_bounds(lines, atoms_start_index):
    n_atoms = 0
    bounds = [0.0, 0.0, 0.0]
    bounds_flags = ["xlo xhi","ylo yhi","zlo zhi"]
    for line in lines[:atoms_start_index - 1]:
        if "atoms" in line:
            n_atoms = int(line.split()[0])
        for i in range(3):
            if bounds_flags[i] in line:
                bounds[i] = float(line.split()[1]) - float(line.split()[0])
    return n_atoms, bounds

def parse_system_data(dir):
    with open(dir,"r") as f:
        lines = f.readlines()    
    atoms_start_index= next(i for i, line in enumerate(lines) if line.strip() == "Atoms") + 1
    n_atoms, bounds = parse_atoms_and_bounds(lines, atoms_start_index)
    atom_lines = []
    for line in lines[atoms_start_index:]:
        if(line.strip()):
            atom_lines.append(line.strip())
        if(len(atom_lines) == n_atoms):
            break
    data= pd.DataFrame([line.split() for line in atom_lines],columns=["i", "molecule-tag", "a_i", "q", "x", "y", "z"])
    data= data.astype({"i":int, "molecule-tag":int, "a_i":int, "q":float, "x":float, "y":float, "z":float})
    return data,bounds

def check_conf(conf,pot,sdata):
    if len(conf["a"]) != len(sdata["q"]):
        print("Charges in system data and atoms in configuration must have the same length")
        return False
    if not conf["a"].isin(pot["a"]).all():
        print("All atoms in configuration must be included in the LJ list")
        print("Found in conf: ", conf["a"].unique())
        print("Found in LJ list:", pot["a"].to_list())
        return False
    return True

def merge_data(conf, pot, sdata):
    if not check_conf(conf,pot,sdata):
        return None
    conf= conf.merge(pot[["a","e","s"]], on="a", how="left")
    conf["q"]= sdata["q"]
    return conf
