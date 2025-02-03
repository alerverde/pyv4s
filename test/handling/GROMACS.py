import pandas as pd

def read_lines(dir):
    with open(dir,"r") as f:
        return f.readlines()

def read_coordinates(dir):
    lines = read_lines(dir)
    n_atoms= int(lines[1].strip())
    atoms_df= pd.read_fwf(dir, widths= [5,5,5,5,8,8,8], skiprows=2, nrows=n_atoms, names= ["molec_num", "molec", "a", "id", "x", "y", "z"])
    atoms_df[["x","y","z"]] = atoms_df[["x","y","z"]].astype(float)
    atoms_df[["molec_num","id"]] = atoms_df[["molec_num","id"]].astype(int)
    atoms_df[["x","y","z"]]= atoms_df[["x","y","z"]]*10
    bounds= list(map(float,lines[-1].strip().split()))
    bounds= [b*10 for b in bounds]
    return atoms_df, bounds

def clean_line(line):
    line = line.strip()
    if not line or line.startswith("["):
        return None
    return line.split()

def read_LJ(lines, i_atomtypes):
    atomtypes_dict = {}
    for line in lines[i_atomtypes:]:
        values = clean_line(line)
        if not values:
            break
        type= values[0]
        s= float(values[5])
        e= float(values[6])
        atomtypes_dict[type]= {"s":s,"e":e}
    return atomtypes_dict

def read_atom_charge(lines,i_atomtypes):
    atoms_dict = {}
    i_atoms= i_atomtypes+1
    while("[ atoms ]\n" in lines[i_atoms:]):
        i_atoms+= lines[i_atoms:].index("[ atoms ]\n") + 2
        for line in lines[i_atoms:]:
            values = clean_line(line)
            if not values:
                break
            at_type= values[1]
            name= values[4]
            q= float(values[6])
            atoms_dict[name]= {"a":at_type,"q":q}
    return atoms_dict

def combine_dicts(atoms_df,atomtypes_dict,atoms_dict):
    atoms_info= pd.DataFrame.from_dict(atoms_dict,orient="index")
    atoms_es= pd.DataFrame.from_dict(atomtypes_dict,orient="index")
    atoms_info["e"]= atoms_info["a"].map(atoms_es["e"])
    atoms_info["s"]= atoms_info["a"].map(atoms_es["s"])
    atoms_info= atoms_info.drop(columns=["a"]).reset_index()
    atoms_info.rename(columns={"index":"a"},inplace=True)
    df_combined = atoms_df.merge(atoms_info[["a","e","s","q"]], on="a", how="left")
    df_filtered = df_combined.drop(columns=["molec_num","molec","id"])
    return df_filtered

def read_topology(dir, atoms_df):
    lines = read_lines(dir)    
    i_atomtypes = lines.index("[ atomtypes ]\n") + 2
    atomtypes_dict = read_LJ(lines,i_atomtypes)
    atoms_dict = read_atom_charge(lines,i_atomtypes)
    atoms_df_modified = combine_dicts(atoms_df,atomtypes_dict,atoms_dict) 
    return atoms_df_modified
