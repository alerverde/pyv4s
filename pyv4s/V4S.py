import ctypes
import os


class Atom(ctypes.Structure):
	_fields_= [("a", ctypes.c_char_p),("x", ctypes.c_float),("y", ctypes.c_float),("z", ctypes.c_float),
	("e", ctypes.c_float),("s", ctypes.c_float),("q", ctypes.c_float)]

def startV4S():
	f= ctypes.CDLL( os.path.join(os.path.dirname(__file__), "Cpp", "libv4s.so") )
	params= [ctypes.POINTER(Atom),ctypes.c_int,ctypes.POINTER(ctypes.c_float),ctypes.POINTER(ctypes.c_bool),ctypes.c_char_p,ctypes.c_int,ctypes.c_float]
	res= ctypes.POINTER(ctypes.POINTER(ctypes.POINTER(ctypes.c_float)))
	f.V4S.argstype,f.V4S.restype= params,ctypes.POINTER(ctypes.c_float)
	f.tetrahedrons.argstype,f.tetrahedrons.restype= params,res
	f.ViSPoints.argstype,f.ViSPoints.restype= params,res
	f.freeMemory.argstype,f.freeMemory.restype= [ctypes.POINTER(ctypes.c_float)],None
	f.freeMemoryTet.argstype,f.freeMemoryTet.restype= [res,ctypes.c_int,ctypes.c_int],None
	return f

def detectOxygenMolecules(df, type_O="O"):
	return (df["a"] == type_O).tolist()

def __prepareCall(df, bounds, study_molecules, type_O):
	if(study_molecules is None): study_molecules= detectOxygenMolecules(df,type_O)
	atoms= [Atom(row["a"].encode("utf-8"),row["x"],row["y"],row["z"],row["e"],row["s"],row["q"]) for i,row in df.iterrows()]
	N= len(atoms)
	arr_table= (Atom * len(atoms))(*atoms)
	n_wat_molecules= sum(study_molecules)
	return ctypes.byref((ctypes.c_bool * N)(*study_molecules)),N,ctypes.byref(arr_table),n_wat_molecules,ctypes.byref((ctypes.c_float * 3)(*bounds)),ctypes.c_char_p(type_O.encode('utf-8'))

def __passArray(pointer,size):
	return [[[pointer[k][j][i] for i in range(size[0])] for j in range(size[1])] for k in range(size[2])]

def computeV4S(calculator, df, bounds, atoms_per_water_molecule=3, R_CUT_OFF=5.0, study_molecules=None, type_O="O"):
	atoms,N_atoms,arr_table,n_wat_molecules,bounds,type_O= __prepareCall(df,bounds,study_molecules,type_O)
	pointer= calculator.V4S(arr_table,N_atoms,bounds,atoms,type_O,atoms_per_water_molecule,ctypes.c_float(R_CUT_OFF))
	list_V4S= [pointer[i] for i in range(n_wat_molecules)]
	calculator.freeMemory(pointer)
	return list_V4S

def computeTetrahedrons(calculator, df, bounds, atoms_per_water_molecule=3, R_CUT_OFF=5.0, study_molecules=None, type_O="O"):
	atoms,N_atoms,arr_table,n_wat_molecules,bounds,type_O= __prepareCall(df,bounds,study_molecules,type_O)
	pointer= calculator.tetrahedrons(arr_table,N_atoms,bounds,atoms,type_O,atoms_per_water_molecule,ctypes.c_float(R_CUT_OFF))
	list_tetrahedrons= __passArray(pointer,(3,5,n_wat_molecules))
	calculator.freeMemoryTet(pointer,(ctypes.c_int)(n_wat_molecules),(ctypes.c_int)(5))
	return list_tetrahedrons

def computeViSPoints(calculator, df, bounds, atoms_per_water_molecule=3, R_CUT_OFF=5.0, study_molecules=None, type_O="O"):
	atoms,N_atoms,arr_table,n_wat_molecules,bounds,type_O= __prepareCall(df,bounds,study_molecules,type_O)
	pointer= calculator.ViSPoints(arr_table,N_atoms,bounds,atoms,type_O,atoms_per_water_molecule,ctypes.c_float(R_CUT_OFF))
	list_tetrahedrons= __passArray(pointer,(4,4,n_wat_molecules))
	calculator.freeMemoryTet(pointer,(ctypes.c_int)(n_wat_molecules),(ctypes.c_int)(4))
	return list_tetrahedrons
