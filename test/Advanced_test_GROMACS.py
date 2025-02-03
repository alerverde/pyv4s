import handling.GROMACS
import pyv4s.V4S
import numpy as np

def write_PDB_line(i, name, code_molec, molec_num, x, y, z, beta):
	return f"ATOM  {i:5d} {name:<4} {code_molec} {molec_num:5d}    {x:8.3f}{y:8.3f}{z:8.3f}  1.00{beta:7.2f}\n"

conf, bounds= handling.GROMACS.read_coordinates("Conf/GROMACS_Graphene/conf.gro")
conf= handling.GROMACS.read_topology("Conf/GROMACS_Graphene/topol.top", conf)

V4S_calculator= pyv4s.V4S.startV4S()

list_Ox_study= pyv4s.V4S.detectOxygenMolecules(conf,type_O="OW")
list_Ox_identity= []
for i in range(len(conf)):
	if(not list_Ox_study[i]): continue
	if(np.abs(conf["z"][i]-34.42) > 3.5): list_Ox_study[i]= False # 3.5A in z-dir
	if(((conf["x"][i] - 59.02)**2 + (conf["y"][i] - 58.97)**2) > 12**2): list_Ox_study[i]= False # Cylinder of radius 12 center in CoM
	if(list_Ox_study[i]): list_Ox_identity.append(i)

list_tet= pyv4s.V4S.computeViSPoints(V4S_calculator,conf,bounds,atoms_per_water_molecule=4,study_molecules=list_Ox_study,type_O="OW")

with open("output.pdb", "w") as f:
	i_continue= 1
	f.write("HEADER    V4S in Graphene Sheet\n")
	f.write(f"CRYST1 {bounds[0]:9.3f}{bounds[1]:9.3f}{bounds[2]:9.3f}{90:7.2f}{90:7.2f}{90:7.2f}\n")
	for i_wat,i,tet in zip(range(len(list_tet)),list_Ox_identity,list_tet):
		for atom,atom_type in zip(range(i,i+4),["OW","HW1","HW2","MW"]):
			f.write(write_PDB_line(i_continue,atom_type,"WAT",i_wat+1,conf["x"][atom],conf["y"][atom],conf["z"][atom],0))
			i_continue+= 1
		for i_site,site in enumerate(tet):
			x,y,z= site[:3]
			vis= site[3]
			f.write(write_PDB_line(i_continue,"S"+str(i_site+1),"SIT",i_wat+1,x,y,z,vis))
			i_continue+= 1
	f.write("END\n")
