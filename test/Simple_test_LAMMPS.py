import handling.LAMMPS
import handling.VMD
import pyv4s.V4S
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np

conf= handling.LAMMPS.read_coordinates("Conf/LAMMPS/conf.xyz")
pot= handling.LAMMPS.read_LJ("Conf/LAMMPS/pot.txt","kJ")
sdata, bounds= handling.LAMMPS.parse_system_data("Conf/LAMMPS/system.data")
conf= handling.LAMMPS.merge_data(conf, pot, sdata)

stud_molec_numbers= pd.read_csv("Conf/LAMMPS/study_atoms.txt", header=None)[0].to_numpy()
stud_molec= np.zeros(len(conf), dtype=bool)
stud_molec[stud_molec_numbers]= True

V4S_calculator= pyv4s.V4S.startV4S()

list_tet= pyv4s.V4S.computeTetrahedrons(V4S_calculator, conf, bounds, study_molecules=stud_molec)
VMD_tet= ""
for t in list_tet:
	for site in t[1:]:
		VMD_tet+= "draw line "+handling.VMD.point_to_string(t[0])+" "+handling.VMD.point_to_string(site)+"\n"
with open("VMD_tetrahedrons.txt","w") as f:
	f.write(VMD_tet)


list_v4s= pyv4s.V4S.computeV4S(V4S_calculator, conf, bounds, study_molecules=stud_molec, atoms_per_water_molecule=4, type_O="OW")
plt.hist(list_v4s, bins=20, color="blue", edgecolor="black", alpha=0.5)
plt.title(r"$V_{4S}~index~histogram$")
plt.xlabel(r"$V_{4S}/kJ~{mol}^{-1}$")
plt.ylabel("Frequency")
plt.show()
