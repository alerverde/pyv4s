import handling.GROMACS
import handling.VMD
import pyv4s.V4S
import matplotlib.pyplot as plt


conf, bounds= handling.GROMACS.read_coordinates("Conf/GROMACS_Bulk/conf.gro")
conf= handling.GROMACS.read_topology("Conf/GROMACS_Bulk/topol.top", conf)
V4S_calculator= pyv4s.V4S.startV4S()

list_tet= pyv4s.V4S.computeTetrahedrons(V4S_calculator, conf, bounds, type_O="OW", atoms_per_water_molecule=4)
VMD_tet= ""
for t in list_tet:
	for site in t[1:]:
		VMD_tet+= "draw line "+handling.VMD.point_to_string(t[0])+" "+handling.VMD.point_to_string(site)+"\n"
with open("VMD_tetrahedrons.txt","w") as f:
	f.write(VMD_tet)

list_v4s= pyv4s.V4S.computeV4S(V4S_calculator, conf, bounds, atoms_per_water_molecule=4, R_CUT_OFF=5.0, study_molecules=None, type_O="OW")

plt.hist(list_v4s, bins=50, color="blue", edgecolor="black", alpha=0.5)
plt.title(r"$V_{4S}~index~histogram$")
plt.xlabel(r"$V_{4S}/kJ~{mol}^{-1}$")
plt.ylabel("Frequency")
plt.savefig('test_GROMACS.png')
