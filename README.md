# pyv4s

This package provides a Python interface for scientists to seamlessly perform V<sub>4S</sub> calculations using high-performance C++ enabling easier manipulation, analysis, and visualization of the results.

For a detailed explanation of V<sub>4S</sub> index visit https://github.com/nicolas-loubet/V4S.

If you use the program for research, please cite: 
- Eur. Phys. J. E (2024) 47, 61 (https://doi.org/10.1140/epje/s10189-024-00454-3).
- J. Chem. Phys. (2024) 160, 144502 (https://doi.org/10.1063/5.0203989).

## Features

- **High-Performance Backend:** Leverages optimized C++ implementations for efficient computation of tetrahedral structures and molecular interactions.
- **Python Interface:** Easy-to-use Python bindings for interacting with the C++ library.
- **Multi-Software:** Supports simulation data independently of the molecular dynamics platform, including GROMACS, LAMMPS, and others.
- **Customizable Calculations:** Flexible options for studying specific molecular interactions and properties.


## Installation

### Prerequisites
1. Python 3.6 or higher.
2. Required Python packages: `matplotlib`, `numpy`, `pandas`.
3. Build tools for compiling the C++ library.

### Setup
1. Clone the repository:
   ```bash copy
   git clone https://github.com/alerverde/pyv4s.git
   ```
2. Compile the C++ library:
   ```bash
   cd pyv4s
   g++ -shared -fPIC -o pyv4s/Cpp/libv4s.so pyv4s/Cpp/V4S.cpp -Ipyv4s/Cpp -lstdc++ -lm
   ```
3. Install Python package:
   ```bash
   pip install .
   ```
   
## Usage: Python Interface
The Python scripts interact with the C++ backend to process and analyze molecular dynamics data.

### `computeTetrahedrons`
Computes the perfect tetrahedron points for each water molecule.

**Parameters:**
  - **`calculator`**: Instance of the `pyv4s.V4S` class, which provides the interface to the C++ backend.
  - **`df`**: Pandas DataFrame containing molecular system data, including atom positions and types. For each row (atom), the minimal information (columns) is:
    - `a` identifier that would be compared with `type_O`,
    - `x` position,
    - `y` position,
    - `z` position,
    - `e` epsilon value of Lennard Jones (LJ) interaction,
    - `s` sigma value of LJ interaction,
    - `q` Coulombic charge.
  - **`bounds`**: Array defining the system boundaries (`x`, `y`, `z`).
  - **`atoms_per_water_molecule (int, optional)`**: Number of atoms per water molecule (default is 3).
  - **`R_CUT_OFF (float, optional)`**: Cutoff radius for computing interactions, typically 5 Å.
  - **`study_molecules (list, optional)`**: Boolean array specifying which molecules to analyze (default is calculated using `detectOxygenMolecules` with the `type_O` parameter).
  - **`type_O (str, optional)`**: Atom type used to define water oxygens (default is `"O"`).

**Return:**
 - **`list[list[list[float]]]`**: A 3D array-like structure with shape `(N, 5, 3)`, where:  
    - `N` is the number of `True` values in `study_molecules`,  
    - `5` represents the five points of the tetrahedron (being the first one the center, and the other four the vertices),  
    - `3` corresponds to the `x`, `y`, and `z` coordinates of each point.

<br>

### `computeV4S` 
Calculates the V<sub>4S</sub> index for each water molecule.

**Parameters:**
  - **`calculator`**: Instance of the `pyv4s.V4S` class, which provides the interface to the C++ backend.
  - **`df`**: Pandas DataFrame containing molecular system data, including atom positions and types. For each row (atom), the minimal information (columns) is:
    - `a` identifier that would be compared with `type_O`,
    - `x` position,
    - `y` position,
    - `z` position,
    - `e` epsilon value of LJ interaction,
    - `s` sigma value of LJ interaction,
    - `q` Coulombic charge.
  - **`bounds`**: Array defining the system boundaries (`x`, `y`, `z`).
  - **`atoms_per_water_molecule (int, optional)`**: Number of atoms per water molecule (default is 3).
  - **`R_CUT_OFF (float, optional)`**: Cutoff radius for computing interactions, typically 5 Å.
  - **`study_molecules (list, optional)`**: Boolean array specifying which molecules to analyze (default is calculated using `detectOxygenMolecules` with the `type_O` parameter).
  - **`type_O (str, optional)`**: Atom type used to define water oxygens (default is `"O"`).

**Return:**
- **`list[float]`**: An array containing the V<sub>4S</sub> value for each molecule, with size `N`, where `N` is the number of `True` values in `study_molecules`.

<br>

### `computeViSPoints`
Calculates the positions of the four tetrahedron points and their corresponding V<sub>iS</sub> values (i = 1, 2, 3, 4).

**Parameters:**
  - **`calculator`**: Instance of the `pyv4s.V4S` class, which provides the interface to the C++ backend.
  - **`df`**: Pandas DataFrame containing molecular system data, including atom positions and types. For each row (atom), the minimal information (columns) is:
    - `a` identifier that would be compared with `type_O`,
    - `x` position,
    - `y` position,
    - `z` position,
    - `e` epsilon value of LJ interaction,
    - `s` sigma value of LJ interaction,
    - `q` Coulombic charge.
  - **`bounds`**: Array defining the system boundaries (`x`, `y`, `z`).
  - **`atoms_per_water_molecule (int, optional)`**: Number of atoms per water molecule (default is 3).
  - **`R_CUT_OFF (float, optional)`**: Cutoff radius for computing interactions, typically 5 Å.
  - **`study_molecules (list, optional)`**: Boolean array specifying which molecules to analyze (default is calculated using `detectOxygenMolecules` with the `type_O` parameter).
  - **`type_O (str, optional)`**: Atom type used to define water oxygens (default is `"O"`).

**Return:**
 - **`list[list[list[float]]]`**: A 3D array-like structure with shape `(N, 4, 4)`, where:  
    - `N` is the number of `True` values in `study_molecules`,  
    - `4` is the number of energy contributions sorted in descending order,
    - `4` corresponds to the `x`, `y`, and `z` coordinates of each point and the corresponding V<sub>iS</sub> (sum of interactions in that point).

<br>

---

### Example: Computing V<sub>4S</sub> Values

First import py4vs package and matplotlib 

```python
import handling.GROMACS
import handling.VMD
import pyv4s.V4S
import matplotlib.pyplot as plt
```

In this example we will use handmade GROMACS functions to extract the data, calculate V<sub>4S</sub> and save the results
```python
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
```
![Medium points](test/test_GROMACS.png?raw=true "Histogram")

## C++ Library
The core C++ library implements high-performance algorithms, exposed via Python bindings.

### Key C++ Functions
In case you want to interact directly with the C++ library (not recommended).

**Notes for every function:**
  - Hydrogen atoms and others (for models with more than 3 atoms) must follow the corresponding oxygen atom in the `Atom*` list.
  - Only water oxygens can have `atom_list[i].a = type_O`.

<br>

### `tetrahedrons`
Calculates perfect tetrahedron points for every water molecule.

**Parameters:**
  - **`atom_list (Atom*)`**: Pointer to the first member of an array of atoms that make up the system. Use the python class `Atom`.
  - **`N_ATOMS (int)`**: Total number of atoms in the system, indicating the end point of the previous array.
  - **`bounds_f (float*)`**: Array defining the system boundaries (`x`,`y`,`z`).
  - **`study_molecules (bool*)`**: Boolean array specifying which molecules to analyze. We recommend create it using `detectOxygenMolecules` python function, and then filtering the molecules of interest.
  - **`type_O (const char*)`**: Atom type used to define the water' oxygens.
  - **`atoms_per_water_molecule (int)`**: Number of atoms per water molecule (3 for TIP3P, 4 for TIP4P/2005, 5 for TIP5P, and so on...).
  - **`R_CUT_OFF (float)`**: Cutoff radius for computing interactions (normally 5 Å).

**Return:**
  - **`float***`**: A 3D array (`N`, `5`, `3`) containing the coordinates of tetrahedron centers and sites for each molecule. The first dimension represents the number of `True` values in `study_molecules`, the second one are the 5 points of the tetrahedron (the first one is the oxygen used as center, while the other 4 are the different points that make up the tetrahedron), and the third one are the coordinates (`x`,`y`,`z`).
  - For releasing memory, use `freeMemoryTet`.

<br>

### `V4S`
Calculates the V<sub>4S</sub> index for every water molecule.

**Parameters:**
  - **`atom_list (Atom*)`**: Pointer to the first member of an array of atoms that make up the system. Use the python class `Atom`.
  - **`N_ATOMS (int)`**: Total number of atoms in the system, indicating the end point of the previous array.
  - **`bounds_f (float*)`**: Array defining the system boundaries (`x`,`y`,`z`).
  - **`study_molecules (bool*)`**: Boolean array specifying which molecules to analyze. We recommend create it using `detectOxygenMolecules` python function, and then filtering the molecules of interest.
  - **`type_O (const char*)`**: Atom type used to define the water' oxygens.
  - **`atoms_per_water_molecule (int)`**: Number of atoms per water molecule (3 for TIP3P, 4 for TIP4P/2005, 5 for TIP5P, and so on...).
  - **`R_CUT_OFF (float)`**: Cutoff radius for computing interactions (normally 5 Å).

**Return:**
  - **`float*`**: A 1D array containing the values of the V<sub>4S</sub> index for each water molecule flaged in the `study_molecules` array. The array size is the number of `True` values in the `study_molecules` array.
  - For releasing memory, use `freeMemory`.

<br>

### `ViSPoints`
Calculates the position of the four points of the tetrahedron and the corresponding V<sub>iS</sub> ( with i={1,2,3,4} ).

**Parameters:**
  - **`atom_list (Atom*)`**: Pointer to the first member of an array of atoms that make up the system. Use the python class `Atom`.
  - **`N_ATOMS (int)`**: Total number of atoms in the system, indicating the end point of the previous array.
  - **`bounds_f (float*)`**: Array defining the system boundaries (`x`,`y`,`z`).
  - **`study_molecules (bool*)`**: Boolean array specifying which molecules to analyze. We recommend create it using `detectOxygenMolecules` python function, and then filtering the molecules of interest.
  - **`type_O (const char*)`**: Atom type used to define the water' oxygens.
  - **`atoms_per_water_molecule (int)`**: Number of atoms per water molecule (3 for TIP3P, 4 for TIP4P/2005, 5 for TIP5P, and so on...).
  - **`R_CUT_OFF (float)`**: Cutoff radius for computing interactions (normally 5 Å).

**Return:**
  - **`float***`**: A 3D array (`N`, `4`, `4`) containing the coordinates of tetrahedron points and the corresponding V<sub>iS</sub> for each molecule. The first dimension represents the number of `True` values in `study_molecules`, the second one are the 4 points of the tetrahedron, and the third one are the positions and V<sub>iS</sub> value (`x`,`y`,`z`,`V_iS`), where `V_iS` is the V<sub>1S</sub> at the first point in the array, and the V<sub>4S</sub> at the last one.
  - For releasing memory, use `freeMemoryTet`.

<br>

### `freeMemory`
Releases memory allocated for 1D arrays.

**Parameters:**
  - **`pointer (float*)`**: Pointer to the allocated memory.

<br>

### `freeMemoryTet`
Releases memory allocated for 3D arrays.

**Parameters:**
  - **`pointer (float***)`**: Pointer to the allocated 3D array (tetrahedron).
  - **`N_ATOMS (int)`**: Total number of atoms in the system (first dimension).
  - **`N_POINTS (int)`**: Number of points per tetrahedron (second dimension).

<br>

## License
This code repository is licensed under the MIT License. See the LICENSE file for details.


## Contributing
This is an open-source project, and we welcome contributions of all kinds. You’re encouraged to open an issue, submit a pull request, or reach out to us via email.


## Authors
- Alejandro R. Verde - alerverde@gmail.com
- Nicolás A. Loubet - loubetnicolasalfredo@gmail.com
