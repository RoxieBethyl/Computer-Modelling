**Author**: Blythe Fernandes

# Project Description

In the main code file, there are 2 tests.
"Test 1" uses 2 particles to check for errors in mathematical equations and errors in code at the early stages of debugging. The test conditions for this test are:
 	Number of particles = 2
	Step size (time) = 0.01
	numstep = 1000

"Test 2" evaluates data for particles that fit a lattice of size 4n^3 where n is the number of particles. (n could be 4, 32, 108). The parameters for test are given by user-input. The user is promted to input the intial parameters. Typical parameters used are:
	**Solid (Solid Argon):**
		- number of particles = 32
		- density = 1.
		- temperature = 0.1
		- step size (time) = 0.01
		- numstep = 1000		

	**Gas (Gaseous Argon):**
		- number of particles = 30
		- density = 0.05
		- temperature = 1.
		- step size (time) = 1.
		- numstep = 2000

## Project Files
### Files within this folder:
### Code files:
 - `mdutilities.py`
 - `particle3D.py`
 - `pbc.py`
 - `ProjectB.py` (Main code file)

### Data files (that are created when the code is run):
 - `msd.dat` 
	The data is this file is stored as (time msd)
 - `rdf.dat`
	Data is stored as (time \[distance\] \[RDF\])
 - `outfile.xyz` (This file stores the particles position at each timestep to be used in VMD)

### HTML files (for documentation):
 - `particle3D.html`
 - `pbc.html`
 - `ProjectB.html`

**Note:** The files with `particle3D` and `pbc` in the filename (both `.py` and `.html`) are taken previous exercises done for Computer Modelling.

