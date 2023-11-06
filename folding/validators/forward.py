import subprocess

"""
def gromacs(protein_pdb, ff, box = dodecahedron, energy_min_1, energy_min_2, energy_min_3, energy_min_ 4, verbose=False):
The function takes a .pdb file for a known sequence and structure of a protein and calculates the molecular
dyanmics of it folding without further user input using the GROMACS molecular dynamics software.
This is done in the following steps:

1) Creating the environment and solution the protein is folding in
2) Stabilizing the protein
3) Stabilizing the protein in the evironment

Notes:
Gromacs is typically run from the command line with heavy user input, the goal of this function is to either skip
the user input or use different servers to find the optimal input for us based on how "typcial" a run is for each protein

Gromacs works best for proteins folded in water. Other solutions need a different workflow
Given that this function makes generic names we may want to have each make a new folder with the proteins name?


Inputs (main):
protein_pdb: string arg for filepath of the .pdb file obtained from ___________
ff: forcefield, this varies by protein and is the second most important input. This is a great use of distributed computed
    and different servers using different force fields is the most optimal way to find this.
box: constrain the size of the simulated environment to reduce processing speed. Barring vastly irregular shapes,
    the rhombic dodecahedron is most optimal. This could be automated based on known size/shape of the protein or simply relgated to other servers
energy_min_x: a .mdp file describing the energy minimizatoin parameters. This can be iterated upon based on validation results
    in the tutorial, energy_min_1=emin-charmm.mdp and energy_min_2=nvt-charmm.mdp, energy_min_3=npt-charmm.mdp, energy_min_4=md-charmm.mdp


Inputs (optional, future):
solution: Non-water solutions will require different molecular dyanmics functions and workflos
verbose: Boolean descibing outputs for mdrun. If true, progress updates are ouput while dynamics are calculated
output: necessary output depends on method of transition state calculations. This can be easily changed

Final Outputs:
________.xvg



Notes:


First validaiton checkpoint: After mdrun. Although there are a lot of instances where we can step in for validation, this
is the most optimal one for determining run success

Second validation checkpoint: After temperature and pressure runs

Third validation checkpoint: Post analysis on any metric such as RMSD, radius of gyration, etc.

### CODE ###

# Strip out all the atoms in the file that are not the protein itself like water and ligands
# These are labeled with HEATM and this can be done with whatever text file processing method is preferred
# TODO: make
!grep -v HETATM test_protein.pdb > temp_clean_protein.pdb
!grep -v CONECT temp_clean_protein.pdb > clean_protein.pdb

# Input the file into GROMACS and get three output files: topology, position restraint, and a post-processed structure file
# TODO: Input is the forcefield "charmm27"
# TODO: This fails if atoms are missing from the .pdb file.
!gmx pdb2gmx -f clean_protein.pdb -o clean_protein.gro -water tip3p -ff "charmm27"


# Build the "box" to run our simulation of one protein molecule
# TODO: If the molecule isn't spherical the rhombic dodecahedron wont be ideal but this is the most common type
# The distance d (nm) establishes the distance from the edge of the protein to the edge of the simulated box.
!gmx editconf -f clean_protein.gro -o clean_protein_box.gro -c -d 1.0 -bt dodecahedron

# spc216.gro is included in gromacs and used in water tip3p solutions. This is the step that needs to change for non-water solution
!gmx solvate -cp clean_protein_box.gro -cs spc216.gro -o protein_solv.gro -p topol.top


# Now we add a more natural salinity to our water solution (0.15M NaCl)
# You need to create an empty .mdp file to run the next function. This can be added to to alter molecular dynamics equations
!touch ions.mdp
!gmx grompp -f ions.mdp -c protein_solv.gro -p topol.top -o ions.tpr

# The print statement is one of the user inputs mentioned above to interact with the program
!printf "SOL\n" | gmx genion -s ions.tpr -o protein_solv.gro -conc 0.15 -p \
topol.top -pname NA -nname CL -neutral


# This is another input that has increible importance on the protein folding. We are going to run energy minimization
# by first consoldating all our parameters and then running the minimization itself. The default works well in standard
# cases but there are a lot of knobs to turn

!cat energy_min_1
!gmx grompp -f energy_min_1 -c protein_solv.gro -p topol.top -o em.tpr
# -v can be added as a term to make this step print out progress. It normally takes a little while
!gmx mdrun -deffnm em


# From here the first stage of the run is complete. There are two main checks to determine a successful molecular dynamics run
# Potential energy. This should be around -100000kj/mol depending on protein size and amount of water.
# Maximum force. This is determined in the .mdp file as "emtol", typically 1000 kJ/(mol * nm). If Force is greater than this but
# but potential energy is correct than the system might be unstable and minmiation parameters in the .mdp file should be changed accordingly

# This is another instance of command line input being necessary
# The xvg file is pandas readable for analysis if wanted
!printf "Potential\n0\n" | gmx energy -f em.edr -o potential.xvg

# Right now the protein is stable with itself. Next we stabilize it within the actual solution its in, temperature, and pressure



# Temp run
# input second energy minimization parameters. Should include nVT (# particles, volume, temp)
!cat energy_min_2
!gmx grompp -f energy_min_2 -c em.gro -r em.gro -p topol.top -o nvt.tpr
# Feature: as before, -v as a term here outputs text while its running.
!gmx mdrun -ntmpi 1 -deffnm nvt


# Pressure run. Should include nPT (# particles, pressure, temp)
!cat energy_min_3
!gmx grompp -f energy_min_3 -c nvt.gro -r nvt.gro -t nvt.cpt -p topol.top -o npt.tpr
!gmx mdrun -ntmpi 1 -v -deffnm npt


# Validaiton checkpoint 2: Analysis of temperature and pressure
# The following functions make pandas friendly data files for analysis
#!echo "Temperature" | gmx energy -f nvt.edr -o temperature.xvg -b 20
# !echo "Pressure" | gmx energy -f npt.edr -o pressure.xvg -xvg none


# The production run
!cat energy_min_4
!gmx grompp -f energy_min_4 -c npt.gro -t npt.cpt -p topol.top -o md.tpr
# Again -v can be added here
!gmx mdrun -ntmpi 1 -deffnm md


# Technically the molecular dynamics calculations are done here. Further work is for specific analysis purposes


# RMSD data
# This function is an all-purpose post-processing tool that sets up the next function and accounts for periodic variations in the protein
!printf "1\n1\n" | gmx trjconv -s md.tpr -f md.xtc -o md_center.xtc -center -pbc mol
# Quantifies structural stability of the protein using RMSD. Useful for RMSD based models of protein folding
!printf "4\n1\n" | gmx rms -s em.tpr -f md_center.xtc -o rmsd_xray.xvg -tu ns -xvg none
"""


def forward(self,):
    # this is where the action happens


    # 1. Select the molecule from online protein database
    # NOTE: The number of possible inputs should be effectively infinite (or at least very large) so that miners cannot lookup results from earlier runs

    protein_pdb = 'https://files.rcsb.org/download/1UBQ.pdb' # can either be local file path or a url to download
    ff = 'default'
    box = 'dodecahedron'


    # 2. Preprocess the input files, cleaning up files and generating required inputs

    preprocess(protein_pdb, ff, box)

    # 3. Run first step locally

    # 4. Send the preprocessed inputs and other required details to the miners who will carry out the ful MD simulation

    # 5. Receive the results from the miners

    # 6. Validate/score the results


def preprocess(protein_pdb, ff, box):

    # Strip out all the atoms in the file that are not the protein itself like water and ligands
    # These are labeled with HEATM and this can be done with whatever text file processing method is preferred
    # TODO: make
    strip_atoms = 'grep -v HETATM test_protein.pdb > temp_clean_protein.pdb'
    subprocess.run(strip_atoms.split(), shell=True, check=True)

    rejoin_atoms = 'grep -v CONECT temp_clean_protein.pdb > clean_protein.pdb'
    subprocess.run(rejoin_atoms.split(), shell=True, check=True)

    # Input the file into GROMACS and get three output files: topology, position restraint, and a post-processed structure file
    # TODO: Input is the forcefield "charmm27"
    # TODO: This fails if atoms are missing from the .pdb file.
    create_env = 'gmx pdb2gmx -f clean_protein.pdb -o clean_protein.gro -water tip3p -ff "charmm27"'
    subprocess.run(create_env.split(), shell=True, check=True)

    # Build the "box" to run our simulation of one protein molecule
    # TODO: If the molecule isn't spherical the rhombic dodecahedron wont be ideal but this is the most common type
    # The distance d (nm) establishes the distance from the edge of the protein to the edge of the simulated box.
    build_box = 'gmx editconf -f clean_protein.gro -o clean_protein_box.gro -c -d 1.0 -bt dodecahedron'
    subprocess.run(build_box.split(), shell=True, check=True)

    # spc216.gro is included in gromacs and used in water tip3p solutions. This is the step that needs to change for non-water solution
    create_solution = 'gmx solvate -cp clean_protein_box.gro -cs spc216.gro -o protein_solv.gro -p topol.top'
    subprocess.run(create_solution.split(), shell=True, check=True)

    # Now we add a more natural salinity to our water solution (0.15M NaCl)
    # You need to create an empty .mdp file to run the next function. This can be added to to alter molecular dynamics equations
    make_empty_file = 'touch ions.mdp'
    subprocess.run(make_empty_file.split(), shell=True, check=True)

    add_ions = 'gmx grompp -f ions.mdp -c protein_solv.gro -p topol.top -o ions.tpr'
    subprocess.run(add_ions.split(), shell=True, check=True)

    # The print statement is one of the user inputs mentioned above to interact with the program
    #TODO: print directly in python
    print_inputs = 'printf "SOL\n" | gmx genion -s ions.tpr -o protein_solv.gro -conc 0.15 -p topol.top -pname NA -nname CL -neutral'
    subprocess.run(print_inputs.split(), shell=True, check=True)