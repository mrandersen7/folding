
import subprocess


class Protein:
    
    @property
    def name(self):
        return self.protein_pdb.split('.')[0]

    def __init__(self, path=None, ff='default', box='dodecahedron'):

        # can either be local file path or a url to download
        if path is None:
            path = self.select_random_protein()

        self.protein_pdb = path
        self.ff = ff
        self.box = box
        # I don't know what this should look like
        self.energy_min = ['1','2','3','4']
        self.remaining_steps = []

    def __str__(self):
        return f"Protein({self.protein_pdb}, {self.ff}, {self.box})"


    def __repr__(self):
        return self.__str__()

    def select_random_protein(self):
        """This function is really important as its where you select the protein you want to fold
        """
        return 'test_protein.pdb'

    def create_environment(self):
        """This function creates the environment the protein is folding in
        """

    def score(self):
        """This function scores the protein configuration
        """

    def preprocess(self):

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


    def run_first_step(self):
        # This is another input that has increible importance on the protein folding. We are going to run energy minimization
        # by first consoldating all our parameters and then running the minimization itself. The default works well in standard
        # cases but there are a lot of knobs to turn

        setup_first = 'gmx grompp -f energy_min_1 -c protein_solv.gro -p topol.top -o em.tpr'
        subprocess.run(setup_first.split(), shell=True, check=True)

        # -v can be added as a term to make this step print out progress. It normally takes a little while
        solve_first = 'gmx mdrun -deffnm em'
        subprocess.run(solve_first.split(), shell=True, check=True)
        
        self.remaining_steps.remove('1')