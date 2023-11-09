
import os
import tqdm
import subprocess
import requests
# import gmxapi as gmx


class Protein:

    @property
    def name(self):
        return self.protein_pdb.split('.')[0]

    def __init__(self, pdb_id=None, ff='CHARMM27', box='dodecahedron'):

        # can either be local file path or a url to download
        if pdb_id is None:
            pdb_id = self.select_random_pdb_id()

        self.pdb_id = pdb_id

        self.output_directory = os.path.join('./data', self.pdb_id)
        # if directory doesn't exist, download the pdb file and save it to the directory
        if not os.path.exists(self.output_directory):
            os.makedirs(self.output_directory)
            self.download_pdb()

        self.ff = ff
        self.box = box
        # I don't know what this should look like
        self.energy_min = ['1','2','3','4']

        self.remaining_steps = []

    def __str__(self):
        return f"Protein({self.protein_pdb}, {self.ff}, {self.box})"


    def __repr__(self):
        return self.__str__()

    def select_random_pdb_id(self):
        """This function is really important as its where you select the protein you want to fold
        """
        return '1UBQ'

    # Function to download PDB file
    def download_pdb(self):
        url = f'https://files.rcsb.org/download/{self.pdb_id}.pdb'
        r = requests.get(url)
        if r.status_code == 200:
            with open(os.path.join(self.output_directory, f'{self.pdb_id}.pdb'), 'w') as file:
                file.write(r.text)
            print(f'PDB file {self.pdb_id}.pdb downloaded successfully.')
        else:
            print(f'Failed to download PDB file with ID {self.pdb_id}.')

    # Function to generate GROMACS input files
    def generate_gromacs_input_files(self, run_first_step=True):
        # Change to output directory
        os.chdir(self.output_directory)

        # Commands to generate GROMACS input files
        commands = [
            f'gmx pdb2gmx -f {self.pdb_id}.pdb -o processed.gro -water spce', # Input the file into GROMACS and get three output files: topology, position restraint, and a post-processed structure file
            'gmx editconf -f processed.gro -o newbox.gro -c -d 1.0 -bt cubic', # Build the "box" to run our simulation of one protein molecule
            'gmx solvate -cp newbox.gro -cs spc216.gro -o solvated.gro -p topol.top',
            'gmx grompp -f ions.mdp -c solvated.gro -p topol.top -o ions.tpr',
            'echo "13" | gmx genion -s ions.tpr -o solv_ions.gro -p topol.top -pname NA -nname CL -neutral',
        ]
        if run_first_step:
            commands += [
                'gmx grompp -f minim.mdp -c solv_ions.gro -p topol.top -o em.tpr',
                'gmx mdrun -v -deffnm em' # Run energy minimization
            ]

        for cmd in tqdm.tqdm(commands):
            os.system(cmd)
            
        # We want to catch any errors that occur in the above steps and then return the error to the user
        return True
    

    def energy(self, md_output):
        """This is potentailly where we calculate the energy of an md simulation result
        """
        return 42
