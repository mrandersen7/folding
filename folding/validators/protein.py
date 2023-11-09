
import os
import glob
import tqdm
import requests
import pandas as pd
# import gmxapi as gmx

# root level directory for the project (I HATE THIS)
ROOT_DIR = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

class Protein:

    @property
    def name(self):
        return self.protein_pdb.split('.')[0]

    def __init__(self, pdb_id=None, ff='charmm27', box='dodecahedron'):

        # can either be local file path or a url to download
        if pdb_id is None:
            pdb_id = self.select_random_pdb_id()

        self.pdb_id = pdb_id

        self.output_directory = os.path.join(ROOT_DIR,'data', self.pdb_id)
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
        return f"Protein(pdb_id={self.pdb_id}, ff={self.ff}, box={self.box}, output_directory={self.output_directory})"


    def __repr__(self):
        return self.__str__()

    def select_random_pdb_id(self):
        """This function is really important as its where you select the protein you want to fold
        """
        return '1UBQ'

    # Function to download PDB file
    def download_pdb(self):
        url = f'https://files.rcsb.org/download/{self.pdb_id}.pdb'
        path = os.path.join(self.output_directory, f'{self.pdb_id}.pdb')
        r = requests.get(url)
        if r.status_code == 200:
            with open(path, 'w') as file:
                file.write(r.text)
            print(f'PDB file {self.pdb_id}.pdb downloaded successfully to path {path!r}.')
        else:
            print(f'Failed to download PDB file with ID {self.pdb_id}.')

    # Function to generate GROMACS input files
    def generate_input_files(self):
        # Change to output directory
        os.chdir(self.output_directory)

        # Commands to generate GROMACS input files
        commands = [
            f'gmx pdb2gmx -f {self.pdb_id}.pdb -ff {self.ff} -o processed.gro -water spce', # Input the file into GROMACS and get three output files: topology, position restraint, and a post-processed structure file
            f'gmx editconf -f processed.gro -o newbox.gro -c -d 1.0 -bt {self.box}', # Build the "box" to run our simulation of one protein molecule
            'gmx solvate -cp newbox.gro -cs spc216.gro -o solvated.gro -p topol.top',
            'touch ions.mdp', # Create a file to add ions to the system
            'gmx grompp -f ions.mdp -c solvated.gro -p topol.top -o ions.tpr',
            'echo "13" | gmx genion -s ions.tpr -o solv_ions.gro -p topol.top -pname NA -nname CL -neutral',
        ]
        # Run the first step of the simulation
        commands += [
            'gmx grompp -f ../minim.mdp -c solv_ions.gro -p topol.top -o em.tpr',
            'gmx mdrun -v -deffnm em' # Run energy minimization
        ]

        for cmd in tqdm.tqdm(commands):
            os.system(cmd)

        # print(os.listdir('.'))
        # # change back to the parent directory
        # os.chdir('../..')
        # print(os.listdir('.'))
        # print(os.listdir(self.output_directory))

        print(glob.glob('em.*'))
        # read the output files as strings and save to self.gro and self.topol
        with open('em.gro', 'rb') as file:
            self.gro = file.read()
        with open('topol.top', 'rb') as file:
            self.topol = file.read()

        # We want to catch any errors that occur in the above steps and then return the error to the user
        return True


    def reward(self, md_output: dict):
        """Calculates the free energy of the protein folding simulation
        """

        edr_filename = None
        for filename, content in md_output.items():
            if filename.endswith('.edr'):
                edr_filename = filename
            # loop over all of the output files and save to local disk
            with open(os.path.join(self.output_directory, filename), 'wb') as f:
                f.write(content)

        if not edr_filename:
            print('No .edr file found in md_output, so reward is zero!!!!!')
            return 0

        commands = [
            f'gmx energy -f {edr_filename} -o free_energy.xvg'
        ]

        # TODO: we still need to check that the following commands are run successfully
        for cmd in tqdm.tqdm(commands):
            os.system(cmd)

        energy_path = os.path.join(self.output_directory, 'free_energy.xvg')
        free_energy = self.get_average_free_energy(energy_path)

        # return the negative of the free energy so that larger is better
        return -free_energy

    # Function to read the .xvg file and compute the average free energy
    def get_average_free_energy(filename):
        # Read the file, skip the header lines that start with '@' and '&'
        data = pd.read_csv(filename, sep='\s+', comment='@', header=None)

        # The energy values are typically in the second column
        energy_values = data[1]

        # Calculate the average free energy
        average_energy = energy_values.mean()
        return average_energy