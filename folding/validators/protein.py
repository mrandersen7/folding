
import os
import re
import tqdm
import requests
import gromacs

import bittensor as bt

from dataclasses import dataclass


# root level directory for the project (I HATE THIS)
ROOT_DIR = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

@dataclass
class Protein:

    @property
    def name(self):
        return self.protein_pdb.split('.')[0]

    def __init__(self, pdb_id=None, ff='charmm27', box='dodecahedron', max_steps=None):

        # can either be local file path or a url to download
        if pdb_id is None:
            pdb_id = self.select_random_pdb_id()

        self.pdb_id = pdb_id

        self.output_directory = os.path.join(ROOT_DIR,'data', self.pdb_id)
        # if directory doesn't exist, download the pdb file and save it to the directory
        if not os.path.exists(self.output_directory):
            os.makedirs(self.output_directory)
            self.download_pdb()
        else:
            bt.logging.info(f'PDB file {self.pdb_id}.pdb already exists in path {self.output_directory!r}.')

        self.ff = ff
        self.box = box

        self.gro_path = os.path.join(self.output_directory, 'em.gro')
        self.topol_path = os.path.join(self.output_directory, 'topol.top')

        if not os.path.exists(self.gro_path) or not os.path.exists(self.topol_path):
            self.generate_input_files()


        required_files = ['em.gro','topol.top','posre.itp']
        self.md_inputs = {}
        for file in required_files:
            self.md_inputs[file] = open(os.path.join(self.output_directory, file), 'r').read()

        mdp_files = ['nvt.mdp','npt.mdp','md.mdp']
        for file in mdp_files:
            content = open(os.path.join(self.output_directory, file), 'r').read()
            if max_steps is not None:
                content = re.sub('nsteps\\s+=\\s+\\d+',f'nsteps = {max_steps}',content)
            self.md_inputs[file] = content


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
            bt.logging.info(f'PDB file {self.pdb_id}.pdb downloaded successfully from {url} to path {path!r}.')
        else:
            bt.logging.error(f'Failed to download PDB file with ID {self.pdb_id} from {url}')
            raise Exception(f'Failed to download PDB file with ID {self.pdb_id}.')

    # Function to generate GROMACS input files
    def generate_input_files(self):
        # Change to output directory
        os.chdir(self.output_directory)

        # Commands to generate GROMACS input files
        gromacs.pdb2gmx(f={self.pdb_id}.pdb, o="processed.gro", p="topol.top", ff={self.ff} , water="tip3p") # Input the file into GROMACS and get three output files: topology, position restraint, and a post-processed structure file
        gromacs.editconf(f="processed.gro", o="newbox.gro", bt={self.box2}, d=1.0, princ=True, input="Protein")
        gromacs.solvate(cp="newbox.gro", cs="spc216.gro", p="topol.top", o="solvated.gro")
        open("ions.mdp", 'w')
        gromacs.grompp(f="ions.mdp", c="solvated.gro", p="topol.top", o="ions.tpr")
        
        #TODO Make this more pythonic. The input in editconf works, so why not here?
        #gromacs.genion(s="ions.tpr", o="protein_solv.gro", conc = 0.15, p="topol.top", pname="NA", nname="CL", neutral=True, input="13")
        os.system('echo "13" | gmx genion -s ions.tpr -o solv_ions.gro -p topol.top -pname NA -nname CL -neutral')

        # Run the first step of the simulation
        gromacs.grompp(f="../minim.mdp", c="solv_ions.gro", p="topol.top", o="em.tpr")
        gromacs.mdrun(v=True, deffnm="em") # Run energy minimization
        

        # strip away trailing number in forcefield name e.g charmm27 -> charmm
        ff_base = ''.join([c for c in self.ff if not c.isdigit()])
        # Copy mdp template files to output directory
        commands += [
            f'cp ../nvt-{ff_base}.mdp nvt.mdp',
            f'cp ../npt-{ff_base}.mdp npt.mdp',
            f'cp ../md-{ff_base}.mdp  md.mdp '
        ]

        for cmd in tqdm.tqdm(commands):
            os.system(cmd)

        # We want to catch any errors that occur in the above steps and then return the error to the user
        return True

    def gro_id(self, path):
        with open(path, 'r') as f:
            name, length, *lines, _ = f.readlines()
            length = int(length)
            print(name, length, len(lines))
        gro_id = hash(name+''.join([''.join(cols[:3] if len(cols)==6 else cols[:2]) for cols in lines[:length]]))
        return gro_id

    def reward(self, md_output: dict, mode: str='13'):
        """Calculates the free energy of the protein folding simulation
        # TODO: Each miner files should be saved in a unique directory and possibly deleted after the reward is calculated
        """

        filetypes = {}
        for filename, content in md_output.items():
            filetypes[filename.split('.')[-1]] = filename
            # loop over all of the output files and save to local disk
            with open(os.path.join(self.output_directory, filename), 'wb') as f:
                f.write(content)

        edr = filetypes.get('edr')
        if not edr:
            bt.logging.error(f'No .edr file found in md_output ({list(md_output.keys())}), so reward is zero!')
            return 0

        gro = filetypes.get('gro')
        if not gro:
            bt.logging.error(f'No .gro file found in md_output ({list(md_output.keys())}), so reward is zero!')
            return 0
        if self.gro_hash(self.gro_path) != self.gro_hash(gro):
            bt.logging.error(f'The hash for .gro file is incorrect, so reward is zero!')
            return 0

        #TODO Make this more pythonic, finish validation mechanism
        commands = [
            f'echo "13"  | gmx energy -f {edr} -o free_energy.xvg'
        ]

        # TODO: we still need to check that the following commands are run successfully
        for cmd in tqdm.tqdm(commands):
            os.system(cmd)

        energy_path = os.path.join(self.output_directory, 'free_energy.xvg')
        free_energy = self.get_average_free_energy(energy_path)

        # return the negative of the free energy so that larger is better
        return -free_energy

    # Function to read the .xvg file and compute the average free energy
    def get_average_free_energy(self, filename):
        # Read the file, skip the header lines that start with '@' and '&'
        with open(filename) as f:
            last_line = f.readlines()[-1]

        # The energy values are typically in the second column
        last_energy = last_line.split()[-1]

        return float(last_energy)