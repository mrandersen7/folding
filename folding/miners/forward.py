
# import gmxapi as gmx
import os
import tqdm
import glob
import folding.protocol

# This is the core miner function, which decides the miner's response to a valid, high-priority request.
def forward(synapse: folding.protocol.Synapse) -> folding.protocol.Synapse:
    # TODO(developer): Define how miners should process requests.
    # This function runs after the synapse has been deserialized (i.e. after synapse.data is available).
    # This function runs after the blacklist and priority functions have been called.

    protein = synapse.protein
    mdrun_args = synapse.mdrun_args # can be something like '-maxh 0.1'

    synapse.md_output = {}

    # Make sure the output directory exists and if not, create it
    if not os.path.exists(protein.output_directory):
        os.makedirs(protein.output_directory)
    # Change to output directory
    os.chdir(protein.output_directory)

    # Write the topology and coordinate files to the output directory
    with open('em.gro', 'w') as file:
        file.write(protein.gro)

    with open('topol.top', 'w') as file:
        file.write(protein.topol)

    # Commands to run GROMACS simulations
    # TODO: investigate how to run GROMACS simulations in parallel/ on GPUs
    # TODO: investigate the use of python API instead of command line
    commands = [
            'gmx grompp -f nvt.mdp -c em.gro -r em.gro -p topol.top -o nvt.tpr', # Temperature equilibration
            'gmx mdrun -deffnm nvt' +mdrun_args,
            'gmx grompp -f npt.mdp -c nvt.gro -t nvt.cpt -p topol.top -o npt.tpr', # Pressure equilibration
            'gmx mdrun -deffnm npt' + mdrun_args,
            'gmx grompp -f md.mdp -c npt.gro -t npt.cpt -p topol.top -o md_0_1.tpr', # Production run
            'gmx mdrun -deffnm md_0_1' + mdrun_args
    ]

    for cmd in tqdm.tqdm(commands):
        # We want to catch any errors that occur in the above steps and then return the error to the user
        os.system(cmd)

    # load the output files as bytes and add to synapse.md_output
    for filename in glob.glob('md_0_1.*'):
        # file_ext = file.split('.')[-1]
        with open(file, 'rb') as f:
            synapse.md_output[filename] = f.read()

    return synapse


