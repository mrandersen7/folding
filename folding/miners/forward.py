
# import gmxapi as gmx
import os
import tqdm
import glob
import bittensor as bt
import folding.protocol

# This is the core miner function, which decides the miner's response to a valid, high-priority request.
def forward(synapse: folding.protocol.Synapse) -> folding.protocol.Synapse:
    # This function runs after the synapse has been deserialized (i.e. after synapse.data is available).
    # This function runs after the blacklist and priority functions have been called.
    # TODO: Determine how many steps to run based on timeout

    bt.logging.info(f"Running GROMACS simulation for protein: {synapse.pdb_id} with files {synapse.md_inputs.keys()} mdrun_args: {synapse.mdrun_args}")
    synapse.md_output = {}

    output_directory = os.path.join('./data/miner', synapse.pdb_id)
    # Make sure the output directory exists and if not, create it
    if not os.path.exists(output_directory):
        os.makedirs(output_directory)
    # Change to output directory
    os.chdir(output_directory)

    # The following files are required for GROMACS simulations and are recieved from the validator
    for filename, content in synapse.md_inputs.items():
        # Write the file to the output directory
        with open(filename, 'w') as file:
            file.write(content)

    # Commands to run GROMACS simulations
    # TODO: investigate how to run GROMACS simulations in parallel/ on GPUs
    # TODO: rework how we want md_run_args to work/what default values are
    # TODO: bt.logging.info was incorporated before, integrate it here
    
    gromacs.grompp(f="nvt.mdp", c="em.gro", r ="em.gro", p="topol.top", o="nvt.tpr") # Temperature equilibration
    gromacs.mdrun(v=True, deffnm="nvt") #, maxh =1)

    gromacs.grompp(f="npt.mdp", c="nvt.gro", r ="nvt.gro", t='nvt.cpt' p="topol.top", o="npt.tpr") # Pressure equilibration
    gromacs.mdrun(v=True, deffnm="npt") #, maxh =1)

    gromacs.grompp(f="md.mdp", c="npt.gro", r ="npt.gro", p="topol.top", o="md_0_1.tpr") # Production run
    gromacs.mdrun(v=True, deffnm="md_0_1") #, maxh =1)
    

    # load the output files as bytes and add to synapse.md_output
    for filename in glob.glob('md_0_1.*'):
        bt.logging.info(f"Attaching file: {filename} to synapse.md_output")
        with open(filename, 'rb') as f:
            synapse.md_output[filename] = f.read()

    return synapse


