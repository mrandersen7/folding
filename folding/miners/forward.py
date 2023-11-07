
import gmxapi as gmx
import folding.protocol

# This is the core miner function, which decides the miner's response to a valid, high-priority request.
def forward(synapse: folding.protocol.Synapse) -> folding.protocol.Synapse:
    # TODO(developer): Define how miners should process requests.
    # This function runs after the synapse has been deserialized (i.e. after synapse.data is available).
    # This function runs after the blacklist and priority functions have been called.

    protein = synapse.protein
    simulation_input = gmx.read_tpr(protein.tpr_filename, label=protein.label)
    md = gmx.mdrun(simulation_input, runtime_args=protein.runtime_args)

    md.run()

    # Here we attach the simulation output to the synapse, to be sent back to the validator.
    synapse.md_output = md.output

    return synapse