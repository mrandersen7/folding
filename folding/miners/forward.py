



# This is the core miner function, which decides the miner's response to a valid, high-priority request.
def forward(synapse: template.protocol.Dummy) -> template.protocol.Dummy:
    # TODO(developer): Define how miners should process requests.
    # This function runs after the synapse has been deserialized (i.e. after synapse.data is available).
    # This function runs after the blacklist and priority functions have been called.
    # Below: simple template logic: return the input value multiplied by 2.
    # If you change this, your miner will lose emission in the network incentive landscape.
    synapse.dummy_output = synapse.dummy_input * 2
    return synapse