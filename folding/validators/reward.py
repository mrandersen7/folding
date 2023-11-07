
import torch
import bittensor as bt
from typing import List

from folding.validators.protein import Protein




def reward(protein: Protein, response: bt.Synapse) -> float:

    return protein.energy(response)


def get_rewards(
    protein: Protein, responses: List[bt.Synapse]
) -> torch.FloatTensor:
    """Applies the reward model across each call. Unsuccessful responses are zeroed."""
    # Get indices of correctly responding calls.

    successful_response_indices: List[int] = [
        idx
        for idx, resp in enumerate(responses)
        if resp.dendrite.status_code == 200
    ]

    # Get all responses from responding calls.
    successful_responses: List[str] = [
        responses[idx].md_output for idx in successful_response_indices
    ]

    # Reward each response.
    successful_rewards = [reward(protein, response) for response in successful_responses]

    # Softmax rewards across samples.
    successful_rewards_normalized = torch.softmax(torch.tensor(successful_rewards), dim=0)

    # Init zero rewards for all calls.
    filled_rewards = torch.ones(len(responses), dtype=torch.float32) * torch.nan
    filled_rewards_normalized = torch.zeros(len(responses), dtype=torch.float32)

    # Fill reward tensor.
    for idx, reward, reward_normalized in zip(
        successful_response_indices,
        successful_rewards,
        successful_rewards_normalized,
    ):
        filled_rewards[idx] = reward
        filled_rewards_normalized[idx] = reward_normalized

    # Return the filled rewards.
    return filled_rewards, filled_rewards_normalized