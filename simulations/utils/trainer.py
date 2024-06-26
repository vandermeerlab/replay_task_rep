import numpy as np
import torch

import torch.nn as nn
import torch.nn.functional as F
import torch.optim as optim

def train_MLP(model, inputs, targets, n_epochs, lr, wd_lambda=0, noise_lev=0):
  """
  Training function

  Args:
    model: torch nn.Module
      The neural network
    inputs: torch.Tensor
      Features (input) with shape `[batch_size, input_dim]`
    targets: torch.Tensor
      Targets (labels) with shape `[batch_size, output_dim]`
    n_epochs: int
      Number of training epochs (iterations)
    lr: float
      Learning rate
    wd_lambda: float
      Weight-decaying parameters
    noise_lev: float
      Noise level added to inputs

  Returns:
    losses: np.ndarray
      Record (evolution) of training loss
  """
  # Loss records
  losses = np.zeros(n_epochs)
  # Relative changes of weights in input-hidden layer and hidden-output layer
  weight_rel_changes = [np.zeros(n_epochs), np.zeros(n_epochs)]

  optimizer = optim.SGD(model.parameters(), lr=lr, weight_decay=wd_lambda)
  criterion = nn.MSELoss()

  init_weight_norms = []
  for weight in model.parameters():
    init_weight_norms.append(torch.linalg.norm(weight))

  for i in range(n_epochs):
    optimizer.zero_grad()
    noisy_inputs = inputs + torch.randn(inputs.shape) * noise_lev
    predictions, _ = model(noisy_inputs)
    loss = criterion(predictions, targets)
    loss.backward()
    optimizer.step()

    # task_predictions, _ = model(X)
    # task_predictions = task_predictions.detach()
    # task_loss = criterion(task_predictions, y)
    # task_losses[i] = task_loss.item()

    # Logging (recordings)
    losses[i] = loss.item()
    # Relative changes in weights
    for w_i, weight in enumerate(model.parameters()):
      weight_rel_changes[w_i][i] = (torch.linalg.norm(weight) - init_weight_norms[w_i]) / init_weight_norms[w_i]

  return losses, weight_rel_changes
