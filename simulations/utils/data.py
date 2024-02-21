import numpy as np
import torch

def create_xor_dataset(x1, x2, y, n_repeats=50, noise_lev = 0.1):
  # Repeat the numbers from x1, x2, and y n_repeats times
  x1 = np.repeat(x1, n_repeats)
  x2 = np.repeat(x2, n_repeats)
  y =  np.repeat(y,  n_repeats)

  # Add noise to data
  x1 = x1 + np.random.randn(x1.shape[0]) * noise_lev
  x2 = x2 + np.random.randn(x2.shape[0]) * noise_lev

  # Shuffle
  index_shuffle = np.arange(x1.shape[0])
  np.random.shuffle(index_shuffle)

  x1 = x1.astype(np.float32)
  x2 = x2.astype(np.float32)
  y  = y.astype(np.float32)

  x1 = x1[index_shuffle]
  x2 = x2[index_shuffle]
  y  = y [index_shuffle]

  # Convert data to tensors
  x1_torch = torch.from_numpy(x1).clone().view(-1, 1)
  x2_torch = torch.from_numpy(x2).clone().view(-1, 1)
  y_torch = torch.from_numpy(y).clone().view(-1, 1)

  # Combine X1 and X2
  X_torch = torch.hstack([x1_torch, x2_torch])

  return X_torch, y_torch
