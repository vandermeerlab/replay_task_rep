import torch.nn as nn
import torch.nn.functional as F

class MLP(nn.Module):
  """
    A simple feedforward MLP with a single hidden layer
        of relu nonlinearities. trained with SGD on MSE loss.
  """

  def __init__(self, in_dim, hid_dim, out_dim, scale_w1, scale_w2):
    """
    Initialize LNNet parameters

    Args:
      in_dim: int
        Input dimension
      out_dim: int
        Ouput dimension
      hid_dim: int
        Hidden dimension
      scale_w1: float
        Scale for input-to-hidden weights
      scale_w2: float
        Scale for hidden-output weights

    Returns:
      Nothing
    """
    super().__init__()
    self.in_hid = nn.Linear(in_dim, hid_dim, bias=False)
    nn.init.normal_(self.in_hid.weight, mean=0.0, std=scale_w1)

    self.hid_out = nn.Linear(hid_dim, out_dim, bias=False)
    nn.init.normal_(self.hid_out.weight, mean=0.0, std=scale_w2)

  def forward(self, x):
    """
    Forward pass of LNNet

    Args:
      x: torch.Tensor
        Input tensor

    Returns:
      hid: torch.Tensor
        Hidden layer activity
      out: torch.Tensor
        Output/Prediction
    """
    hid = F.relu(self.in_hid(x))  # Hidden activity
    out = self.hid_out(hid)  # Output (prediction)
    return out, hid
