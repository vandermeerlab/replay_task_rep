import numpy as np
import matplotlib.pyplot as plt
from .constants import COLORS

def plot_var(var_array,
            param_labels=['rich', 'lazy'],
            colors = [COLORS['rich'], COLORS['lazy']],
            y_ticks=None,
            y_label="error",
            fig_size=(10, 5)):
  """
  Plot time-varying variables function

  Args:
    var_array: np.ndarray
      Log of variables (ex: MSE loss) per epoch
    param_labels: list
      List of parameter labels

  Returns:
    Nothing
  """
  fig, axis = plt.subplots()
  fig.set_size_inches(fig_size)
  for p_i, param in enumerate(param_labels):
    axis.plot(np.mean(var_array[p_i, :, :], 0), color=colors[p_i], label=param)

  axis.set_xlabel("epoch")
  axis.set_ylabel(y_label)
  axis.set_xticks([0, 2500, 5000])
  if y_ticks:
      axis.set_yticks(y_ticks)
  axis.legend()
  return fig, axis

def plot_task_rep(X, y, colors=[], marker_size=200):
  fig, axis = plt.subplots()
  for i in range(y.shape[0]):
      if colors:
         color = colors[i]
      else:
        if X[i, 0] == -1. and X[i, 1] == 1.:
           color = COLORS['left_fr']
        elif X[i, 0] == 1. and X[i, 1] == 1.:
           color = COLORS['right_fr']
        elif X[i, 0] == -1. and X[i, 1] == -1.:
           color = COLORS['left_wr']
        elif X[i, 0] == 1. and X[i, 1] == -1.:
           color = COLORS['right_wr']
      marker = '^' if y[i] == 1 else 'v'
      axis.scatter(X[i, 0], X[i, 1], c = color, marker=marker, s=marker_size)

  axis.set_xlim(-1.5, 1.5)
  axis.set_xticks([-1, 1])
  axis.set_xticklabels(['left', 'right'])
  axis.set_ylim(-1.5, 1.5)
  axis.set_yticks([-1, 1])
  axis.set_yticklabels(['wr', 'fr'])

  return fig, axis

def plot_hidden_in(X, axis):
  axis.scatter(X[:, 0], X[:, 1], color=COLORS['hidden_in'])
  axis.set_xlabel('input weights $\it{i}$')
  axis.set_ylabel('input weights $\it{j}$')
  axis.set_xlim(-0.5, 0.5)
  axis.set_ylim(-0.5, 0.5)
  axis.set_xticks([-0.4, 0, 0.4])
  axis.set_yticks([-0.4, 0, 0.4])

  return axis

def plot_hidden_out(X, axis):
  axis.scatter(np.arange(X.shape[1]), X[0], color=COLORS['hidden_out'])
  axis.set_xlabel('unit')
  axis.set_ylabel('output weights')
  axis.set_ylim(-0.2, 0.2)
  axis.set_yticks([-0.2, 0, 0.2])

  return axis
