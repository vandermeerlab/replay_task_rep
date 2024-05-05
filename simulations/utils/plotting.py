import numpy as np
import matplotlib.pyplot as plt
from .constants import COLORS

def plot_var(var_array,
            running_avg = 0,
            epoch_start = 0,
            param_labels=['rich', 'lazy'],
            colors = [COLORS['rich'], COLORS['lazy']],
            y_ticks=None,
            y_tick_labels=[],
            y_label="error",
            fig_size=(8, 5)):
  """
  Plotting time-varying variables

  Args:
    var_array: np.ndarray
      Log of variables (ex: MSE loss) per epoch
    running_avg: int
      window size for running average calculation
    param_labels: list
      List of parameter labels
    other parameters:
      plt figure parameters.

  Returns:
    Nothing
  """
  fig, axis = plt.subplots()
  fig.set_size_inches(fig_size)
  for p_i, param in enumerate(param_labels):
    if running_avg > 0:
      var = np.convolve(np.mean(var_array[p_i, :, :], 0), np.ones(running_avg)/running_avg, mode='valid')
    else:
      var = np.mean(var_array[p_i, :, :], 0)
    axis.plot(np.arange(epoch_start, var.shape[0]+epoch_start), var,
              color=colors[p_i], linewidth=2, label=param)

  axis.set_xlabel("epoch")
  axis.set_ylabel(y_label)
  axis.set_xticks(np.array([0, 1250, 2500, 3750, 5000]) + epoch_start)
  axis.set_xticklabels(np.array([0+epoch_start, None, 2500+epoch_start, None, 5000+epoch_start]))
  if y_ticks:
      axis.set_yticks(y_ticks)
  if y_tick_labels:
    axis.set_yticklabels(y_tick_labels)
  axis.legend()
  return fig, axis

def plot_task_inputs(X, y, colors=[], markers=[], marker_size=200):
  """
  Plotting data inputs in the task representation.
  """
  fig, axis = plt.subplots()
  for i in range(y.shape[0]):
      if colors:
         color = colors[i]
      else:
        color = COLORS['pos'] if y[i] == 1 else COLORS['neg']
      if markers:
         marker = markers[i]
      else:
        marker = 'o' if X[i, 1] == 1 else '^'
      axis.scatter(X[i, 0], X[i, 1], c = color, marker=marker, s=marker_size)

  axis.set_xlim(-1.5, 1.5)
  axis.set_xticks([-1, 1])
  axis.set_xticklabels(['left', 'right'])
  axis.set_ylim(-1.5, 1.5)
  axis.set_yticks([-1, 1])
  axis.set_yticklabels(['ctx2', 'ctx1'])

  return fig, axis

def plot_hidden_in(X, axis, color):
  """
  Plotting input weights of the hidden layer units.
  """
  axis.scatter(X[:, 0], X[:, 1], color=color)
  axis.set_xlabel('input weights $W^{in}_{i}}$')
  axis.set_ylabel('input weights $W^{in}_{j}}$')
  axis.set_xlim(-0.5, 0.5)
  axis.set_ylim(-0.5, 0.5)
  axis.set_xticks([-0.4, -0.2, 0, 0.2, 0.4])
  axis.set_xticklabels([-0.4, None, 0, None, 0.4])
  axis.set_yticks([-0.4, -0.2, 0, 0.2, 0.4])
  axis.set_yticklabels([-0.4, None, 0, None, 0.4])

  return axis

def plot_hidden_out(X, axis, color):
  """
  Plotting output weights of the hidden layer units.
  """
  axis.scatter(np.arange(X.shape[1]), X[0], color=color)
  axis.set_xlabel('unit')
  axis.set_xlim(0, 500)
  axis.set_xticks([0, 125, 250, 375, 500])
  axis.set_xticklabels([0, None, 250, None, 500])
  axis.set_ylabel('output weights $W^{out}$')
  axis.set_ylim(-0.25, 0.25)
  axis.set_yticks([-0.2, -0.1, 0, 0.1, 0.2])
  axis.set_yticklabels([-0.2, None, 0, None, 0.2])

  return axis

def plot_3d_embedding(X, y, embedding):
  fig = plt.figure()
  fig.set_size_inches(8, 8)
  axis = fig.add_subplot(111, projection='3d')

  for i in range(X.shape[0]):
      color = COLORS['pos'] if y[i] == 1 else COLORS['neg']
      marker = 'o' if X[i, 1] == 1 else '^'
      axis.scatter(embedding[i, 0], embedding[i, 1], embedding[i, 2], c=color, marker=marker, s=200)

  axis.set_xlabel('PC1')
  axis.set_ylabel('PC2')
  axis.set_zlabel('PC3')
  axis.set_xticks([-15, 0, 15])
  axis.set_yticks([-15, 0, 15])
  axis.set_zticks([-15, 0, 15])

  return fig, axis
