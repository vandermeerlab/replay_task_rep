import numpy as np
import matplotlib.pyplot as plt

def plot_var(var_array,
            param_labels=['rich', 'lazy'],
            y_ticks=None,
            y_label="MSE",
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
  fig, ax = plt.subplots()
  fig.set_size_inches(fig_size)
  for p_i, param in enumerate(param_labels):
    ax.plot(np.mean(var_array[p_i, :, :], 0), label=param)
    ax.set_xlabel("epoch")
    ax.set_ylabel(y_label)
    ax.set_xticks([0, 2500, 5000])
    if y_ticks:
        ax.set_yticks(y_ticks)
    ax.legend()
  return fig
