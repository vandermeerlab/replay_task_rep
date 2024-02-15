import numpy as np
import matplotlib.pyplot as plt

def plot_var(var_array,
            parameters,
            y_label="MSE",
            title="Training loss (Mean Squared Error)"):
  """
  Plot time-varying variables function

  Args:
    title: string
      Specifies plot title
    var_array: np.ndarray
      Log of variables (ex: MSE loss) per epoch

  Returns:
    Nothing
  """
  plt.figure(figsize=(10, 5))
  for p_i, param in enumerate(parameters):
    plt.plot(np.mean(var_array[p_i, :, :], 0), label=param)
    plt.xlabel("Epoch")
    plt.ylabel(y_label)
    plt.title(title)
  plt.legend()
  plt.show()
