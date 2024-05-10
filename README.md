# replay_task_rep

This repo contains research code for the paper: [Paradoxical replay can protect contextual task representations from destructive interference when experience is unbalanced](https://www.biorxiv.org/content/10.1101/2024.05.09.593332v1).

`/simulations` folder contains code for neural networks simulation (`Fig.1 - 3`).

`/analysis` folder contains code for neural data analysis (`Fig.4`). Raw data can be found at [here](https://datasets.datalad.org/?dir=/workshops/mind-2017/MotivationalT).

## Getting Started

Either download or clone the repo:
```
git clone https://github.com/vandermeerlab/replay_task_rep
```
Then navigate to the downloaded folder:
```
cd /path/to/replay_task_rep
```
Install the package and requirements:
```
pip install -r requirements.txt
```

## Requirements
```
# Neural network simulations
Python 3.10

## Packages
pytorch == 2.1.2
sklearn == 1.3.1
numpy == 1.23.5
scipy == 1.9.3
matplotlib == 3.7.0

# Neural data analysis
MATLAB 2022b (Mathworks)
```
