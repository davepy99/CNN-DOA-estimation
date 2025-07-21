Author: David Pauly

CNN project to improve the speed of Direction-of-Arrival (DOA) algorithms for linear phase arrays through a CNN, MUSIC/MVDR hybrid approach. The project consists of 3 main parts:

- A CNN (DOA_Estimation_CNN.ipynb) that is trained on a generated dataset with covariance matrices labeled with the corresponding DOAs
- A Matlab Skripts that Generates the datasets for training and testing (/Dataset Generation)
- A Matlab skript that takes the predicted ranges from the CNN and resolves the precise DOA within these ranges to reach full MUSIC/MVDR precision (/DOA Algorithms)

