# Improved MCQD

We improve on maximum clique dynamic algorithm (MCQD) introduced by Konc 2007.
We use GNN to predict the best `Tlimit` parameter value.
The value of `Tlimit` is empiricaly estimated in the original paper and is set to `0.0025`.
The `Tlimit` parameter is used to decide whether or not to perform vertex sorting.

**Code has been tested on OSx and Linux. You might need to install additional libraries.**

### Run on Linux and OSx:
```
brew update
brew install boost openssl

conda create imcqd
conda activate imcqd
conda install -n imcqd spektral
conda install -n imcqd sklearn numpy pandas
```

```
# run this command first to compile the code
./compile.sh

# to generate training and testing data run (might take some time):
./build/imcqd

# to run just python script
cd python
python -W ignore main.py
```

# References
[1] Konc, Janez, and Dušanka Janezic. "An improved branch and bound algorithm for the maximum clique problem." proteins 4.5 (2007).

[2] Zhou, Jie, et al. "Graph neural networks: A review of methods and applications." arXiv preprint arXiv:1812.08434 (2018).