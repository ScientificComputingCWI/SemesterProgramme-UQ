# Install instructions

### Cloning this repository

1) Install [Git](https://git-scm.com/downloads). You can check if it is already installed by typing git in a terminal.
2) Type the following

   ```
    git clone https://github.com/ScientificComputingCWI/SemesterProgramme-UQ.git
   ```

### Installing EasyVVUQ

EasyVVUQ depends upon the chaospy library, which is currently undergoing work to be up to date with Numpy 2. To be on the safe side, we therefore recommend to install EasyVVUQ in a virtual environment where we enforce Numpy 1. That said, if you currently have a Numpy version < 2 (check with `np.__version__`) you can simply install EasyVVUQ with `pip install EasyVVUQ`. If you do have Numpy 2, or if you prefer a separate virtual environment voor the autumn school, follow the steps below.

1) Choose a directory where you wish to install the virtual environment, e.g. the home directory `~/`
2) Type the following to create a virtual environment:

```
cd  ~/
python3 -m venv myenv
```

3) activate the environment (on Linux):

```
source myenv/bin/activate
```

3) activate the environment (on Windows), run the `activate` file in:

```
myenv/Scripts/activate
```

4) Install EasyVVUQ, Numpy (< 2.0), and jupyterlab in `venv` (Linux)

```
pip install easyvvuq
pip install "numpy>=1.0,<2.0"
pip install jupyter-lab
```

4) Install EasyVVUQ, Numpy (< 2.0), and jupyterlab in `venv` (Windows)

```
pip install easyvvuq
pip install "numpy>=1.0,<2.0"
pip install jupyter lab
```

5) Run jupyterlab (Linux):

```
myenv/bin/jupyterlab
```

5) Run jupyterlab (Windows):

```
jupyter lab
```

6) Run the `forward_UQ/exercises/EasyVVUQ_install_check.ipynb` notebook.

### Installing FabSim3 with FabUQCampaign

For detailed FabSim3 installation instructions, click [here](https://fabsim3.readthedocs.io/en/latest/installation/). Briefly, the installation entails the following steps, assuming you will install FabSim3 in your home directory:

```
pip3 install ruamel.yaml rich
cd
git clone https://github.com/djgroen/FabSim3.git
cd FabSim3
python3 configure_fabsim.py
```

At this point you will see some instructions on adding FabSim3 to your environment variables.

Finally, you can check if the install was successfull using the `FabDummy` plugin:

```
fabsim localhost install_plugin:FabDummy
fabsim localhost dummy:dummy_test
```

which should execute without error. Similarly, you can install the `FabUQCampaign` [plugin](https://github.com/wedeling/FabUQCampaign) via

```
fabsim localhost install_plugin:FabUQCampaign
```

This will install the plugin in `FabSim3/plugins`.
