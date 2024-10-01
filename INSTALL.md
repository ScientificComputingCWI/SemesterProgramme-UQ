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
python -m venv myenv
```

3) activate the environment (on Linux):

```
source myenv/bin/activate
```

3) activate the environment (on Windows), run the `activate.bat` file in:

```
myenv/Scripts/activate.bat
```

4) Install EasyVVUQ, Numpy (< 2.0), and jupyterlab in `venv`

```
pip install easyvvuq
pip install "numpy>=1.0,<2.0"
pip install jupyterlab
```

5) Run jupyterlab:

```
myenv/bin/jupyterlab
```

6) Run the `forward_UQ/exercises/EasyVVUQ_install_check.ipynb` notebook.
