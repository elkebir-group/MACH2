# MACH2

## Installation

These packages are necessary to run MACH2.

        graphviz
        gurobi>=9.0.0,<10.0.0
        numpy==1.24.3
        pandas==2.0.1
        pyomo==6.5.0
        networkx==3.1

Can be installed using conda

        $ conda create -n mach2 numpy pandas pyomo networkx gurobi=9 python-graphviz jupyterlab python=3
        $ conda activate mach2

The package can be loaded to Jupyter Lab. Since `setup.py` is yet to set up properly, one needs to run Jupyter Notebook from package folder. 