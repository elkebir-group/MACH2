# MACH2

A mathematical framework for inferring migration histories of metastatic cancer from clone phylogeny and the location of extant clones.

## Table of contents

1. [Installation](#installation)  
        - [Prerequisite](#prerequisite)  
        - [Install using `pip`](#install-using-pip)
2. [Usage instruction](#usage-instruction)  
        - [From `JupyterLab`](#from-jupyterlab)  
        - [From Command Line](#from-command-line)  
                <!-- - [I/O formats](#i-o-formats)   
                - [Usage](#usage)   -->

## 1. Installation

### Prerequisites

- **Python** - `MACH2` requires Python 3.7 or newer.
- **ILP solver** - `MACH2` requires an ILP solver installed to solve **PMH-TR**. Currently `MACH2` only supports `Gurobi optimizer`, but we are going to add support for more ILP solvers in the future. `MACH2` requires a valid Gurobi installation and license key. The location of Gurobi should be present in `LD_LIBRARY_PATH` (linux) or `DYLD_LIBRARY_PATH` (macOS) the license key should be saved in the environment variable `GRB_LICENSE_KEY`.

<!-- ### Install using `pip`

`MACH2` can easily be installed using `pip`, the package installer for Python. Open a terminal or command prompt and run the following command:

                $ pip install mach2

If you want to use `MACH2` with `JupyterLab`, you'll need additional dependencies. To install these optional dependencies, you can run the following command:

                $ pip install mach2[jupyter] 

PyPI version is currently broken. We'll update it ASAP.-->

### Install using `conda`

`MACH2` can be installed using `conda`. We advise to create a new environment in `conda`. If creating a new environment, dependencies can be installed simultaneously.

                $ conda create -n mach2 python=3 pandas networkx gurobi jupyterlab graphviz pygraphviz -c conda-forge -c gurobi
                $ conda activate mach2

If using existing conda environment, the following command installs the dependencies.

                $ conda install -c conda-forge -c gurobi pandas networkx gurobi jupyterlab graphviz pygraphviz

Next, we install `MACH2`. To that end, we download `MACH2` repository from GitHub and install it.

                $ git clone https://github.com/elkebir-group/MACH2.git
                $ cd MACH2
                $ pip install . --no-deps

## Usage Instruction

### I/O formats

We describe various formats used by `MACH2`.

1. **Tree file** : The tree file contains a list of edges that define the structure of a tree. Each line in the file represents an edge, and the edges should be in the format: `node1 node2`. For example:

        1   2
        2   3
        2   4 
        3   5

2. **Tree file with timing/comigrations** : Tree file with timestamps. Edges with the same timestamp belong to the same comigration, and a timestamp with `-1` represents non-migration. Each line corresponds to an edge in the format: `node1 node2 timestamp`. For example:

        1   2   -1
        2   3   1
        2   4   1
        3   5   2

3. **Observed labeling file** : The observed labeling file contains zero or more location labels assigned to each node of the input clonal tree. Each line in the file corresponds to a node and the labels assigned to it in the format: `node label1 label2 ...`. If a node is not observed anywhere, it may be skipped. For example:

        1   A   B
        3   B
        4   A   C
        5   C

4. **Location labeling file** : The location labeling file contains the unique location label of origin assigned to each node. Each line in the file corresponds to a node and the location label of origin are in the format: `node label`. For example:

        1   A
        2   B
        3   C

5. **Node of origin file** : The node of origin file maps the nodes of the refined tree to the nodes of the input tree. Each line in the file corresponds to a vertex and the labels are in the format: `leaf label`. For example:

        1   A
        2   B
        3   C
        
Additionaly, `MACH2` can output files in Graphviz DOT format or JSON format.

### Usage

`MACH2` takes as input two files - 

1. **Tree file** : Tree file describing the input clone tree.
2. **Observed labeling file** : Labeling file describing the observed labeling of input clone tree.


`MACH2` Can be run using command line, or can be directly accessed from `JupyterLab`.

#### From JupyterLab

The following code snippet imports `MACH2`, runs it for input tree file `input.tree` and input observed labeling file `input.observed.labeling`, and saves the solutions to a variable `solutions`.


                import mach2
                tree = mach2.MultiLabeledTree.from_files('input.tree', 'input.observed.labeling')
                solutions = mach2.MACH2(tree, primary_location='primary', criteria_ordering='UMC').solve()
                print(len(solutions))
                solutions.summary()


`solutions` is a `SolutionSet` object that behaves as a set. `print(len(solutions))` prints the number of retrived solutions.
The last line draws the summary graph for `solutions.
It is possible to inspect individual solutions too.


                solution1 = [sol for sol in solutions][0]
                solution1.draw()
                solution1.migration_graph().draw()

The second line draws the tree with node labeling, and the third line draws the corresponding migration graph.
For more details, check the documentation for each function.


#### From Terminal


For each solution, `MACH2` can output three types of files.

1. **Tree file with timing/comigrations** : Refined tree file with timestamps/comigrations.
2. **Location labeling file** : Location labeling file describing the location labeling of the refined tree.
3. **Node of origin file** : Node of origin file mapping refined tree nodes to input tree nodes.

Additionaly `MACH2` can return JSON file encoding all the solutions. The JSON file can be directly passed to [MACH2-viz](https://github.com/elkebir-group/mach2-viz). The exact format of the JSON file is described [here](https://github.com/elkebir-group/mach2-viz/README.md).

`MACH2` also prints `<primary location> <number of migrations> <number of comigrations> Optimal <running time (in seconds)>` on console.

`MACH2` can be run using python.

                usage: mach2 [-h] [-c CRITERIA] [-s] [-p PRIMARY] [--colormap COLORMAP] [--log] [-o OUTPUT] [--max_solutions MAX_SOLUTIONS] [-t THREADS] [--viz]
                        clonal_tree observed_labeling

                MACH2

                positional arguments:
                clonal_tree           Input clonal tree
                observed_labeling     Input observed labeling

                optional arguments:
                -h, --help            show this help message and exit
                -c CRITERIA, --criteria CRITERIA
                                        Criteria ordering
                -s, --seeding_locations
                                        Prioritize solutions with the least number of seeding locations (default=False)
                -p PRIMARY, --primary PRIMARY
                                        Primary anatomical location
                --colormap COLORMAP   Color map file
                --log                 Outputs Gurobi logging (default=False)
                -o OUTPUT, --output OUTPUT
                                        Output folder (default=current folder)
                --max_solutions MAX_SOLUTIONS
                                        Maximum number of solutions retained (default=37888)
                -t THREADS, --threads THREADS
                                        Number of threads
                --viz, --open_in_viz  Open the locations on MACH2-viz (default=False) 

An example execution

        $ mach2 data/ovarian/patient2.tree data/ovarian/patient2.observed.labeling --colormap data/ovarian/coloring.txt
        

