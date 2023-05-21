# MACH2

A mathematical framework for inferring migration histories of metastatic cancer from clone phylogeny and the location of extant clones.

## Table of contents

1. Installation
        - Prerequisite
        - Install using `pip`
2. Usage instruction
        - I/O formats
        - Input
        - Output
        - Usage

## 1. Installation

### Prerequisites

- **Python** - `MACH2` requires Python 3.7 or newer.
- **ILP solver** - `MACH2` requires an ILP solver installed to solve **PMH-TR**. Currently `MACH2` only supports `Gurobi optimizer`, but we are going to add support for more ILP solvers in the future. `MACH2` requires a valid Gurobi installation and license key. The location of Gurobi should be present in `LD_LIBRARY_PATH` (linux) and `DYLD_LIBRARY_PATH` (macOS) the license key should be saved in the environment variable `GRB_LICENSE_KEY`.

### Install using `pip`

`MACH2` can easily be installed using `pip`, the package installer for Python. Open a terminal or command prompt and run the following command:

                $ pip install mach2

If you want to use `MACH2` with `JupyterLab`, you'll need additional dependencies. To install these optional dependencies, you can run the following command:

                $ pip install mach2[jupyter]


## Usage Instruction

### I/O formats

We describe various formats used by `MACH2`.

1. **Tree file** : The tree file contains a list of edges that define the structure of a tree. Each line in the file represents an edge, and the edges should be in the format: `vertex1 vertex2`. For example:

        1   2
        2   3
        2   4 
        1   3

2. **Labeling file** : The labeling file contains the labels assigned to a set of vertices. Each line in the file corresponds to a vertex and the labels are in the format: `leaf label`. For example:

        1   A
        2   B
        3   C

1. **Multi-graph file** : Like tree file, the multi-graph file contains a list of edges that define the structure of the multi-graph. Each line in the file represents the set of edges between any pair of vertices __vertex1__ and __vertex2__ in the format: `vertex1 vertex2 #edges`. For example:

        A   B   2
        B   A   1
        A   C   1

Additionaly, `MACH2` can output files in Graphviz DOT format or JSON format.

### Input

`MACH2` takes as input two files - 

1. **Input tree file** : Tree file describing the input clone tree.
2. **Leaf labeling file** : Labeling file describing the leaf labeling of input clone tree.

### Output

For each solution, `MACH2` can output five types of files.

1. **Refined tree file** : Tree file describing the output refined tree. 
2. **Vertex labeling file** : Labeling file describing the vertex labeling of the refined tree.
3. **Refined tree DOT** : Refined tree with vertex labeling in DOT format.
4. **Migration graph file** : Multi-graph file describing the migration graph.
5. **Migration graph DOT** : Migration graph in DOT format.

Additionaly `MACH2` can return JSON file encoding all the solutions (refined tree with vertex labeling and migration graph). The JSON file can be directly passed to [MACH2-viz](https://github.com/elkebir-group/mach2-viz). The exact format of the JSON file is described [here](https://github.com/elkebir-group/mach2-viz/README.md).

`MACH2` also prints `<primary location> <number of migrations> <number of comigrations> Optimal <running time (in seconds)>` on console.

### Usage

`MACH2` can be run using python.

                usage: mach2 [-h] [-p PRIMARY] [-c COLORMAP] [--log] [-o OUTPUT] [-N NSOLUTIONS] [-C] [-t THREADS] [-s] [-S] clone_tree leaf_labeling

                MACH2

                positional arguments:
                clone_tree            Input clone tree
                leaf_labeling         Input leaf labeling

                options:
                -h, --help            show this help message and exit
                -p PRIMARY, --primary PRIMARY
                                        Primary anatomical site
                -c COLORMAP, --colormap COLORMAP
                                        Color map file
                --log                 Outputs Gurobi logging
                -o OUTPUT, --output OUTPUT
                                        Output folder
                -N NSOLUTIONS, --nsolutions NSOLUTIONS
                                        Maximum number of solutions retained
                -C, --count_solutions
                                        Only prints the number of solutions (default=False)
                -t THREADS, --threads THREADS
                                        Number of threads
                -s, --suboptimal      Returns suboptimal solutions without duplicates, may be slow (default=False)
                -S, --seeding_sites   Minimizes the number of seeding sites too (default=False)

An example execution

        $ mach2 data/mcpherson_2016/patient1.tree data/mcpherson_2016/patient1.labeling -c data/mcpherson_2016/coloring.txt
        

