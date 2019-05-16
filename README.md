# BioDendro

**A package to cluster metabolomics data, and plot dendrograms**

[![Build Status](https://travis-ci.org/ccdmb/BioDendro.svg?branch=master)](https://travis-ci.org/ccdmb/BioDendro)

## Background

- Project owner: Catherine Rawlinson (PhD candidate)
- Email: catherine.rawlinson@postgrad.curtin.edu.au

Converts MGF format and component list into non-redundant list.
Component-analyte list is converted into a data matrix and analytes are dynamically binned and clustered.

## Install on Linux or Mac from bash

Installing is easiest with pip.
Assuming you have python3 installed you can run the following to install.

```bash
python3 -m pip install --user biodendro

# or

git clone git@github.com:CurtinIC/BioDendro.git && cd BioDendro
python3 -m pip install --user biodendro
```

The `--user` flag tells pip to install to a user directory rather than a system directory.
Generally this will be under `~/.local` for Mac and Linux.
Make sure that `~/.local/bin` is added to your `$PATH` if this is the case [see here](https://stackoverflow.com/questions/14637979/how-to-permanently-set-path-on-linux-unix).

For more advanced users, we recommend using a [virtual environment](https://virtualenv.pypa.io/en/stable/) or [conda environment](https://docs.conda.io/projects/conda/en/latest/user-guide/tasks/manage-environments.html).

To install as root, you can omit `--user`, though this is generally discouraged.

```bash
sudo python3 -m pip install biodendro
```

To install the latest and greatest version, you can use [git](https://git-scm.com/), to install directly from the repository.

```bash
python3 -m pip install --user git+https://github.com/CurtinIC/BioDendro.git

# or

git clone git@github.com:CurtinIC/BioDendro.git && cd BioDendro
python3 -m pip install --user .
```

Both the `BioDendro` script and the python package will now be available to use (assuming Python is configured correctly).


## Quick Start Example - command line

The quickest way to run is using the command-line interface.

A list of options can be obtained with the `--help` flag.

```bash
BioDendro --help
```

To run the basic pipeline using the example MGF and components file do:

```bash
BioDendro --results-dir my_results_dir MSMS.mgf component_list.txt
```


## Quick Start Example - Python library

The pipeline is also available as a python function/library.
The command above would be equivalent to the following in python.

```python
import BioDendro

tree = BioDendro.pipeline("MSMS.mgf", "component_list.txt", results_dir="my_results_dir")
```

From there you could analyse the results stored in `tree`.
The example jupyter notebooks contain more detailed explanations of different parameters.

[quick-start-example.ipynb](quick-start-example.ipynb) contains basic information about running the pipelines.

[longer-example.ipynb](longer-example.ipynb) contains more detailed information about how the pipeline works, and how you can modify parameters.


![Scheme](cluster-d10.png "Clustering")


## Command line API

The pipeline can also be run from a bash or bash-like terminal.
This is useful if you're not planning on tweaking the parameters much and just want to run the darn thing.

For these examples, we're using the ipython magic command `%%bash` to run the commands in bash.
You can omit the %%bash bit if you're running straight in the terminal.

To get a list of all options available use the `--help` (or `-h`) flag.

```bash
%%bash
BioDendro --help
```
The minimum options to run the pipeline are the MGF file and a components list.

Using the example data in the BioDendro repo we could run...

```bash
%%bash

BioDendro MSMS.mgf component_list.txt
```

As before, the results will be stored in a directory with the current date and the current time added to the end of it.

You can change the parameters to use by supplying additional flags, however, this will run the whole pipeline again, so it you just need to adjust the cutoff or decide to use braycurtis instead of jaccard distances, you might be better off using the python API.

```bash
%%bash

BioDendro --scaling --cluster-method braycurtis --cutoff 0.5 MSMS.mgf component_list.txt
```

would be equivalent to running the following in python

```python
tree = BioDendro.pipeline("MSMS.mgf", "component_list.txt", clustering_method="braycurtis", scaling=True, cutoff=0.5)
```
