# BioDendro

**A package to cluster metabolomics data, and plot dendrograms**

## Background

- Project owner: Catherine Rawlinson (PhD candidate)
- Email: catherine.rawlinson@postgrad.curtin.edu.au

Converts MGF format and component list into non-redundant list.
Component-analyte list is converted into a data matrix and analytes are dynamically binned and clustered.

## Install

Installing is easiest with pip.
Assuming you have python3 installed you can run the following to install.

```bash
pip3 install --user git+https://github.com/CurtinIC/BioDendro.git

# or

git clone git@github.com:CurtinIC/BioDendro.git && cd BioDendro
pip install --user .
```

To install as root, you can omit `--user`. I.E.

```bash
sudo pip3 install git+https://github.com/CurtinIC/BioDendro.git
```

Both the `BioDendro` script and the python package will now be installed on your path
(assuming Python is configured correctly).


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
