# BioDendro

**A package to cluster and visualise MS/MS spectral data**

[![Build Status](https://travis-ci.org/ccdmb/BioDendro.svg?branch=master)](https://travis-ci.org/ccdmb/BioDendro)
[![PyPI version](https://badge.fury.io/py/BioDendro.svg)](https://badge.fury.io/py/BioDendro)
[![Anaconda-Server Badge](https://anaconda.org/darcyabjones/biodendro/badges/version.svg)](https://anaconda.org/darcyabjones/biodendro)

<div align="center">
    <img src="https://github.com/ccdmb/BioDendro/blob/master/images/banner.png" alt="" width="300" />
</div>

## Background

- Project owner: Catherine Rawlinson (PhD candidate)
- Email: catherine.rawlinson@postgrad.curtin.edu.au


BioDendro is a metabolomics package and workflow that enables analysts to flexibly cluster
and interrogate thousands of MS/MS spectra and quickly identify the core fragment
patterns causing groupings.
This helps identify potential functional properties of components based on core
chemical backbones of a larger class, even when the individual metabolite of
interest is not found in public databases.

BioDendro takes raw MS/MS data in MGF format, and a component list.
The components list is the total of analytes within your sample set in the following format...

```
SampleID_userinfo_userinfo_m/z_RT
```

With retention time (RT) in units of minutes.
For example:

```
Sample1_pos_C18_123.1234_5.60
Sample2_pos_C19_321.4321_10.60
```

This can be generated using [XCMS](https://xcmsonline.scripps.edu/landing_page.php?pgcontent=mainPage), [MZmine2](http://mzmine.github.io/), or proprietary instument vendor software.

Converts MGF format and component list into non-redundant list.
Component-analyte list is converted into a data matrix and analytes are dynamically binned and clustered.


## Requirements

- Python version 3.5 or more recent.
- The python packages numpy, pandas, scipy, matplotlib, plotly, xlrd, xlsxwriter, and pillow (Installed automatically).
- We recommend running the pipeline in [Jupyter notebooks](https://jupyter.org/), and provide example notebooks.

BioDendro is tested to run with Python 3.5-3.7, Plotly 3.8 and 3.9, and Pandas 0.23 and 0.24.
Other versions may work.


## Install on Windows

A detailed guide to installing Anaconda Python, Jupyter, and BioDendro on Windows operating systems is provided in a [separate pdf file](https://github.com/ccdmb/BioDendro/blob/master/Download%20and%20install%20instructions%20for%20BioDendro%20using%20Windows%2010.pdf).
Advanced users with some knowledge of Python may also use the command line installation instructions below.


## Install from the command line

BioDendro can be installed from [PyPI](https://pypi.org/project/BioDendro) using [pip](https://pip.pypa.io/en/stable/), or from [anaconda](https://anaconda.org/darcyabjones/biodendro) using [conda](https://docs.conda.io/en/latest/).

Users that are less familiar with Python and pip are recommended to read our [INSTALLING_WITH_PIP.md](https://github.com/ccdmb/BioDendro/blob/master/INSTALLING_WITH_PIP.md) document which explains things in more detail, including where things will be installed and how to use virtual environments.
For details on installing and using conda, see their [getting-started guide](https://docs.conda.io/projects/conda/en/latest/user-guide/overview.html).


Assuming you have Python 3 installed you can install BioDendro and its dependencies using the pip:

```bash
python3 -m pip install --user BioDendro

# Only required if you're using the provided notebooks
python3 -m pip install --user jupyter
```


To install BioDendro and dependencies using conda (assuming you have installed Anaconda):

```bash
conda install -c darcyabjones biodendro

# Only required if you're using the provided notebooks
conda install jupyter
```


Both the `BioDendro` script and the python package (which can be used with the notebooks) should now be available to use.

Note that the above commands will not download the example notebooks or data.
You can download those files separately, or download the whole repository as recommended in the windows install guide.



## Quick Start - Python library

The pipeline available as a python package.
To run the full pipeline in python.

```python
import BioDendro

tree = BioDendro.pipeline(
    "Fireflies_MSMS.mgf",
    "Fireflies_feature_list.txt",
    results_dir="my_results_dir"
)
```

From there you could analyse the results stored in the `tree` object.
The example jupyter notebooks contain more detailed explanations of different parameters.

[quick-start-example.ipynb](https://github.com/ccdmb/BioDendro/blob/master/quick-start-example.ipynb) contains basic information about running the pipelines.

[longer-workflow.ipynb](https://github.com/ccdmb/BioDendro/blob/master/longer-workflow.ipynb) contains more detailed information about how the pipeline works, and how you can modify parameters.

We suggest that beginners download the [quick-start notebook available here](https://github.com/ccdmb/BioDendro/raw/master/quick-start-example.ipynb) (Right-click, save-as) and modify parameters and files as necessary.


## Quick Start - command line

The pipeline is also available as a command line script.
This is useful if you're not planning on tweaking the parameters much and just want to run the darn thing.

A list of options can be obtained with the `--help` (or `-h`) flag.

```bash
BioDendro --help
```

The minimum options to run the pipeline are the MGF file and a components list.
To run the basic pipeline with the same parameters as in the python quickstart:

```bash
BioDendro --results-dir my_results_dir MSMS.mgf component_list.txt
```

If the `--results-dir` parameter is ommitted, the results will be stored in a directory with the current date and the current time added to the end of it.

You can change the parameters to use by supplying additional flags, however, this will run the whole pipeline again, so it you just need to adjust the cutoff or decide to use braycurtis instead of jaccard distances, you might be better off using the python API.

```bash
BioDendro --scaling --cluster-method braycurtis --cutoff 0.5 MSMS.mgf component_list.txt
```

would be equivalent to running the following in python

```python
tree = BioDendro.pipeline("MSMS.mgf", "component_list.txt", clustering_method="braycurtis", scaling=True, cutoff=0.5)
```
