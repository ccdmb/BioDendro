# BioDendro - package to cluster metabolomics data, and plot dendrograms

## Background
Project owner: Catherine Rawlinson (PhD candidate)
Email: catherine.rawlinson@postgrad.curtin.edu.au

Converts MGF format and component list into non-redundant list.
Component-analyte list is converted into a data matrix and analytes are dynamically binned and clustered.

## Install

```
pip install git+https://github.com/CurtinIC/BioDendro.git

# or

git clone git@github.com:CurtinIC/BioDendro.git && cd BioDendro
pip install .
```

## Quick Start Example

### Jupyter Notebook file = quick-start-example.ipynb
```
#Load modules
import os
import plotly
import BioDendro

# Run the complete BioDendro pipeline
tree = BioDendro.pipeline("./MSMS.mgf", "./component_list.txt", clustering_method="braycurtis")

# View the new dendrogram cutoff inline
plotly.offline.init_notebook_mode(connected=True) # for visualising plot inline
plotly.offline.iplot(k, filename='simple_dendrogram')

```
### Explore BioDendro pipeline arguments 

```
# listed are all BioDendro pipeline arguments
!BioDendro -h
```

![Scheme](cluster-d10.png "Clustering")

