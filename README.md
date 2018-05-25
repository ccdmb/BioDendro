# BioDendro - package to cluster metabolomics data, and plot dendrograms

## Background
Project owner: Catherine Rawlinson (PhD candidate)
Email: catherine.rawlinson@curtin.edu.au

Converts MGF format and component list into non-redundant list.
Component-analyte list is converted into a data matrix and analytes are dynamically binned and clustered.

## Example


### Prepare data 
```
msms-prepare-input.py -i MSMS_CAT.mgf -j component_list_from_12CVs.txt 
```
### Import BioDendro functions
```
from BioDendro import *

data=Dendrogram('out.xlsx') #Load the data file
k=data.visualize(cutoff=5.0) # Set cluster distance 
```
### To visualise a plot inline
```
plotly.offline.init_notebook_mode(connected=True) # for visualising plot inline
plotly.offline.iplot(k,filename='simple_dendrogram')
```

![Scheme](cluster-d10.png "Clustering")

## BioDendro workflow

### Sample analyte de-replication

option neutral loss set option -n
```
msms-prepare-input.py MSMS.mgf component_list.txt <-n>
```

### Analyte dynamic binning

```
data.clusterize(bin_threshold=0.8E-3)
```

### Sample linkage and clustering

```
data.generate_linkage(self,cutoff=0.5,clustering_method='jaccard')
```

### Analysis results and visualization

```
data.generate_out()
data.visualize(cutoff=0.5,x=900,y=400)
```
