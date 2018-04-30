<<<<<<< HEAD
# BioDendro - package to cluster analyte data, and plot dendrograms
=======
# BioDendro - package to cluster metabolomics data, and plot dendrograms
>>>>>>> d9527d31d45651208ca8c6ffc0f977849b5a5c1d

## Background

## Example

```
#Import BioDendro functions
from BioDendro import *

data=Dendrogram('out.xlsx') #Load the data file
data.clusterize() 
data.generate_linkage() #takes long
data.generate_out()  #plots and .txt saved to results

```

