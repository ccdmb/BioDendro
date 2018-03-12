# BioDendro - package to cluster genome data, and plot dendrograms

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

