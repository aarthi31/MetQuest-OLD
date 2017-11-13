This folder contains all the input data used to identify sub-networks in E. coli iJO1366 model,i.e.,

1. Genome-scale metabolic model *i*JO1366 in *SBML* format
2. Directed Bipartite graph of iJO1366 as a ".gpickle" file. This file can be read in python by using the following command:

```python
from networkx as import read_gpickle
	G = read_gpickle("iJO1366.gpickle")
```
3. List of seed metabolites including the source and the target metabolites
4. Cutoff - 15
