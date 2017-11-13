# MetQuest

MetQuest is a dynamic programming based algorithm for exhaustively identifying sub-networks from metabolic networks. Given a set of source and target set of metabolites, the algorithm seeks to identify the sub-networks of reactions that aid in their conversions.

## Input
The input to the algorithm are the SBML file(s) of the metabolic networks of organisms, _seed metabolites_, cut-off &beta, number of such sub-networks required. The implementation will automatically convert the SBML files to a directed bipartite graph representation. When more than one metabolic network needs to be used, the implementation will introduce a common extracellular environment through which the metabolic interactions happen. Note that the metabolites are renamed according to the model names, in order to easily distinguish the metabolites between species.

## Output

MetQuest first calculates the scope _M_c_ of the source metabolites, i.e.,the metabolites that can be produced from the given set of seed metabolites. The algorithm then initializes and maintains a table of size |_M_c_|,   
