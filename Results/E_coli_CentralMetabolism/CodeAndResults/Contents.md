#Contents of the files in current directory

1. DAGFinder_ecoliiJO-1.py - Main python file that contains the implementation of SubNetHunter. Consists of the list of source metabolites, target metabolites, seed metabolites, cutoff to generate the subnetworks and the number of subnetworks to be generated

2. forward_pass.py - Python implementation of the Guided BFS (Refer main manuscript for the details)

3. generate_partitions.py - Python function file that generates parititions of numbers for generating a given sum. This function is called in the main implementation (1)

4. ecoliiJO1366_15_iJO1366 pyr_c.pickle - Pickle file consisting of all the sub-networks found for all the metabolites for different path-lengths

5. glc_DASH_D_e_iJO1366 pyr_c_15.txt - Text file consisting of all the sub-networks that lead to pyruvate from the starting seed metabolites. Note that the file follows the format given below:
SubNetwork number
Reaction Name	List of reactants	List of products

6. onlyglc_DASH_D_e_iJO1366 pyr_c_15.txt - Text file consisting of all the subnetworks that lead to pyruvate from the glucose. Note that the file follows the format given below:
SubNetwork number
Reaction Name	List of reactants	List of products
Glycolysis in this file - SubNetwork 108
(Note that both the forward and the reverse reactions are assigned the same reaction IDs in the subnetwork)
