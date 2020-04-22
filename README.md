# Qubit-Mapping-Subgraph-Isomorphism
This repository is the Python implementation for algorithms introduced in "Qubit Mapping Based on Subgraph Isomorphism and Filtered Limited-Depth Search" (https://arxiv.org/abs/2004.07138). Readers can run "QTokyo_run.py" to recur results in the paper.

We provide the following selections:

A. QFilter_type (line 60)
	One can select type from {'0', '1', '12', '12x', '2x'}
		 '0 :: if no filter is used'
		 '1 :: if Q0 is used for all layers'
		'12 :: if Q0 is used for the first layer and Q1 is used for the other layers'			
	       '12x :: if Q0 is used for the first layer and Q0+Q1 is used for the other layers'
		'2x :: if Q0+Q1 is used for all layers '
	where Q0 and Q1 are the set of qubits in the front layer and the first look-ahead layer of the current logical circuit.
	The default filter type is '12'. 

B. Initial mapping (line 67)
	One can select mapping from {'topgraph', 'wgtgraph', 'empty', 'naive'}
	The default mapping is 'topograph'

C. Size of the circuits (line 69)
	One can select size from {'small', 'medium', 'large', 'all'}
	The default size is 'small'

Our results show that Filter_type '12' is the fastest and the quality of Filter_type '2x' is very close to that of '0' (i.e., search without filter).
