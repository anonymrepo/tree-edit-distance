This code implements computation of the tree edit distance using the algorithm proposed by [Zhang and Shasha][1].
 
The tree edit distance is defined as the minimal cost operation sequence of transforming one tree graph into another tree graph by applying one of the following operations: 
- Insertion of an additional node
- Deletion of a node
- Relabeling of a node

The costs of each operation can be provided as a parameter once the distance computation is invoked or be assumed to be unit for each of the operations.

Note that, the trees need to be ordered to ensure the algorithm yields sound results. 

[1]: https://doi.org/10.1137/0218082
