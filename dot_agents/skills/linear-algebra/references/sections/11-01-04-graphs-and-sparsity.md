# 11.1.4 Graphs and Sparsity

Here is a sparse symmetric matrix A and its adjacency graph $\mathcal { G } _ { A }$ :

![](images/golub_600_649__9f155fb8547debf46269e7261c73f66fe33c481dc95018906c7c0408aedd5236.jpg)

In an adjacency graph for a symmetric matrix, there is a node for each row, numbered by the row number, and there is an edge between node i and node j if the off-diagonal entry $a _ { i j }$ is nonzero. In general, a graph $\mathcal { G } ( V , E )$ is a set of labeled nodes V together with a set of edges E, e.g.,

$$
V = \{1, 2, 3, 4, 5, 6, 7, 8, 9 \},
$$

$$
E = \{(1, 4), (1, 6), (1, 7), (2, 5), (2, 8), (3, 4), (3, 5), (4, 6), (4, 7), (4, 9), (5, 8), (7, 8) \}.
$$

Adjacency graphs for symmetric matrices are undirected. This means there is no difference between edge $( i , j )$ and edge $( j , i )$ . If P is a permutation matrix, then, except for vertex labeling, the adjacency graphs for A and $\bar { P } A P ^ { T }$ “look the same.”

Node i and node j are neighbors if there is an edge between them. The adjacency set for a node is the set of its neighbors and the cardinality of that set is the degree of the node. For the above example we have

<table><tr><td>Node</td><td>1</td><td>2</td><td>3</td><td>4</td><td>5</td><td>6</td><td>7</td><td>8</td><td>9</td></tr><tr><td>Degree</td><td>3</td><td>2</td><td>2</td><td>5</td><td>3</td><td>2</td><td>3</td><td>3</td><td>1</td></tr></table>

Graph theory is a very powerful language that facilitates reasoning about sparse matrix factorizations. Of particular importance is the use of graphs to predict structure, something that is critical to the design of efficient implementations. For a much deeper appreciation of these issues than what we offer below, see George and Liu (1981), Duff, Erisman, and Reid (1986), and Davis (2006).
