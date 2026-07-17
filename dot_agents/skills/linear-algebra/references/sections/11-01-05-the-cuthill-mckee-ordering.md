# 11.1.5 The Cuthill-McKee Ordering

Because bandedness is such a tractable form of sparsity, it is natural to approach the Sparse Cholesky challenge by making $\tilde { A } = P A P ^ { \hat { T } }$ as “banded as possible” subject to cost constraints. However, this is too restrictive as Example 2 in §11.1.3 shows. Profile minimization is a better way to induce good sparsity in G. The profile of a symmetric $A \in \mathbb { R } ^ { n \times n }$ is defined by

$$
\operatorname{profile} (A) = n + \sum_ {i = 1} ^ {n} \left(i - f _ {i} (A)\right)
$$

where the profile indices $f _ { 1 } ( A ) , \ldots , f _ { n } ( A )$ are given by

$$
f _ {i} (A) = \min \{j: 1 \leq j \leq i, a _ {i j} \neq 0 \}. \tag {11.1.6}
$$

For the 9-by-9 example in (11.1.5), profile(A) = 37. We use that matrix to illustrate a heuristic method for approximate profile minimization. The first step is to choose a “starting node” and to relabel it as node 1. For reasons that are given later, we choose node 2 and set $S _ { 0 } = \{ 2 \}$ :

![](images/golub_600_649__0f8b7dfa4545846b14f08779418ea0097f457598db174b9f8cacfda46a6dae8e.jpg)

<details>
<summary>flowchart</summary>

```mermaid
graph TD
    7 --> 8
    7 --> 4
    8 --> 2
    4 --> 1
    4 --> 6
    4 --> 3
    1 --> 5
    5 --> 2
    6 --> 4
    3 --> 1
```
</details>

Original $\mathcal { G } _ { A }$

![](images/golub_600_649__a11742984f25d11cd19f6f958b1bc69731db60c9ea431b87c5069fd1283c6813.jpg)

<details>
<summary>flowchart</summary>

Simple undirected graph diagram with 8 nodes and 1 label, showing connections and a central node.
</details>

Labeled: $S _ { 0 }$

We then proceed to label the remaining nodes as follows:

Label the neighbors of $S _ { 0 }$ . Those neighbors make up $S _ { 1 }$ .

Label the unlabeled neighbors of nodes in $S _ { 1 }$ . Those neighbors make up $S _ { 2 }$

Label the unlabeled neighbors of nodes in $S _ { 2 }$ . Those neighbors make up $S _ { 3 }$ . etc.

If we follow this plan for the example, then $S _ { 1 } = \{ 8 , 5 \} , S _ { 2 } = \{ 7 , 3 \} , S _ { 3 } = \{ 1 , 4 \}$ , and $S _ { 4 } = \{ 6 , 9 \}$ . These are the level sets of node 2 and here is how they are determined one after the other:

![](images/golub_600_649__8b095ada894dac3554823f342a54785915d4bd6de2f0e861fdcd6670536d4c72.jpg)

<details>
<summary>flowchart</summary>

```mermaid
graph TD
    A[" "] ---_B["2"]
    A_---_C["1"]
    A_---_D["3"]
    A_---_E[" "]
    A_---_F[" "]
    A_---_G[" "]
    A_---_H[" "]
    A_---_I[" "]
    A_---_J[" "]
```
</details>

Labeled: S0, S1

![](images/golub_600_649__a0ccc5390cddb9b224ce0727783abebac80c729abb1973ff800b93903e59623a.jpg)

<details>
<summary>flowchart</summary>

```mermaid
graph TD
    A[" "] --> B[" "]
    A --> C[" "]
    A --> D[" "]
    B --> E[" "]
    B --> F[" "]
    C --> G[" "]
    C --> H[" "]
    D --> I[" "]
    D --> J[" "]
    E --> K[" "]
    F --> L[" "]
    G --> M[" "]
    H --> N[" "]
    I --> O[" "]
    J --> P[" "]
```
</details>

Labeled: S0, S1, S2

![](images/golub_600_649__c727550509de30a02361267e4958a0ac728604a1358207b32db8d2168015dc04.jpg)

<details>
<summary>flowchart</summary>

```mermaid
graph TD
    4 --> 2
    4 --> 7
    2 --> 1
    2 --> 3
    3 --> 6
    3 --> 5
    5 --> 6
    5 --> 7
    7 --> 4
    7 --> 6
    7 --> 5
    6 --> 7
    7 --> 6
    7 --> 5
```
</details>

$\mathrm { L a b e l e d } \colon S _ { 0 } , S _ { 1 } , S _ { 2 } , S _ { 3 }$

![](images/golub_600_649__5cb63412fb3515101cbbe7b629568fd8b5bf899c10f3a959e6205f97a6168e1f.jpg)

<details>
<summary>flowchart</summary>

```mermaid
graph TD
    4 --> 2
    4 --> 7
    2 --> 1
    2 --> 6
    1 --> 3
    3 --> 5
    3 --> 6
    3 --> 7
    7 --> 4
    7 --> 6
    7 --> 8
    7 --> 9
    6 --> 5
    6 --> 8
    6 --> 9
```
</details>

$_ { \mathrm { L a b e l e d : } ~ S _ { 0 } , ~ S _ { 1 } , ~ S _ { 2 } , ~ S _ { 3 } , ~ S _ { 4 } }$

By “concatenating” the level sets we obtain the Cuthill-McKee reordering :

![](images/golub_600_649__d8f15050a73ddb868c4d173366b94ec0e19d47270f22c60a636a9969de764448.jpg)

<details>
<summary>text_image</summary>

p: 2 | 8 | 5 | 7 | 3 | 1 | 4 | 6 | 9 .
\underbrace{S_0} _{S_1} \underbrace{S_2} _{S_3} \underbrace{S_4} _{S_4}
</details>

Observe the band structure that is induced by this ordering:

![](images/golub_600_649__454ffc8e924d33cbc38fb6007825552cf3434ca3ca130d078f98c99a5ae24390.jpg)

<details>
<summary>text_image</summary>

A(p,p) = 
⑨
④
2
1
7
6
3
8
5
(11.1.7)
</details>

Note that profile $\left( A ( p , p ) \right) = 2 5$ . Moreover, $A ( p , p )$ is a 5-by-5 block tridiagonal matrix with square diagonal blocks that have dimension equal to the cardinality of the level sets $S _ { 0 } , \ldots , S _ { 4 }$ . This suggests why a good choice for $S _ { 0 }$ is a node that has “far away” neighbors. Such a node will have a relatively large number of level sets and that means the resulting block tridiagonal matrix $A ( p , p )$ will have more diagonal blocks. Heuristically, these blocks will be smaller and that implies a tighter profile. See George and Liu (1981, Chap. 4) for a discussion of this topic and why the reverse Cuthill-McKee ordering $p ( n { : } { - } 1 { : } 1 )$ typically results in less fill-in during the Cholesky process.
