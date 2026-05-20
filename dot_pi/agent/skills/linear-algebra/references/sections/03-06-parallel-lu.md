# 3.6 Parallel LU

In §3.2.11 we show how to organize a block version of Gaussian elimination (without pivoting) so that the overwhelming majority of flops occur in the context of matrix multiplication. It is possible to incorporate partial pivoting and maintain the same level-3 fraction. After stepping through the derivation we proceed to show how the process can be effectively parallelized using the block-cyclic distribution ideas that were presented in §1.6.

# 3.6.1 Block LU with Pivoting

Throughout this section assume $A \in \mathbb { R } ^ { n \times n }$ and for clarity that $n = r N$

$$
A = \left[ \begin{array}{c c c} A _ {1 1} & \dots & A _ {1 N} \\ \vdots & \ddots & \vdots \\ A _ {N 1} & \dots & A _ {N N} \end{array} \right] \quad A _ {i, j} \in \mathbb {R} ^ {r \times r}. \tag {3.6.1}
$$

We revisit Algorithm 3.2.4 (nonrecursive block LU) and show how to incorporate partial pivoting.

The first step starts by applying scalar Gaussian elimination with partial pivoting to the first block column. Using an obvious rectangular matrix version of Algorithm 3.4.1 we obtain the following factorization:

$$
P _ {1} \left[ \begin{array}{c} A _ {1 1} \\ A _ {2 1} \\ \vdots \\ A _ {N 1} \end{array} \right] = \left[ \begin{array}{c} L _ {1 1} \\ L _ {2 1} \\ \vdots \\ L _ {N 1} \end{array} \right] U _ {1 1}. \tag {3.6.2}
$$

In this equation, $P _ { 1 } \in \mathbb { R } ^ { n \times n }$ is a permutation, $L _ { 1 1 } \in \mathbb { R } ^ { r \times r }$ is unit lower triangular, and $U _ { 1 1 } \in \mathbb { R } ^ { r \times r }$ is upper triangular.

The next task is to compute the first block row of U . To do this we set

$$
P _ {1} A = \left[ \begin{array}{c c c} \tilde {A} _ {1 1} & \dots & \tilde {A} _ {1 N} \\ \vdots & \ddots & \vdots \\ \tilde {A} _ {N 1} & \dots & \tilde {A} _ {N N} \end{array} \right], \quad \tilde {A} _ {i, j} \in \mathbb {R} ^ {r \times r}, \tag {3.6.3}
$$

and solve the lower triangular multiple-right-hand-side problem

$$
L _ {1 1} \left[ U _ {1 2} \mid \dots \mid U _ {1 N} \right] = \left[ \tilde {A} _ {1 2} \mid \dots \mid \tilde {A} _ {1 N} \right] \tag {3.6.4}
$$

for $U _ { 1 2 } , \dots , U _ { 1 N } \in \mathbb { R } ^ { r \times r }$ . At this stage it is easy to show that we have the partial factorization

$$
P _ {1} A = \left[ \begin{array}{c c c c} L _ {1 1} & 0 & \dots & 0 \\ \hline L _ {2 1} & I _ {r} & \dots & 0 \\ \vdots & \vdots & \ddots & \vdots \\ L _ {N 1} & 0 & \dots & I _ {r} \end{array} \right] \left[ \begin{array}{c c} I _ {r} & 0 \\ \hline 0 & A ^ {(\mathrm{new})} \end{array} \right] \left[ \begin{array}{c c c c} U _ {1 1} & U _ {1 2} & \dots & U _ {1 N} \\ \hline 0 & I _ {r} & \dots & 0 \\ \vdots & \vdots & \ddots & \vdots \\ 0 & 0 & \dots & I _ {r} \end{array} \right]
$$

where

$$
A ^ {(\text { new })} = \left[ \begin{array}{c c c} \tilde {A} _ {2 2} & \dots & \tilde {A} _ {2 N} \\ \vdots & \ddots & \vdots \\ \tilde {A} _ {N 2} & \dots & \tilde {A} _ {N N} \end{array} \right] - \left[ \begin{array}{c} L _ {2 1} \\ \vdots \\ L _ {N 1} \end{array} \right] [ U _ {1 2} | \dots | U _ {1 N} ]. \tag {3.6.5}
$$

Note that the computation of $A ^ { \mathrm { ( n e w ) } }$ is a level-3 operation as it involves one matrix multiplication per A-block.

The remaining task is to compute the pivoted LU factorization of $A ^ { ( \mathrm { n e w } ) }$ . Indeed, if

$$
P ^ {(\mathrm{new})} A ^ {(\mathrm{new})} = L ^ {(\mathrm{new})} U ^ {(\mathrm{new})}
$$

and

$$
P ^ {(\mathrm{new})} \left[ \begin{array}{c} L _ {2 1} \\ \vdots \\ L _ {N 1} \end{array} \right] = \left[ \begin{array}{c} \tilde {L} _ {2 1} \\ \vdots \\ \tilde {L} _ {N 1} \end{array} \right],
$$

then

$$
P A = \left[ \begin{array}{c c} L _ {1 1} & 0 \quad \dots \quad 0 \\ \hline \tilde {L} _ {2 1} & \\ \vdots & L ^ {(\text {new})} \\ \tilde {L} _ {N 1} & \end{array} \right] \left[ \begin{array}{c c} U _ {1 1} & U _ {1 2} \dots U _ {1 N} \\ \hline 0 & \\ \vdots & U ^ {(\text {new})} \\ 0 & \end{array} \right]
$$

is the pivoted block LU factorization of A with

$$
P = \left[ \begin{array}{c c} I _ {r} & 0 \\ 0 & P ^ {(\mathrm{new})} \end{array} \right] P _ {1}.
$$

In general, the processing of each block column in A is a four-part calculation:

Part A. Apply rectangular Gaussian Elimination with partial pivoting to a block column of A. This produces a permutation, a block column of $L _ { : }$ and a diagonal block of U . See (3.6.2).

Part B. Apply the Part A permutation to the “rest of $A . ^ { \mathfrak { n } }$ See (3.6.3).

Part C. Complete the computation of $U \mathrm { { ^ { * } s } }$ next block row by solving a lower triangular multiple right-hand-side problem. See (3.6.4).

Part D. Using the freshly computed L-blocks and U-blocks, update the “rest of $A . ^ { \mathfrak { n } }$ See (3.6.5).

The precise formulation of the method with overwriting is similar to Algorithm 3.2.4 and is left as an exercise.

# 3.6.2 Parallelizing the Pivoted Block LU Algorithm

Recall the discussion of the block-cyclic distribution in §1.6.2 where the parallel computation of the matrix multiplication update $C = C + A B$ was outlined. To provide insight into how the pivoted block LU algorithm can be parallelized, we examine a representative step in a small example that also makes use of the block-cyclic distribution.

Assume that $N = 8$ in (3.6.1) and that we have a $p _ { \mathrm { r o w } } { - } \mathrm { b y } { - } p _ { \mathrm { c o l } }$ processor network with $p _ { \mathrm { r o w } } = 2$ and $p _ { \mathrm { c o l } } ~ = ~ 2$ . At the start, the blocks of $A \ = \ ( A _ { i j } )$ are cyclically distributed as shown in Figure 3.6.1. Assume that we have carried out two steps of block $L U$ and that the computed $L _ { i j }$ and $U _ { i j }$ have overwritten the corresponding $A -$ blocks. Figure 3.6.2 displays the situation at the start of the third step. Blocks that are to participate in the Part A factorization

$$
P _ {3} \left[ \begin{array}{c} A _ {3 3} \\ \vdots \\ A _ {8 3} \end{array} \right] = \left[ \begin{array}{c} L _ {3 3} \\ \vdots \\ L _ {8 3} \end{array} \right] U _ {3 3}
$$

are highlighted. Typically, $p _ { \mathrm { r o w } }$ processors are involved and since the blocks are each $r { \mathrm { - } } \mathrm { b y } { \mathrm { - } } r$ , there are r steps as shown in (3.6.6).

<table><tr><td>Proc(0,0)</td><td>Proc(0,1)</td><td>Proc(0,0)</td><td>Proc(0,1)</td><td>Proc(0,0)</td><td>Proc(0,1)</td><td>Proc(0,0)</td><td>Proc(0,1)</td></tr><tr><td>A11</td><td>A12</td><td>A13</td><td>A14</td><td>A15</td><td>A16</td><td>A17</td><td>A18</td></tr><tr><td>Proc(1,0)</td><td>Proc(1,1)</td><td>Proc(1,0)</td><td>Proc(1,1)</td><td>Proc(1,0)</td><td>Proc(1,1)</td><td>Proc(1,0)</td><td>Proc(1,1)</td></tr><tr><td>A21</td><td>A22</td><td>A23</td><td>A24</td><td>A25</td><td>A26</td><td>A27</td><td>A28</td></tr><tr><td>Proc(0,0)</td><td>Proc(0,1)</td><td>Proc(0,0)</td><td>Proc(0,1)</td><td>Proc(0,0)</td><td>Proc(0,1)</td><td>Proc(0,0)</td><td>Proc(0,1)</td></tr><tr><td>A31</td><td>A32</td><td>A33</td><td>A34</td><td>A35</td><td>A36</td><td>A37</td><td>A38</td></tr><tr><td>Proc(1,0)</td><td>Proc(1,1)</td><td>Proc(1,0)</td><td>Proc(1,1)</td><td>Proc(1,0)</td><td>Proc(1,1)</td><td>Proc(1,0)</td><td>Proc(1,1)</td></tr><tr><td>A41</td><td>A42</td><td>A43</td><td>A44</td><td>A45</td><td>A46</td><td>A47</td><td>A48</td></tr><tr><td>Proc(0,0)</td><td>Proc(0,1)</td><td>Proc(0,0)</td><td>Proc(0,1)</td><td>Proc(0,0)</td><td>Proc(0,1)</td><td>Proc(0,0)</td><td>Proc(0,1)</td></tr><tr><td>A51</td><td>A52</td><td>A53</td><td>A54</td><td>A55</td><td>A56</td><td>A57</td><td>A58</td></tr><tr><td>Proc(1,0)</td><td>Proc(1,1)</td><td>Proc(1,0)</td><td>Proc(1,1)</td><td>Proc(1,0)</td><td>Proc(1,1)</td><td>Proc(1,0)</td><td>Proc(1,1)</td></tr><tr><td>A61</td><td>A62</td><td>A63</td><td>A64</td><td>A65</td><td>A66</td><td>A67</td><td>A68</td></tr><tr><td>Proc(0,0)</td><td>Proc(0,1)</td><td>Proc(0,0)</td><td>Proc(0,1)</td><td>Proc(0,0)</td><td>Proc(0,1)</td><td>Proc(0,0)</td><td>Proc(0,1)</td></tr><tr><td>A71</td><td>A72</td><td>A73</td><td>A74</td><td>A75</td><td>A76</td><td>A77</td><td>A78</td></tr><tr><td>Proc(1,0)</td><td>Proc(1,1)</td><td>Proc(1,0)</td><td>Proc(1,1)</td><td>Proc(1,0)</td><td>Proc(1,1)</td><td>Proc(1,0)</td><td>Proc(1,1)</td></tr><tr><td>A81</td><td>A82</td><td>A83</td><td>A84</td><td>A85</td><td>A86</td><td>A87</td><td>A88</td></tr></table>

Figure 3.6.1.

Part A:

<table><tr><td>Proc(0,0)</td><td>Proc(0,1)</td><td>Proc(0,0)</td><td>Proc(0,1)</td><td>Proc(0,0)</td><td>Proc(0,1)</td><td>Proc(0,0)</td><td>Proc(0,1)</td></tr><tr><td>U11</td><td>U12</td><td>U13</td><td>U14</td><td>U15</td><td>U16</td><td>U17</td><td>U18</td></tr><tr><td>Proc(1,0)</td><td>Proc(1,1)</td><td>Proc(1,0)</td><td>Proc(1,1)</td><td>Proc(1,0)</td><td>Proc(1,1)</td><td>Proc(1,0)</td><td>Proc(1,1)</td></tr><tr><td>L21</td><td>U22</td><td>U23</td><td>U24</td><td>U25</td><td>U26</td><td>U27</td><td>U28</td></tr><tr><td>Proc(0,0)</td><td>Proc(0,1)</td><td>Proc(0,0)</td><td>Proc(0,1)</td><td>Proc(0,0)</td><td>Proc(0,1)</td><td>Proc(0,0)</td><td>Proc(0,1)</td></tr><tr><td>L31</td><td>L32</td><td>A33</td><td>A34</td><td>A35</td><td>A36</td><td>A37</td><td>A38</td></tr><tr><td>Proc(1,0)</td><td>Proc(1,1)</td><td>Proc(1,0)</td><td>Proc(1,1)</td><td>Proc(1,0)</td><td>Proc(1,1)</td><td>Proc(1,0)</td><td>Proc(1,1)</td></tr><tr><td>L41</td><td>L42</td><td>A43</td><td>A44</td><td>A45</td><td>A46</td><td>A47</td><td>A48</td></tr><tr><td>Proc(0,0)</td><td>Proc(0,1)</td><td>Proc(0,0)</td><td>Proc(0,1)</td><td>Proc(0,0)</td><td>Proc(0,1)</td><td>Proc(0,0)</td><td>Proc(0,1)</td></tr><tr><td>L51</td><td>L52</td><td>A53</td><td>A54</td><td>A55</td><td>A56</td><td>A57</td><td>A58</td></tr><tr><td>Proc(1,0)</td><td>Proc(1,1)</td><td>Proc(1,0)</td><td>Proc(1,1)</td><td>Proc(1,0)</td><td>Proc(1,1)</td><td>Proc(1,0)</td><td>Proc(1,1)</td></tr><tr><td>L61</td><td>L62</td><td>A63</td><td>A64</td><td>A65</td><td>A66</td><td>A67</td><td>A68</td></tr><tr><td>Proc(0,0)</td><td>Proc(0,1)</td><td>Proc(0,0)</td><td>Proc(0,1)</td><td>Proc(0,0)</td><td>Proc(0,1)</td><td>Proc(0,0)</td><td>Proc(0,1)</td></tr><tr><td>L71</td><td>L72</td><td>A73</td><td>A74</td><td>A75</td><td>A76</td><td>A77</td><td>A78</td></tr><tr><td>Proc(1,0)</td><td>Proc(1,1)</td><td>Proc(1,0)</td><td>Proc(1,1)</td><td>Proc(1,0)</td><td>Proc(1,1)</td><td>Proc(1,0)</td><td>Proc(1,1)</td></tr><tr><td>L81</td><td>L82</td><td>A83</td><td>A84</td><td>A85</td><td>A86</td><td>A87</td><td>A88</td></tr></table>

Figure 3.6.2.

for $j = 1 { : } r$

Columns $A _ { k k } ( : , j ) , \dotsc , A _ { N , k } ( : , j )$ are assembled in the processor housing $A _ { k k }$ , the “pivot processor”

The pivot processor determines the required row interchange and the Gauss transform vector

The swapping of the two A-rows may require the involvement of two processors in the network

The appropriate part of the Gauss vector together with (3.6.6) $A _ { k k } ( j , j { : } r )$ is sent by the pivot processor to the processors that house $A _ { k + 1 , k } , \dotsc , A _ { N , k }$

The processors that house $A _ { k k } , \ldots , A _ { N , k }$ carry out their share of the update, a local computation

# end

Upon completion, the parallel execution of Parts B and C follow. In the Part B computation, those blocks that may be involved in the row swapping have been highlighted. See Figure 3.6.3. This overhead generally engages the entire processor network, although communication is local to each processor column.

# Part B:

<table><tr><td>Proc(0,0)</td><td>Proc(0,1)</td><td>Proc(0,0)</td><td>Proc(0,1)</td><td>Proc(0,0)</td><td>Proc(0,1)</td><td>Proc(0,0)</td><td>Proc(0,1)</td></tr><tr><td>U11</td><td>U12</td><td>U13</td><td>U14</td><td>U15</td><td>U16</td><td>U17</td><td>U18</td></tr><tr><td>Proc(1,0)</td><td>Proc(1,1)</td><td>Proc(1,0)</td><td>Proc(1,1)</td><td>Proc(1,0)</td><td>Proc(1,1)</td><td>Proc(1,0)</td><td>Proc(1,1)</td></tr><tr><td>L21</td><td>U22</td><td>U23</td><td>U24</td><td>U25</td><td>U26</td><td>U27</td><td>U28</td></tr><tr><td>Proc(0,0)</td><td>Proc(0,1)</td><td>Proc(0,0)</td><td>Proc(0,1)</td><td>Proc(0,0)</td><td>Proc(0,1)</td><td>Proc(0,0)</td><td>Proc(0,1)</td></tr><tr><td>L31</td><td>L32</td><td>U33</td><td>A34</td><td>A35</td><td>A36</td><td>A37</td><td>A38</td></tr><tr><td>Proc(1,0)</td><td>Proc(1,1)</td><td>Proc(1,0)</td><td>Proc(1,1)</td><td>Proc(1,0)</td><td>Proc(1,1)</td><td>Proc(1,0)</td><td>Proc(1,1)</td></tr><tr><td>L41</td><td>L42</td><td>L43</td><td>A44</td><td>A45</td><td>A46</td><td>A47</td><td>A48</td></tr><tr><td>Proc(0,0)</td><td>Proc(0,1)</td><td>Proc(0,0)</td><td>Proc(0,1)</td><td>Proc(0,0)</td><td>Proc(0,1)</td><td>Proc(0,0)</td><td>Proc(0,1)</td></tr><tr><td>L51</td><td>L52</td><td>L53</td><td>A54</td><td>A55</td><td>A56</td><td>A57</td><td>A58</td></tr><tr><td>Proc(1,0)</td><td>Proc(1,1)</td><td>Proc(1,0)</td><td>Proc(1,1)</td><td>Proc(1,0)</td><td>Proc(1,1)</td><td>Proc(1,0)</td><td>Proc(1,1)</td></tr><tr><td>L61</td><td>L62</td><td>L63</td><td>A64</td><td>A65</td><td>A66</td><td>A67</td><td>A68</td></tr><tr><td>Proc(0,0)</td><td>Proc(0,1)</td><td>Proc(0,0)</td><td>Proc(0,1)</td><td>Proc(0,0)</td><td>Proc(0,1)</td><td>Proc(0,0)</td><td>Proc(0,1)</td></tr><tr><td>L71</td><td>L72</td><td>L73</td><td>A74</td><td>A75</td><td>A76</td><td>A77</td><td>A78</td></tr><tr><td>Proc(1,0)</td><td>Proc(1,1)</td><td>Proc(1,0)</td><td>Proc(1,1)</td><td>Proc(1,0)</td><td>Proc(1,1)</td><td>Proc(1,0)</td><td>Proc(1,1)</td></tr><tr><td>L81</td><td>L82</td><td>L83</td><td>A84</td><td>A85</td><td>A86</td><td>A87</td><td>A88</td></tr></table>

Figure 3.6.3.

Note that Part C involves just a single processor row while the “big” level-three update that follows typically involves the entire processor network. See Figures 3.6.4 and 3.6.5.

Part C: 

<table><tr><td>Proc(0,0)</td><td>Proc(0,1)</td><td>Proc(0,0)</td><td>Proc(0,1)</td><td>Proc(0,0)</td><td>Proc(0,1)</td><td>Proc(0,0)</td><td>Proc(0,1)</td></tr><tr><td>U11</td><td>U12</td><td>U13</td><td>U14</td><td>U15</td><td>U16</td><td>U17</td><td>U18</td></tr><tr><td>Proc(1,0)</td><td>Proc(1,1)</td><td>Proc(1,0)</td><td>Proc(1,1)</td><td>Proc(1,0)</td><td>Proc(1,1)</td><td>Proc(1,0)</td><td>Proc(1,1)</td></tr><tr><td>L21</td><td>U22</td><td>U23</td><td>U24</td><td>U25</td><td>U26</td><td>U27</td><td>U28</td></tr><tr><td>Proc(0,0)</td><td>Proc(0,1)</td><td>Proc(0,0)</td><td>Proc(0,1)</td><td>Proc(0,0)</td><td>Proc(0,1)</td><td>Proc(0,0)</td><td>Proc(0,1)</td></tr><tr><td>L31</td><td>L32</td><td>U33</td><td>A34</td><td>A35</td><td>A36</td><td>A37</td><td>A38</td></tr><tr><td>Proc(1,0)</td><td>Proc(1,1)</td><td>Proc(1,0)</td><td>Proc(1,1)</td><td>Proc(1,0)</td><td>Proc(1,1)</td><td>Proc(1,0)</td><td>Proc(1,1)</td></tr><tr><td>L41</td><td>L42</td><td>L43</td><td>A44</td><td>A45</td><td>A46</td><td>A47</td><td>A48</td></tr><tr><td>Proc(0,0)</td><td>Proc(0,1)</td><td>Proc(0,0)</td><td>Proc(0,1)</td><td>Proc(0,0)</td><td>Proc(0,1)</td><td>Proc(0,0)</td><td>Proc(0,1)</td></tr><tr><td>L51</td><td>L52</td><td>L53</td><td>A54</td><td>A55</td><td>A56</td><td>A57</td><td>A58</td></tr><tr><td>Proc(1,0)</td><td>Proc(1,1)</td><td>Proc(1,0)</td><td>Proc(1,1)</td><td>Proc(1,0)</td><td>Proc(1,1)</td><td>Proc(1,0)</td><td>Proc(1,1)</td></tr><tr><td>L61</td><td>L62</td><td>L63</td><td>A64</td><td>A65</td><td>A66</td><td>A67</td><td>A68</td></tr><tr><td>Proc(0,0)</td><td>Proc(0,1)</td><td>Proc(0,0)</td><td>Proc(0,1)</td><td>Proc(0,0)</td><td>Proc(0,1)</td><td>Proc(0,0)</td><td>Proc(0,1)</td></tr><tr><td>L71</td><td>L72</td><td>L73</td><td>A74</td><td>A75</td><td>A76</td><td>A77</td><td>A78</td></tr><tr><td>Proc(1,0)</td><td>Proc(1,1)</td><td>Proc(1,0)</td><td>Proc(1,1)</td><td>Proc(1,0)</td><td>Proc(1,1)</td><td>Proc(1,0)</td><td>Proc(1,1)</td></tr><tr><td>L81</td><td>L82</td><td>L83</td><td>A84</td><td>A85</td><td>A86</td><td>A87</td><td>A88</td></tr></table>

Figure 3.6.4.

Part D: 

<table><tr><td>Proc(0,0)</td><td>Proc(0,1)</td><td>Proc(0,0)</td><td>Proc(0,1)</td><td>Proc(0,0)</td><td>Proc(0,1)</td><td>Proc(0,0)</td><td>Proc(0,1)</td></tr><tr><td>U11</td><td>U12</td><td>U13</td><td>U14</td><td>U15</td><td>U16</td><td>U17</td><td>U18</td></tr><tr><td>Proc(1,0)</td><td>Proc(1,1)</td><td>Proc(1,0)</td><td>Proc(1,1)</td><td>Proc(1,0)</td><td>Proc(1,1)</td><td>Proc(1,0)</td><td>Proc(1,1)</td></tr><tr><td>L21</td><td>U22</td><td>U23</td><td>U24</td><td>U25</td><td>U26</td><td>U27</td><td>U28</td></tr><tr><td>Proc(0,0)</td><td>Proc(0,1)</td><td>Proc(0,0)</td><td>Proc(0,1)</td><td>Proc(0,0)</td><td>Proc(0,1)</td><td>Proc(0,0)</td><td>Proc(0,1)</td></tr><tr><td>L31</td><td>L32</td><td>U33</td><td>A34</td><td>A35</td><td>A36</td><td>A37</td><td>A38</td></tr><tr><td>Proc(1,0)</td><td>Proc(1,1)</td><td>Proc(1,0)</td><td>Proc(1,1)</td><td>Proc(1,0)</td><td>Proc(1,1)</td><td>Proc(1,0)</td><td>Proc(1,1)</td></tr><tr><td>L41</td><td>L42</td><td>L43</td><td>A44</td><td>A45</td><td>A46</td><td>A47</td><td>A48</td></tr><tr><td>Proc(0,0)</td><td>Proc(0,1)</td><td>Proc(0,0)</td><td>Proc(0,1)</td><td>Proc(0,0)</td><td>Proc(0,1)</td><td>Proc(0,0)</td><td>Proc(0,1)</td></tr><tr><td>L51</td><td>L52</td><td>L53</td><td>A54</td><td>A55</td><td>A56</td><td>A57</td><td>A58</td></tr><tr><td>Proc(1,0)</td><td>Proc(1,1)</td><td>Proc(1,0)</td><td>Proc(1,1)</td><td>Proc(1,0)</td><td>Proc(1,1)</td><td>Proc(1,0)</td><td>Proc(1,1)</td></tr><tr><td>L61</td><td>L62</td><td>L63</td><td>A64</td><td>A65</td><td>A66</td><td>A67</td><td>A68</td></tr><tr><td>Proc(0,0)</td><td>Proc(0,1)</td><td>Proc(0,0)</td><td>Proc(0,1)</td><td>Proc(0,0)</td><td>Proc(0,1)</td><td>Proc(0,0)</td><td>Proc(0,1)</td></tr><tr><td>L71</td><td>L72</td><td>L73</td><td>A74</td><td>A75</td><td>A76</td><td>A77</td><td>A78</td></tr><tr><td>Proc(1,0)</td><td>Proc(1,1)</td><td>Proc(1,0)</td><td>Proc(1,1)</td><td>Proc(1,0)</td><td>Proc(1,1)</td><td>Proc(1,0)</td><td>Proc(1,1)</td></tr><tr><td>L81</td><td>L82</td><td>L83</td><td>A84</td><td>A85</td><td>A86</td><td>A87</td><td>A88</td></tr></table>

Figure 3.6.5.

The communication overhead associated with Part D is masked by the matrix multiplications that are performed on each processor.

This completes the k = 3 step of parallel block LU with partial pivoting. The process can obviously be repeated on the trailing 5-by-5 block matrix. The virtues of the block-cyclic distribution are revealed through the schematics. In particular, the dominating level-3 step (Part D) is load balanced for all but the last few values of k. Subsets of the processor grid are used for the “smaller,” level-2 portions of the computation.

We shall not attempt to predict the fraction of time that is devoted to these computations or the propagation of the interchange permutations. Enlightenment in this direction requires benchmarking.

# 3.6.3 Tournament Pivoting

The decomposition via partial pivoting in Step A requires a lot of communication. An alternative that addresses this issue involves a strategy called tournament pivoting. Here is the main idea. Suppose we want to compute $P W = L U$ where the blocks of

$$
W = \left[ \begin{array}{l} W _ {1} \\ W _ {2} \\ W _ {3} \\ W _ {4} \end{array} \right] \in \mathbb {R} ^ {n \times r}
$$

are distributed around some network of processors. Assume that each $W _ { i }$ has many more rows than columns. The goal is to choose r rows from W that can serve as pivot rows. If we compute the “local” factorizations

$$
P _ {1} W _ {1} = L _ {1} U _ {1}, \qquad P _ {2} W _ {2} = L _ {2} U _ {2}, \qquad P _ {3} W _ {3} = L _ {3} U _ {3}, \qquad P _ {4} W _ {4} = L _ {4} U _ {4},
$$

via Gaussian elimination with partial pivoting, then the top r rows of the matrices $P _ { 1 } W _ { 1 } , P _ { 2 } W _ { 2 } , P _ { 3 } W _ { 3 }$ , are $P _ { 4 } W _ { 4 }$ are pivot row candidates. Call these square matrices $W _ { 1 } ^ { \prime } , W _ { 2 } ^ { \prime } , W _ { 3 } ^ { \prime } ,$ , and $W _ { 4 } ^ { \prime }$ and note that we have reduced the number of possible pivot rows from n to 4r.

Next we compute the factorizations

$$
P _ {1 2} W _ {1 2} ^ {\prime} = P _ {1 2} \left[ \begin{array}{c} W _ {1} ^ {\prime} \\ W _ {2} ^ {\prime} \end{array} \right] = L _ {1 2} U _ {1 2},
$$

$$
P _ {3 4} W _ {3 4} ^ {\prime} = P _ {3 4} \left[ \begin{array}{c} W _ {3} ^ {\prime} \\ W _ {4} ^ {\prime} \end{array} \right] = L _ {3 4} U _ {3 4},
$$

and recognize that the top r rows of $P _ { 1 2 } W _ { 1 2 } ^ { \prime }$ and the top r rows of $P _ { 3 4 } W _ { 3 4 } ^ { \prime }$ are even better pivot row candidates. Assemble these 2r rows into a matrix $W _ { 1 2 3 4 }$ and compute

$$
P _ {1 2 3 4} W _ {1 2 3 4} = L _ {1 2 3 4} U _ {1 2 3 4}.
$$

The top r rows of $P _ { 1 2 3 4 } W _ { 1 2 3 4 }$ are then the chosen pivot rows for the LU reduction of $W$ .

Of course, there are communication overheads associated with each round of the “tournament,” but the volume of interprocessor data transfers is much reduced. See Demmel, Grigori, and Xiang (2010).

# Problems

P3.6.1 In §3.6.1 we outlined a single step of block LU with partial pivoting. Specify a complete version of the algorithm.

P3.6.2 Regarding parallel block LU with partial pivoting, why is it better to “collect” all the permutations in Part A before applying them across the remaining block columns? In other words, why not propagate the Part A permutations as they are produced instead of having Part B, a separate permutation application step?

P3.6.3 Review the discussion about parallel shared memory computing in §1.6.5 and §1.6.6. Develop a shared memory version of Algorithm 3.2.1. Designate one processor for computation of the multipliers and a load-balanced scheme for the rank-1 update in which all the processors participate. A barrier is necessary because the rank-1 update cannot proceed until the multipliers are available. What if partial pivoting is incorporated?

# Notes and References for §3.6

See the scaLAPACK manual for a discussion of parallel Gaussian elimination as well as:

J. Ortega (1988). Introduction to Parallel and Vector Solution of Linear Systems, Plenum Press, New York.   
K. Gallivan, W. Jalby, U. Meier, and A.H. Sameh (1988). “Impact of Hierarchical Memory Systems on Linear Algebra Algorithm Design,” Int. J. Supercomput. Applic. 2, 12–48.   
J. Dongarra, I. Duff, D. Sorensen, and H. van der Vorst (1990). Solving Linear Systems on Vector and Shared Memory Computers, SIAM Publications, Philadelphia, PA.   
Y. Robert (1990). The Impact of Vector and Parallel Architectures on the Gaussian Elimination Algorithm, Halsted Press, New York.   
J. Choi, J.J. Dongarra, L.S. Osttrouchov, A.P. Petitet, D.W. Walker, and R.C. Whaley (1996). “Design and Implementation of the ScaLAPACK LU, QR, and Cholesky Factorization Routines,” Scientific Programming, 5, 173–184.   
X.S. Li (2005). “An Overview of SuperLU: Algorithms, Implementation, and User Interface,” ACM Trans. Math. Softw. 31, 302–325.   
S. Tomov, J. Dongarra, and M. Baboulin (2010). “Towards Dense Linear Algebra for Hybrid GPU Accelerated Manycore Systems,” Parallel Comput. 36, 232–240.

The tournament pivoting strategy is a central feature of the optimized LU implementation discussed in:

J. Demmel, L. Grigori, and H. Xiang (2011). “CALU: A Communication Optimal LU Factorization Algorithm,” SIAM J. Matrix Anal. Applic. 32, 1317-1350.

E. Solomonik and J. Demmel (2011). “Communication-Optimal Parallel 2.5D Matrix Multiplication and LU Factorization Algorithms,” Euro-Par 2011 Parallel Processing Lecture Notes in Computer Science, 2011, Volume 6853/2011, 90–109.

This page intentionally left blank
