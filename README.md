# Testing Independence via Graph Correlation

## Simulation - preliminaries

* [Network(Graph) vs.Node attributes](https://rawgit.com/neurodata/youjin/master/report/network_attribute.html) 
* [Network 1 vs Network 2](http://htmlpreview.github.io/?https://github.com/neurodata/youjin/blob/master/report/network_network.html)
* [Multiple Networks](http://htmlpreview.github.io/?https://github.com/neurodata/youjin/blob/master/report/multinetworks.html)
* [Local Dependence](http://htmlpreview.github.io/?https://github.com/neurodata/youjin/blob/master/report/local_dependence.html) 
* [Optimal Time in Diffusion Distance](http://htmlpreview.github.io/?https://github.com/neurodata/youjin/blob/master/report/local_time.html)

## Simulation - block models
* [Two Blocks(B2.1)](http://htmlpreview.github.io/?https://github.com/neurodata/youjin/blob/master/report/twoblocks.html)
* [Three Blocks(B3.1)](http://htmlpreview.github.io/?https://github.com/neurodata/youjin/blob/master/report/B3_1.html)
* [Three Blocks(B3.2)](https://rawgit.com/neurodata/youjin/master/report/B3_2.html)
* [Three Blocks(B3.3)](http://htmlpreview.github.io/?https://github.com/neurodata/youjin/blob/master/report/threeblocks.html)
* [Three Blocks(B3.4)](http://htmlpreview.github.io/?https://github.com/neurodata/youjin/blob/master/report/B3_4.html)
* [Comments](http://htmlpreview.github.io/?https://github.com/neurodata/youjin/blob/master/report/comments.html)

## Simulation - hierarchical latent models 
* [preliminaries](http://htmlpreview.github.io/?https://github.com/neurodata/youjin/blob/master/report/statistics.html)
* [Two Blocks](https://rawgit.com/neurodata/youjin/master/report/latent_two.html)
* [Three Blocks](https://rawgit.com/neurodata/youjin/master/report/latent_three.html)

## Simulation - distance metrics
* [diffusion maps](https://rawgit.com/neurodata/youjin/master/report/diffusion.html)

## Real Data Analysis
* [C.elegans Neuronal Network](http://htmlpreview.github.io/?https://github.com/neurodata/youjin/blob/master/report/RealData.html)


## Comparisons 
* [Random Dot Product Graph](https://rawgit.com/neurodata/youjin/master/report/RDPG.html)
* [Degree-Corrected Stochastic Block Model](https://rawgit.com/neurodata/youjin/master/report/DCSBM.html)





## Reference
* celegans : exploratory data analysis on neuronal network of C. elegans.

  - http://www.jstor.org/stable/pdf/2528823.pdf?_=1460126515800 - distance measures for mixed data types

*  Diffusion Distance : find an appropriate distance for network(graph) between two extremes - adjacency matrix vs. geodesic distance.

 - http://www.ams.jhu.edu/~priebe/.FILES/MinhTang-jhu_nov11.pdf
 - using graph diffusion distance for comparing two graphs : http://ieeexplore.ieee.org/xpl/login.jsp?tp=&arnumber=6736904&url=http%3A%2F%2Fieeexplore.ieee.org%2Fxpls%2Fabs_all.jsp%3Farnumber%3D6736904
 
* Local Dependence : assumes that each edge variable depends on a finite subset of other edge variables, e.g. finite neighborhoods, M-dependence
 - http://www.stat.rice.edu/~ms88/publications/h.ergm.pdf
