## Distributed consent

This repository contains the codes associated with the article

Limits of individual consent and models of distributed consent in online social networks<br/>
[J. Lovato], [A. Allard], [R. Harp], [J. Onaolapo] and [L. Hébert-Dufresne]<br/>
[arXiv:2006.16140]


#### Facebook100 dataset

The [Facebook100 dataset](http://doi.org/10.1016/j.physa.2011.12.021) must be [downloaded](https://archive.org/details/oxford-2005-facebook-matrix) and the `.mat` files must be be placed in the `Facebook100` directory. The individual edgelists are then extracted by executing the python script `extract_edgelists.py`.


#### Numerical simulations and generating the figures

Figures 2 and 3 were generated by running (from the respective folders)

1. the numerical simulation by executing the bash script `run_simulations.sh`
2. the python script `plot_figure.py`



[arXiv:2006.16140]: https://arxiv.org/abs/2006.16140
[J. Lovato]: http://juniperlovato.com/
[A. Allard]: http://antoineallard.info
[R. Harp]: http://www.uvm.edu/~rharp/
[J. Onaolapo]: https://www.uvm.edu/~jonaolap/
[L. Hébert-Dufresne]: http://laurenthebertdufresne.github.io/
