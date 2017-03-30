# Kolakoski
Various c++ algorithms for computation with the Kolakoski sequence.

Includes some joint work with Yongyi Chen and Michael Yan from Spring 2015.

See the following papers, which discuss the Kolakoski sequence and various related algorithms. Admittedly, the code is still lacking in documentation.

https://www.overleaf.com/4401114cgstjy#/13131479/

arXiv:1702.08156
 
arXiv:1703.00180 (maybe the url will be updated; not sure)
 
The procedues are largely just straightforward recursion, but there are a few tricky elements. These are cycle decomposition of a permutation in order to exponentiate permutations efficiently, constructing Eulerian cycles from 2-regular graphs, and FFT multiplication.

The programs were executed on a Latitude E7450, 64-bit, 8 GB Ram, 2.60 GHz, Visual Studio 2015. For some reason, VS15's diagnostics claim that "correlations.cpp" uses 13 GB of RAM. I honestly have no idea what that means. My estimate was closer to 6.5 GB.

------------------------------------------------------------------------

Brief description of files, as I remember.

