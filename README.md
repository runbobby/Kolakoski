# Kolakoski
Various c++ programs for computation with the Kolakoski sequence.

Includes some joint work with Yongyi Chen and Michael Yan from Spring 2016. In particular. correlations.cpp is work from 2017.

See the following papers, which discuss the Kolakoski sequence and various related algorithms. Admittedly, the code is still lacking in documentation.

https://www.overleaf.com/4401114cgstjy#/13131479/

arXiv:1702.08156
 
arXiv:1703.00180 (maybe the url will be updated; not sure)
 
The procedues are largely just straightforward recursion, but there are a few tricky elements. These are cycle decomposition of a permutation in order to exponentiate permutations efficiently, constructing Eulerian cycles from 2-regular graphs, and FFT multiplication.

The programs were executed on a Latitude E7450, 64-bit, 8 GB Ram, 2.60 GHz, Visual Studio 2015. For some reason, VS15's diagnostics claim that "correlations.cpp" uses 13 GB of RAM. I honestly have no idea what that means. My estimate was closer to 6.5 GB.

------------------------------------------------------------------------

Brief description of files, as I remember.

Startsdist.cpp:  So there is one big commented block in main(). The purpose of this is to informally test the Generalized Uniformness Conjecture (GUC) (see the arXiv papers) using a chi^2 test. Also in this program, we look for sequences whose orbit length under the functions C\_{1,n}(1, -)  or C\_{-1,n}(-1,-) match the "easy upperbound"; in this sense, they have the maximum possible orbit lengths. This is related to the rather funny conjectures 3.4 and 3.7 in arXiv 1702.08156.

OrbitChecker.cpp: umm probably badly named. Purpose is to verify Conjecture 3.4 for n = 2, j<=23.

KolakoskiMNOrbit.cpp: Purpose is to verify Conjecture 3.4 for 0 < N < 4096, j <= 13. In arXiv:1702.08156, we proved that this implies the conjecture for j<=13 and all N

KolakoskiMNLength.cpp: I think this was a temporary program for debugging; the main stuff was moved to GraphAlgorithm.cpp

GraphAlgorithm.cpp: The main purpose is to compmute bounds on the (supremum and infimum) density of 1 in the Kolakoski sequence. See section 6 of the Overleaf paper. It is an old open problem to prove/disprove that this density equals 1/2. We used a Bellman-Ford-type algorithm and heuristically guess that negative cycles are found after 20(k) iterations. This is obviously very optimistic. However, the guess can be verified, because if there are no negative cycles, then the algorithm will stop updating weights, and if this happens after 20(k) iterations, that's fast. I must admit that there's probably some bugs for k >= 26.

correlations.cpp: This program computes the correlation function tf(m,n,d) for small pairs (m,n), and heuristically, it does so near exactly for very large d: d ~ 5000 for m+n = 3, d ~ 150000 for m+n = 5, and d ~ 700000 for m+n = 7. The GUC essentially states that the sequence K(m,n) induces a walk on G\_{m,n,k} that looks random; in particular, it uses each edge asymptotically the same proportion of the time. Well, an Eulerian cycle also does this. We prove in arXiv:1703.00180 that using the Eulerian cycle gives the correct correlation frequencies for rather large d.

Of course, we push the limits and decide to use an Eulerian cycle that induces a periodic sequence with period ~ 10^8.5 How do we compunte correlations in a sequence this long for d = 1, 2, ..., 10^6? Of course, we could do one loop for each d through the entire sequence. This would be mathematically correct and unfortunately too slow for my machine. Fortunately, the correlations can be computed concurrently by transforming the problem into polynomial multiplication, and this can be done using FFT. The output goes into a text file (separated by commas) so that it can be imported into an interactive analytical environment. In particular, I used mathematica.
