# Rank Aggregation Problem

This code was developed as part of the course ***Design and Analysis of Algorithms*** of the Graduate Program in Computer Science of the Federal University of Minas Gerais.

## The work

In this work is presented some heuristics and approximation algorithms to solve the Rank Aggregation Problem based on [Kameny Young](https://en.wikipedia.org/wiki/Kemeny%E2%80%93Young_method) method. ***Rank aggregation problem***
can be defined as: *"The problem of computing a consensus ranking of the alternatives, given the individual ranking preferences of several judges"[[1](http://www10.org/cdrom/papers/pdf/p577.pdf)]*.

## Brief motivation
Rank aggragation has been extensively studied in the context of social choice theory, where many vote systems were proposed.
The Kemeny-Young method is one of the most classical and natural. The problem of computing Kemeny-Young ranking is NP-hard[[2](http://www.ime.usp.br/~rbrito/docs/voting/BF00303169.pdf)] even if there are 4 voters.

## Heuristics and Approximation Algorithms
The implmented algorithms are the followings:
- Pick_A_List: 2-approximation algorithm
- kPick_A_List: 2-approximation algorithm with locally kamenization(guarantees Condorcet winner)
- FAS_pivot: 2-approximation algorithm
- kFas_pivot: 2-approximation algorithm with locally kamenization(guarantees Condorcet winner)
- FHP_greedy: greedy algorithm
- kFHP_greedy: greedy algorithm with locally kamenization(guarantees Condorcet winner)
- Mixed: meta-heuristic, returns the best result among the results of the above algorithms

For further details consult the file Final_Project_PAA.pdf which is found in this repository(Unfortunately, it is portuguese because I had no time to translate it to English - sorry)

### References
- [[3](http://dimacs.rutgers.edu/~alantha/papers2/acn05conf.pdf)] N.  Ailon,  M.  Charikar,  and  A.  Newman.   Aggregation  inconsistent  information:  Ranking and clustering.  2005.
- [3] A. Bar-Noy and J. Naor. Sorting, minimal feedback sets and hamilton paths in tournaments. 1988.
- [5] C. Dwork, R. Kumar, M. Naor, and D. Sivakumar. Rank Aggregation Methods for the Web. WWW10, 2001.
- [6] C.  Dwork,  R.  Kumar,  M.  Naor,  and  D.  Sivakumar. Rank  Aggregation  Revisited. 10th International World Wide Web Conference, May 2001.
- [7] J. Kleinberg and E. Tardos. Algorithm Design .  Addison Wesley, 2006.
- [8] G. Lv. An Analysis of Rank Aggregation Algorithms.arxiv.org, 2014.
- [9] H. P. Young. Condorcet's theory of voting. AMERICAN POLITICAL SCIENCE REVIEW,84(4), December 1988
