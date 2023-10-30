# Most Frugal Explanations

LibDAI-based C++ implementation of the Most Frugal Explanation heuristic as well as the Annealed MAP approximation algorithm, together with a tool that translates .bif Bayesian networks into .fg factor graphs that can be read by LibDAI and an implementation of the MAP independence algorithm. References:

J. Kwisthout (2015). Most Frugal Explanations in Bayesian Networks. Artificial Intelligence, 218, 56 - 73. 

J. Kwisthout (2021). Explainable AI using MAP-independence. Proceedings of the 16th European Conference on Symbolic and Quantitative Approaches to Reasoning with Uncertainty (ECSQARU'21). Springer LNCS, vol 12897.

C. Yuan, T. Lu, and M. J. Druzdzel. Annealed MAP. In D. Chickering and J. Halpern,
editors, Proceedings of the Twentieth Conference in Uncertainty in Articial Intelligence,
pages 628{635. AUA, 2004.

Correct compilation requires a compiled libdai library. These are not included here, but libdai.a or libdai32.a can be requested from the author if you cannot build those libraries yourself.

Make sure to make directories /object and /release before compiling. Depending on your version of the STL, you may need to change <experimental/random> in <random> to avoid compilation errors (like "'mt19937' is not a member of 'std'").
