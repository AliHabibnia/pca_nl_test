# pca_nl_test
A Nonlinearity Test for Principal Component Analysis

** Attached MATLAB code is developed to test whether the underlying structure within the recorded data is linear or nonlinear. The nonlinearity measure introduced in Kruger et al (2005) performs a multivariate analysis assessing the underlying relationship within a given variable set by dividing the data series into smaller regions, calculating the sum of the discarded eigenvalues and the accuracy bounds of them for each region and then benchmarking them against each other. This can help researchers to better define the nature of the relationship among their data series and to utilize linear and nonlinear PCA methods in a more effective fashion.

Developers:
               Dr. Ali Habibnia
               London School of Economics
               alihabibnia@gmail.com

               Eghbal Rahimikia
               Iran University of Science and Technology
               erahimi@ut.ac.ir

               Hossein Mahdikhah
               University of Tehran
               hmkh@outlook.com

Files:

* Main function: pca_nl_test.m
* Core (related) function: pca_nl_core_test.m
* Demo: Demo.m
* Otimization functions: f_1.m and f_2.m
* README: readme (this file)


Reference papers:

* Kruger, U., Antory, D., Hahn, J., Irwin, G. and McCullough, G. (2005). Introduction of a nonlinearity measure for principal component models. Computers & Chemical Engineering, 29(11-12), pp.2355-2362.
* Kruger U., Zhang J., Xie L. (2008) Developments and Applications of Nonlinear Principal Component Analysis – a Review. In: Gorban A.N., Kégl B., Wunsch D.C., Zinovyev A.Y. (eds) Principal Manifolds for Data Visualization and Dimension Reduction. Lecture Notes in Computational Science and Enginee, vol 58. Springer, Berlin, Heidelberglications of nonlinear principal component analysis – a review

** Copyright 2018
** Cite: Habibnia, A., Rahimikia, E., & Mahdikhah, H. (2018). A Nonlinearity Test for Principal Component Analysis, MATLAB Central File Exchange. Retrieved Feb, 2018.
** Last Revision: Feb-2018 (Version 1.0)
