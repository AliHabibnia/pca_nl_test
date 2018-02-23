%{ 
A Nonlinearity Test for Principal Component Analysis (V 1.0)

A test developed by Kruger et al. (2005) to determine whether the underlying structure within
the recorded data is linear or nonlinear. 

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

See also pca_nl_test and pca_nl_core_test.

  Copyright 2018
  Cite: Habibnia, A., Rahimikia, E., & Mahdikhah, H. (2018). A Nonlinearity Test for Principal Component Analysis, MATLAB Central File Exchange. Retrieved Feb, 2018.
  Last Revision: Feb-2018 (Version 1.0)
  Feel free to send your feedback.
%}

%%
clc
clear
close all

t = 8.*rand(1002, 1)-4;
e1 = normrnd(0, sqrt(0.005), [1002,1]);
e2 = normrnd(0, sqrt(0.005), [1002,1]);
e3 = normrnd(0, sqrt(0.005), [1002,1]);

%% Linear transformations of t
z1_l = (t+e1)';
z2_l = (t+e2)';
z3_l = (t+e3)';

%% Nonlinear transformations of t
z1_nl = ((t.^3)+e1)';
z2_nl = ((t.^3)+e2)';
z3_nl = ((t.^3)+e3)';

%% Paper Examples
x = 6.*rand(3000, 1)-3;
e1 = normrnd(0, 1, [3000,1]);
e2 = normrnd(0, 1, [3000,1]);
z1_p_l = (x+(0.05*e1))';
z2_p_l = (x+(0.05*e2))';
z2_p_nl = (sin(x)+(0.05*e2))';

%% 
%{
Brief review of options

'optimization' : 'greedy' / 'genetic'
>> default: 'greedy'
Toolbox requirement: Global Optimization Toolbox for 'genetic'.

'region_selection' : 'auto' / 'manual' / 'uniform' -- 'uniform'.
>> default: {'uniform', 4}
>>>> auto default: [3, 6]
>>>> uniform default: 4
>>>> manual default: *NAN*

'series_fig' : 'on' / 'off'
>> default: 'off'

'scatter_fig' : 'on' / 'off'
>> default: 'off'

'output_fig' : 'on' / 'off'
>> default: 'on'

'wait_bar' : 'on' / 'off'
>> default: 'on'

'print_results' : 'on' / 'off'
>> default: 'on'

'PC_number' : Number of principal components.
>> default: 1

'single_figure' : Open all figures in a single window (valid for 2-8 figures).
>> default: 'off'
%}

%%%%%%%%%%%%%%%%%%%%%%%%% Examples %%%%%%%%%%%%%%%%%%%%%%%%%

%% Example 1: Two series with linear interrelationship.
%{
the "linear example" in Ref: Kruger U., Zhang J., Xie L. (2008) Developments and Applications of Nonlinear Principal Component Analysis – a Review. In: Gorban A.N., Kégl B., Wunsch D.C., Zinovyev A.Y. (eds) Principal Manifolds for Data Visualization and Dimension Reduction. Lecture Notes in Computational Science and Enginee, vol 58. Springer, Berlin, Heidelberglications of nonlinear principal component analysis – a review

The example generates two series with linear relation breaking them down to 4 even regions (using the
default value of the 'region_selection' which is {'uniform', 4}).

This will result in a Scatter diagram of the two series and four plots of benchmarking the residual
variances against accuracy bounds of each disjunct region. These plots and also the 'InRegion_per' in
the printed results (which shows the percentage of regions falling inside the accuracy bounds and
can be turned off using the 'print_results' option) yield that no violation of the accuracy bounds arise,
which leads to the acceptance of the hypothesis that the underlying relationship between the two series is linear.
%}

[tbx1_e1, tbx2_e1, opt_e1] = pca_nl_test([z1_l; z2_l], 'scatter_fig', 'on', 'single_fig', 'on');


%% Example 2: Two series with nonlinear interrelationship.
%{
the "nonlinear example" in Ref: Kruger U., Zhang J., Xie L. (2008) Developments and Applications of Nonlinear Principal Component Analysis – a Review. In: Gorban A.N., Kégl B., Wunsch D.C., Zinovyev A.Y. (eds) Principal Manifolds for Data Visualization and Dimension Reduction. Lecture Notes in Computational Science and Enginee, vol 58. Springer, Berlin, Heidelberglications of nonlinear principal component analysis – a review

The example generates two series with nonlinear relation breaking them down to 4 even regions.

This will result in a Scatter diagram of the two series and a single figure containing four plots of
benchmarking the residual variances against accuracy bounds of each disjunct region. These plots and also
the 'InRegion_per' in the printed results yield that in contrast to the linear example, the residual variance
of the reconstructed data for each of the disjunct regions other than the benchmark, fall outside the accuracy bounds,
leading to the rejection of the hypothesis that the underlying relationship between the two series is linear.
%}

[tbx1_e2, tbx2_e2, opt_e2] = pca_nl_test([z1_l; z2_nl], 'scatter_fig', 'on', 'single_figure', 'on');



%% Example 3: Three nonlinear transformations with linear interrelationship.
%{
In this example three nonlinear transformations are considered which despite each of them being a nonlinear
transformation of the series t, have a linear interrelationship amongst themselves. This can be visually shown
using 'scatter_fig' option which results in a 3D scatter plot showing the linear relation between these three series.
for this example we use 'auto' option of region selection, dividing series into 3 to 5 regions. Also, since there are more
than two series involved, we are going to use the genetic algorithm as the optimization option which leads to a more
accurate calculation of accuracy bounds in a faster fashion.

The outputs here, just like the first example, are going to show that the null hypothesis of
linear underlying relationship between the series cannot be rejected.
%}

[tbx1_e3, tbx2_e3, opt_e3] = pca_nl_test([z1_nl; z2_nl; z3_nl], 'scatter_fig', 'on', 'single_figure', 'on',...
    'region_selection', {'auto', [3,5]}, 'optimization', 'genetic');



%% Example 4: Two series with nonlinear interrelationship (sin transformation).
%{
the "nonlinear example" in Ref: Kruger, U., Antory, D., Hahn, J., Irwin, G. and McCullough, G. (2005). Introduction of a nonlinearity measure for principal component models. Computers & Chemical Engineering, 29(11-12), pp.2355-2362.

Same as the second example, only using sin as the nonlinear transformation and dividing the series into 3 regions.
The results will be same as the second example.
%}

[tbx1_e4, tbx2_e4, opt_e4] = pca_nl_test([z1_p_l; z2_p_nl], 'scatter_fig', 'on', 'single_figure', 'on', 'region_selection', {'uniform', 3});

