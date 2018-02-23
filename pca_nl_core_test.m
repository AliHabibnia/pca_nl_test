function [summ, bounds, M] = pca_nl_core_test(Z, N, K, opt, pc_n)

%{
A Nonlinearity Test for Principal Component Analysis (V 1.0)

  [summ] = pca_nl_core_test(Z, N, K, opt, pc_n) returns
  Residual variances (discarded eigenvalues) for all the disjunct regions.

  A test developed by Kruger et al (2005) to determine whether
  the underlying structure within the recorded data is linear or nonlinear.

  N series of data each with the total of K observations are broken down
  to several regions and stored in Z .
  The correlation/covariance matrices for each disjunct region is calculated,
  carring out singular value decomposition in order to find out the sum of
  discarded eigenvalues by leaving out the first pc_n.
  These values are compared against the accuracy bounds for the discarded eigenvalues
  of the first region using the opt as the optimization algorithm.

  [summ, bounds] = pca_nl_core_test(...) returns the accuracy bounds
  for the first region.

  [summ, bounds, M] = pca_nl_core_test(...) returns the number of disjunct regions.

  % Example:
  %   Generate two series with a nonlinear underlying relation
  t = 8.*rand(1000, 1)-4;
  e1 = normrnd(0, sqrt(0.005), [1000,1]);
  e2 = normrnd(0, sqrt(0.005), [1000,1]);
  z1_l = (t+e1)';
  z2_nl = ((t.^3)+e2)';

  % Divide the series into 4 disjunct regions
  dataset = [z1_l; z2_nl];
  [~, InX] = sort(abs(dataset(1,:)));
  dataset = dataset(:,InX);
  Z = mat2cell(dataset, 2, diff(round(linspace(0, 1000, 4+1))));

  % Scale and center all regions based on the first one
  [~, std_info] = mapstd(Z{1});
  for v = 1:1:4
  Z{v} = mapstd('apply', Z{v}, std_info);
  end

  % Compute the nonlinearity test
  [summ, bounds, M] = pca_nl_core_test(Z, 2, 1000, 'greedy', 1);

  See also pca_nl_test.

Reference papers:

* Kruger, U., Antory, D., Hahn, J., Irwin, G. and McCullough, G. (2005). Introduction of a nonlinearity measure for principal component models. Computers & Chemical Engineering, 29(11-12), pp.2355-2362.
* Kruger U., Zhang J., Xie L. (2008) Developments and Applications of Nonlinear Principal Component Analysis – a Review. In: Gorban A.N., Kégl B., Wunsch D.C., Zinovyev A.Y. (eds) Principal Manifolds for Data Visualization and Dimension Reduction. Lecture Notes in Computational Science and Enginee, vol 58. Springer, Berlin, Heidelberglications of nonlinear principal component analysis – a review

  Copyright 2018
  Cite: Habibnia, A., Rahimikia, E., & Mahdikhah, H. (2018). A Nonlinearity Test for Principal Component Analysis, MATLAB Central File Exchange. Retrieved Feb, 2018.
  Last Revision: Feb-2018 (Version 1.0)
  Feel free to send your feedback.
%}

%%
M = numel(Z); % Number of regions

%% Calculating the confidence limits and the accuracy bounds for the first disjunct region.
S1zz = corr(Z{1}'); % Correlation matrix of the disjunct region for which the accuracy bounds are to be determined (Rzz).
Eee = sqrt(1/(round(K/M)-3))*1.96; % Epsilon.

%% Upper (S1u) and lower (S1l) confidence limit matrices (Rzz_u and Rzz_l).
S1l = zeros(N);
S1u = zeros(N);
for i=1:N-1
    for j=i+1:N
        
        if abs(S1zz(i,j))<1
            x = log((1+S1zz(i,j))/(1-S1zz(i,j)))/2;
        else
            x = log(1.9/.1)/2;
        end
        
        y1 = x+Eee;
        y2 = x-Eee;
        
        S1u(i,j) = (exp(2*y1)-1)/(exp(2*y1)+1);
        S1u(j,i) = S1u(i,j);
        
        S1l(i,j) = (exp(2*y2)-1)/(exp(2*y2)+1);
        S1l(j,i) = S1l(i,j);
    end
    S1l(i,i) = 1;
    S1u(i,i) = 1;
end

S1l(N,N) = 1;
S1u(N,N) = 1;

Er = S1u - S1l;
Ss = S1l;

%% Using greedy or genetic algorithem to find the accuracy bounds (rho_max and rho_min).
n_o_t = 0;

if ~license('test', 'optimization_toolbox')
    sprintf('The "Optimization Toolbox" is not available. The optimization option is changed to "greedy".');
    n_o_t = 1;
end

if n_o_t == 1 || strcmpi(opt, 'greedy')
    vect = zeros(20000, N);
    for i=1:20000
        vect(i,:)= svd(Ss);
        Ss=Ss+Er/20000;
    end
    
    bounds(1)=min(sum(vect(:,(pc_n+1):end),2));
    bounds(2)=max(sum(vect(:,(pc_n+1):end),2));
else
    
    num_of_var = 1;
    lb = 1;
    up = 20000000;
    
    options = optimoptions('ga');
    options.Display = 'iter';
    
    min_f = @(w) f_1(Ss,w,Er,pc_n);
    max_f = @(w) f_2(Ss,w,Er,pc_n);
    [~, fval_1] = ga(min_f, num_of_var, [], [], [], [], lb, up,[],[], options); bounds(1) = fval_1;
    [~, fval_2] = ga(max_f, num_of_var, [], [], [], [], lb, up,[],[], options); bounds(2) = -fval_2;
    
end

%% Calculating the sum of discarded eigenvalues for all regions.
summ=zeros(M,1);

for i=1:M
    Snzz = cov(Z{i}');
    Sn=svd(Snzz);
    summ(i)=sum( Sn( (pc_n+1):end ,1) ,1);
end
end