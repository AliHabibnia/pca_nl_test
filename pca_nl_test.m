function [tbx_1, tbx_2, options] = pca_nl_test(data, varargin)

%{
A Nonlinearity Test for Principal Component Analysis (V 1.0)
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

Manual
  [tbx_1] = pca_nl_test(data, ...) returns Residual variances (discarded eigenvalues)
  for all the disjunct regions and the accuracy bounds for the benchmark region, considering
  each region as the benchmark once.

  A test developed by Kruger et al (2005) to determine whether
  the underlying structure within the recorded data is linear or nonlinear.

  'data' is a N*K matrix of N series of K observations. These series are broken down
  to several regions. The correlation/covariance matrices for each disjunct region is calculated,
  carrying out singular value decomposition in order to find out the sum of
  discarded eigenvalues by leaving out the first PC_number.
  These values are compared against the accuracy bounds for the discarded eigenvalues
  of the benchmark region using genetic (toolbox requirement: Global Optimization Toolbox) or greedy optimization algorithms.

  [tbx_1, tbx_2] = pca_nl_test(data, ...) returns a table checking whether the value for each region falls
  outside of the accouracy bounds and hence rejecting the null hypothesis of linear structure.

  [tbx_1, tbx_2, options] = pca_nl_test(data, ...) returns a struct containing a list of all the input options.

Input Options
  [tbx_1, tbx_2, options] = pca_nl_test(data, 'optimization', 'greedy', ...) choosing the optimization algorithm.
  the available options are 'greedy' and 'genetic'. default: 'greedy'

  [tbx_1, tbx_2, options] = pca_nl_test(data, 'region_selection', {'uniform', M}, ...) choosing the how data is
  going to be divided into several regions for the tests. The input is a 1*2 cell. The first item of the cell
  defines the approach for the division and it can be 'uniform', 'auto' or 'manual'. Choosing 'uniform' will divide
  data into M equal regions (if K is not divisible by M, observations are left out from the end of series). Default value
  for M in case of 'uniform' is 4. 'auto' will act the same as 'uniform' with the only difference being that M contains the
  minimum and the maximum number of regions. Default value for M in case of 'auto' is [3, 6]. Another option for the region
  selection is 'manual' in which case M provides indices for breaking the data. default: {'uniform', 4}


  [tbx_1, tbx_2, options] = pca_nl_test(data, 'series_fig', 'off', ...) option to turn the plot of the series 'on'
  or 'off'. default: 'off'

  [tbx_1, tbx_2, options] = pca_nl_test(data, 'scatter_fig', 'off', ...) option to turn the scatter plot of the series
  to show the underlying relation of them, 'on' or 'off'. This option is only available if N<4. default: 'off'

  [tbx_1, tbx_2, options] = pca_nl_test(data, 'output_fig', 'on', ...) option to turn the output plots of the regions
  and accuracy bounds 'on' or 'off'. default: 'on'

  [tbx_1, tbx_2, options] = pca_nl_test(data, 'wait_bar', 'on', ...) option to turn the progress bar 'on' or 'off'.
  default: 'on'

  [tbx_1, tbx_2, options] = pca_nl_test(data, 'print_results', 'on', ...) option to turn printing results 'on' or 'off'.
  default: 'on'

  [tbx_1, tbx_2, options] = pca_nl_test(data, 'PC_number', 1, ...) choosing the number of principal components to keep.
  default: 1

  [tbx_1, tbx_2, options] = pca_nl_test(data, 'single_figure', 'off', ...) option to plotting up to 8 plots in a single
  figure. default: 'off'


Example:

Generate two series with a nonlinear underlying relation
  t = 8.*rand(1000, 1)-4;
  e1 = normrnd(0, sqrt(0.005), [1000,1]);
  e2 = normrnd(0, sqrt(0.005), [1000,1]);
  z1_l = (t+e1)';
  z2_nl = ((t.^3)+e2)';

Compute the nonlinearity test
  [tbx1_e1, tbx2_e1, opt_e1] = pca_nl_test([z1_l; z2_l]);

See also pca_nl_core_test.

  Ref: Uwe Kruger, David Antory, Juergen Hahn, George W. Irwin, Geoff McCullough
       		Introduction of a nonlinearity measure for principal component models.
       Uwe Kruger, Junping Zhang, Lei Xie
       		Developments and applications of nonlinear principal component analysis – a review


  Copyright 2018
  Cite: Habibnia, A., Rahimikia, E., & Mahdikhah, H. (2018). A Nonlinearity Test for Principal Component Analysis, MATLAB Central File Exchange. Retrieved Feb, 2018.
  Last Revision: Feb-2018 (Version 1.0)


Version History:
Initial release: 1.0
Feel free to send your feedback.

%}

%% Parser
p = inputParser;
default_on = 'on';
default_off = 'off';
default_region_selection = {'uniform', 4};
default_pc_num = 1;
default_optimization = 'greedy';

    function TF = check_data(x)
        TF = false;
        if ~isa(x, 'double')
            error('The data type must be double.');
        elseif size(x, 1) < 2 || size(x, 2) <= 2
            error('The number of rows/columns of dataset is not sufficient.');
        elseif ~isequal(sum(sum(isnan(x))), 0)
            error('All elements of dataset must be numeric.');
        else
            TF = true;
        end
    end

    function TF = check_pcn(x)
        TF = false;
        if ~isa(x, 'double') || ~isequal(numel(x),1) || x<=0
            error('The number of princial components should be a positive numeric value.');
        elseif x >= size(data, 1)
            error('The number of principal components should be less than the number of the series.');
        else
            TF = true;
        end
    end

    function TF = check_rgs(x)
        TF = false;
        
        if ischar(x) && (~strcmpi(x,'uniform') && ~strcmpi(x,'auto'))
            error('In this format, region selection criterion can be "uniform" or "auto" / Use cell if you have pair inputs for "region_selection".');
            
        elseif ischar(x) && strcmpi(x,'manual')
            error('There isn`t any default value for this criterion. Please check the documentation.');
            
        elseif iscell(x) && (~strcmpi(x{1},'uniform') && ~strcmpi(x{1},'auto') && ~strcmpi(x{1},'manual'))
            error('The first element of the cell only accepts "uniform", "auto", or "manual".');
            
        elseif iscell(x) && strcmpi(x{1},'manual') && ~isnumeric(x{2})
            error('The second element should be numeric (double).');
            
        elseif iscell(x) && strcmpi(x{1},'uniform') && ~isnumeric(x{2})
            error('The second element (number of regions) should be numeric.');
            
        elseif iscell(x) && strcmpi(x{1},'uniform') && (~isequal(round(x{2}), x{2}) || numel(x{2})>1 || x{2}<2)
            error(['The second element (number of regions) should be a positive integer between 2 and ', num2str(round(size(data, 2)/2)), '.']);
            
        elseif iscell(x) && strcmpi(x{1},'auto') && (~isnumeric(x{2}) || size(x{2},1) > 1 || size(x{2},2) > 2 || ...
                size(x{2},2) < 2 || x{2}(2) <= x{2}(1) || x{2}(2) <=0 || x{2}(1) <=0)
            error('The second element structure is ([lower_bond, upper_bond]). Current implementation is not valid.');
            
        else
            TF = true;
        end
    end

addRequired(p, 'data', @check_data);
expected_on_off = {'on', 'off'};
addParameter(p, 'region_selection', default_region_selection, @check_rgs);
expected_greedy_op = {'greedy', 'genetic'};
addParameter(p, 'optimization', default_optimization, @(x) any(validatestring(x, expected_greedy_op)));
addParameter(p, 'series_fig', default_off, @(x) any(validatestring(x, expected_on_off)));
addParameter(p, 'scatter_fig', default_off, @(x) any(validatestring(x, expected_on_off)));
addParameter(p, 'output_fig', default_on, @(x) any(validatestring(x, expected_on_off)));
addParameter(p, 'wait_bar', default_on, @(x) any(validatestring(x, expected_on_off)));
addParameter(p, 'print_results', default_on, @(x) any(validatestring(x, expected_on_off)));
addParameter(p, 'PC_number', default_pc_num, @check_pcn);
addParameter(p, 'single_figure', default_off, @(x) any(validatestring(x, expected_on_off)));

parse(p, data, varargin{:});

%% DATA processing
dataset = p.Results.data;
sam_size = (1:1:size(dataset, 2))';
[M, N] = size(dataset);
[~, InX] = sort(abs(dataset(1,:))); dataset = dataset(:,InX);

if strcmpi(p.Results.series_fig, 'on')
    figure
    plot(p.Results.data');
    title('Series plot');
    xlabel('Column'); ylabel('Row');
end

if strcmpi(p.Results.scatter_fig, 'on')
    
    if isequal(size(p.Results.data,1),2)
        figure
        scatter(p.Results.data(1,:), p.Results.data(2,:));
        title('Scatter plot');
        xlabel('Var 1'); ylabel('Var 2');
    elseif isequal(size(p.Results.data,1),3)
        figure
        scatter3(p.Results.data(1,:), p.Results.data(2,:), p.Results.data(3,:));
        title('Scatter plot');
        xlabel('Var 1'); ylabel('Var 2'); zlabel('Var 3');
    else
        sprintf('The scatter plot cannot be displayed because of the number of variables.');
    end
end

rs_val = p.Results.region_selection;

if (ischar(rs_val) && strcmpi(rs_val, 'auto')) || (iscell(rs_val) && strcmpi(rs_val{1}, 'auto'))
    
    if isequal(size(rs_val), [1,2])
        optm_range = rs_val{2};
    else
        optm_range = [3, 6]; %% default value of range for "auto" option.
    end
    
    for r = optm_range(1):1:optm_range(2)
        num_gen = diff(round(linspace(0, N, r+1)));
        dataset_rv{r} = mat2cell(dataset, M, num_gen);
        
        min_num = min(num_gen);
        if ~all(num_gen == min_num)
            max_numbs = find(num_gen == max(num_gen));
            for nc = max_numbs
                dataset_rv{r}{nc} = dataset_rv{r}{nc}(:,1:min_num);
            end
        end
        
    end
    
    dt_process = 'auto';
elseif (ischar(rs_val) && strcmpi(rs_val, 'uniform')) || (iscell(rs_val) && strcmpi(rs_val{1}, 'uniform'))
    
    if isequal(size(rs_val), [1,2])
        uni_val = rs_val{2};
    else
        uni_val = 4; %% default value for "uniform" option.
    end
    
    num_gen = diff(round(linspace(0, N, uni_val+1)));
    dataset_rv = mat2cell(dataset, M, num_gen);
    
    min_num = min(num_gen);
    if ~all(num_gen == min_num)
        max_numbs = find(num_gen == max(num_gen));
        for nc = max_numbs
            dataset_rv{nc} = dataset_rv{nc}(:,1:min_num);
        end
    end
    
    dt_process = 'uniform';
    optm_range = [1, 1];
    
elseif strcmpi(rs_val{1}, 'manual')
    manual_range = unique([1, rs_val{2}, N]);
    
    if min(manual_range) < 1 || max(manual_range) > N
        error(['Manual cut(s) is(are) out of range. Valid range: 1:', num2str(N)]);
    elseif manual_range(2) < 3
        error('Number of columns in one of the regions is lower than 2.');
    else
        out_inx = discretize(sam_size, manual_range);
        
        for m = 1:1:out_inx(end)
            dataset_rv{m} = dataset(:,out_inx == m);
        end
        
        optm_range = repmat(numel(manual_range)-1, [1, 2]);
        dt_process = 'manual';
    end
    
end

%% Standardization and Model
if strcmpi(p.Results.wait_bar, 'on')
    wb = waitbar(0, 'Please wait...');
end

if ~strcmpi(dt_process, 'auto')
    rg_size = 1:1:numel(dataset_rv);
    rg_mat = [rg_size.' flipud(nchoosek(rg_size, numel(rg_size)-1))];
    
    loop_en = size(rg_mat,1);
    
    dataset_rv_base = dataset_rv;
    for lr = 1:1:loop_en
        
        % Standardization
        [~, std_info] = mapstd(dataset_rv_base{rg_mat(lr,1)});
        
        for v = 1:1:numel(dataset_rv)
            dataset_rv{v} = mapstd('apply', dataset_rv_base{v}, std_info);
        end
        
        [summ, bounds, rig_num] = pca_nl_core_test(dataset_rv(:,rg_mat(lr,:)), M, N , p.Results.optimization, p.Results.PC_number);
        results.output(lr,:) = summ; results.bounds(lr,:) = bounds; results.rig_num(lr,:) = rig_num;
        
        if strcmpi(p.Results.wait_bar, 'on')
            waitbar(lr / loop_en)
        end
    end
    
else
    
    for xr = optm_range(1):1:optm_range(2)
        rg_size = 1:1:numel(dataset_rv{xr});
        rg_mat{xr} = [rg_size.' flipud(nchoosek(rg_size, numel(rg_size)-1))];
        
        dataset_rv_base = dataset_rv;
        for lr = 1:1:size(rg_mat{xr},1)
            % Standardization
            [~, std_info] = mapstd(dataset_rv_base{xr}{rg_mat{xr}(lr,1)});
            
            for v = 1:1:numel(dataset_rv{xr})
                dataset_rv{xr}{v} = mapstd('apply', dataset_rv_base{xr}{v}, std_info);
            end
            
            [summ, bounds, rig_num] = pca_nl_core_test(dataset_rv{xr}(:,rg_mat{xr}(lr,:)), M, N, p.Results.optimization, p.Results.PC_number);
            results{xr}.output(lr,:) = summ;
            results{xr}.bounds(lr,:) = bounds;
            results{xr}.rig_num(lr,:) = rig_num;
        end
        
        if strcmpi(p.Results.wait_bar, 'on')
            waitbar(xr / optm_range(2))
        end
    end
end

%% Results
if strcmpi(p.Results.wait_bar, 'on')
    close(wb);
end

if ~strcmpi(dt_process, 'auto')
    if strcmpi(p.Results.output_fig, 'on')
        for xr = optm_range(1):1:optm_range(2)
            num_of_plots = size(rg_mat,1);
            for lr = 1:1:num_of_plots
                x = 0:.1:rig_num+1;
                
                if strcmpi(p.Results.single_figure, 'on') && num_of_plots>=2 && num_of_plots<=8
                    
                    if lr == 1
                        figure
                    end
                    
                    if num_of_plots == 2
                        subplot(1,2,lr);
                        set(gcf,'units','normalized','outerposition',[0 0 1 1])
                        single_figure_check = 0;
                    elseif num_of_plots == 3 || num_of_plots == 4
                        subplot(2,2,lr);
                        set(gcf,'units','normalized','outerposition',[0 0 1 1])
                        single_figure_check = 0;
                    elseif num_of_plots == 5 || num_of_plots == 6
                        subplot(3,2,lr);
                        set(gcf,'units','normalized','outerposition',[0 0 1 1])
                        single_figure_check = 0;
                    elseif num_of_plots == 7 || num_of_plots == 8
                        subplot(4,2,lr);
                        set(gcf,'units','normalized','outerposition',[0 0 1 1])
                        single_figure_check = 0;
                    end
                else
                    figure
                    single_figure_check = 1;
                end
                
                plot(x, results.bounds(lr,1)*ones(1,length(x)))
                hold on
                plot(x, results.bounds(lr,2)*ones(1,length(x)))
                
                if single_figure_check == 1 || lr == 1 || lr == 2
                    ylabel('Residual variances');
                    title(['Accuracy bounds and results for region #', num2str(rg_mat(lr,1))]);
                else
                    title(['for region #', num2str(rg_mat(lr,1))]);
                end
                
                for i = 1:rig_num
                    hold on
                    plot(i, results.output(lr,i), '+', 'LineWidth', 1.3, 'MarkerSize', 7, 'MarkerEdgeColor','b')
                end
                hold on
                plot(1, results.output(lr,1), 'o', 'LineWidth', 1.3, 'MarkerSize', 7, 'MarkerEdgeColor', 'r');
                
            end
        end
    end
    
else
    
    if strcmpi(p.Results.output_fig, 'on')
        for xr = optm_range(1):1:optm_range(2)
            loop_cx = 1;
            num_of_plots = size(rg_mat{xr},1);
            for lr = 1:1:num_of_plots
                x = 0:.1:xr+1;
                
                if strcmpi(p.Results.single_figure, 'on') && num_of_plots>=2 && num_of_plots<=8
                    
                    if loop_cx == 1
                        figure
                    end
                    
                    if num_of_plots == 2
                        subplot(1,2,lr);
                        set(gcf,'units','normalized','outerposition',[0 0 1 1])
                        single_figure_check = 0;
                    elseif num_of_plots == 3 || num_of_plots == 4
                        subplot(2,2,lr);
                        set(gcf,'units','normalized','outerposition',[0 0 1 1])
                        single_figure_check = 0;
                    elseif num_of_plots == 5 || num_of_plots == 6
                        subplot(3,2,lr);
                        set(gcf,'units','normalized','outerposition',[0 0 1 1])
                        single_figure_check = 0;
                    elseif num_of_plots == 7 || num_of_plots == 8
                        subplot(4,2,lr);
                        set(gcf,'units','normalized','outerposition',[0 0 1 1])
                        single_figure_check = 0;
                    end
                else
                    figure
                    single_figure_check = 1;
                end
                
                plot(x, results{xr}.bounds(lr,1)*ones(1,length(x)))
                hold on
                plot(x, results{xr}.bounds(lr,2)*ones(1,length(x)))
                
                if single_figure_check == 1 || lr == 1 || lr == 2
                    ylabel('Residual variances');
                    title(['Accuracy bounds and results for region #', num2str(rg_mat{xr}(lr,1)), ' - ','num of regions: ', num2str(xr)]);
                else
                    title(['for region #', num2str(rg_mat{xr}(lr,1)), ' - ','num of regions: ', num2str(xr)]);
                end
                
                for i = 1:xr
                    hold on
                    plot(i, results{xr}.output(lr,i), '+', 'LineWidth', 1.3, 'MarkerSize', 7, 'MarkerEdgeColor','b')
                end
                
                hold on
                plot(1, results{xr}.output(lr,1), 'o', 'LineWidth', 1.3, 'MarkerSize', 7, 'MarkerEdgeColor', 'r');
                
                loop_cx = 0;
            end
        end
    end
end

for or = optm_range(1):1:optm_range(2)
    
    if strcmpi(dt_process, 'auto')
        T_bounds = results{or}.bounds; LowerBound = T_bounds(:,1); UpperBound = T_bounds(:,2);
        num_of_regions = numel(results{or}.rig_num);
        Regions_output = results{or}.output;
    else
        T_bounds = results.bounds; LowerBound = T_bounds(:,1); UpperBound = T_bounds(:,2);
        num_of_regions = numel(results.rig_num);
        Regions_output = results.output;
    end
    
    rg_str = cell(num_of_regions, 1);
    
    for rp = 1:1:num_of_regions
        add_rg_str = ['Benchmark - region ', num2str(rp)];
        rg_str{rp} = add_rg_str;
    end
    
    InRegion = Regions_output>kron(LowerBound, ones(1, num_of_regions)) & Regions_output<kron(UpperBound, ones(1, num_of_regions));
    
    InRegion_per = (sum(InRegion, 2)/size(InRegion, 2))*100;
    
    if strcmpi(dt_process, 'auto')
        tbx_1{or} = table(Regions_output, LowerBound, UpperBound, 'RowNames', rg_str);
        tbx_2{or} = table(InRegion, InRegion_per, 'RowNames', rg_str);
        
    else
        tbx_1 = table(Regions_output, LowerBound, UpperBound, 'RowNames', rg_str);
        tbx_2= table(InRegion, InRegion_per, 'RowNames', rg_str);
    end
end

if  strcmpi(p.Results.print_results, 'on') && ~strcmpi(dt_process, 'auto')
    disp('------------------------');
    disp('Numerical report: '); disp(tbx_1); disp('------------------------');
    disp('In/Out of region report: '); disp(tbx_2);
end

options.raw_data = p.Results.data;
options.region_selection = p.Results.region_selection;
options.optimization = p.Results.optimization;
options.output_fig = p.Results.output_fig;
options.series_fig = p.Results.series_fig;
options.scatter_fig = p.Results.scatter_fig;
options.wait_bar = p.Results.wait_bar;
options.print_results = p.Results.print_results;
options.number_of_principal_components = p.Results.PC_number;
options.all_data_size = size(data);
options.results.output_table = tbx_1;
options.results.in_out_table = tbx_2;

end