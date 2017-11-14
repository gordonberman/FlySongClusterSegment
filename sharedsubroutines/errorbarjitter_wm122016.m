function errorbarjitter(data,h,varargin)

% SYNTAX
%
%   errorbarjitter(data)
%   errorbarjitter(data,param1,val1,param2,val2,...)
%
%   e.g.
%   errorbarjitter(data,'barends','yes','colors',color_array)
%
%
% Plots mean (or median)±SD plus jitter plot of raw data
% Columns are categories, rows are individual samples.
%
% PARAMETERS (all optional)
%
% OFFSET: plots with user defined offset of mean±SD and jittered raw data
%
% AVERAGE: plot with either mean ('mean') or median ('median') or none
%   ('none')
%
% SORT: plot with either sorted data ('sort') or user entered 
%   order ('nosort'). Default is nosort.
%
% COLHEADERS: plot with user defined column headers. 
%   It is a bad idea to sort without providing colheaders,
%   since it may not be easy to track the source of the data. 
%
% FACTOR: plot with user defined factor for scaling the jitter
%
% LEFT_OR_RIGHT: when LEFT_OR_RIGHT = 'left', plot with mean±SD line 
%   on left, when = 'right', plot with bar on right (default)
%
% BARENDS: when barends = 'yes', plot with capped ends of SD
%   bars, when = 'no' (default), plot without barends
%
% STD: Provide a separate array containing the standard deviations
%   (or any other estimates of distribution) for each datum in d, and an error
%   bar equal to ± these deviations will be plotted on each datum
%
% COLORS: Plot data with user defined colors, using standard Matlab nomenclature, 
%   specified as cell array equal in length to the number of data columns.
%
% COLOR_OPTS: A three element logical vector that indicates whether to 
%   color the raw data (position 1), the average point (position 2), and/or 
%   the bar representing the standard deviation estimate (position 3). e.g. [0 1 1] 
%   would color the average point and the variance line
%
% COLOR_MAP: Provide a matrix equal in dimension to the data matrix to plot
%   each datum a color from the scale of the current colormap
%
% DATA_MARKER: Use Matlab Marker specifier for data points ('o' default)
%
% DATA_MARKER_SIZE: Size of markers used for data points
%
% AVE_MARKER: Use Matlab Marker specifier for average points ('o' default)
%
% AVE_MARKER_SIZE: Size of markers used for average points
%
% DELETE_OUTLIERS: Delete outliers using Grubb's test, as implemented by
%   Brett Shoelson (www.mathworks.com/matlabcentral/fileexchange/3961-deleteoutliers)
%   Provide 2 element array containing alpha and rep: e.g. [.05 1]
%
% PLOT_SDS: Plot horizontal lines corresponding to the # of standard
%   deviations of the data in the first column. e.g. 2
%
% SAVE_DIRECT: Save figure directly to specified format without plotting in
%   Matlab window. To use, specify output filename, which must include
%   legitimate suffix for Matlab saveas command. e.g. 'plot.png'. NOTE:
%   If you save in Matlab .fig format (i.e. 'plot.fig'), the figure is saved in 
%   'invisible' mode. You can open the plot in visible mode with
%   openfig('plot.fig','new','visible') or, after opening, set(gcf,'Visible','on')
%
% ACKNOWLEDGEMENT: This function depends on jitter.m, written by Richie
%   Cotton 2006/03/21.
% 
%
%
% $	Author: David Stern	$   $   Date :2011/11/02   $
%
% Bug Fixes and Improvements
% 2012/09/27    Added ability to plot one column of data.
% 2012/10/16    Fixed bug that caused crash when tried to sort multiple 
%               samples with the same mean
%               Add ability to plot error bars on each sample
%               Added ability to plot different groups with different
%               colors
% 2012/10/26    Bug fixes to std and colors options. Plot failed when
%               options missing.
% 2012/10/27    Added option to save figure directly without plotting
% 2012/10/30    Changed input format to handle many variables
%               Control colors of point, averages, and lines separately 
%               Plot using regular Matlab symbols
% 2012/12/10    Added ability to define average and data marker sizes
% 2012/12/12    Plot with tight X axis
%               Added ability to pass axis handle
% 2013/05/10    Added ability to delete outliers and to plot horizontal SDs
% 2013/05/14    Bug fix fpr plotting horizontal SDs
% 2013/12/20    Added option to specify matrix of values for datapoint
%               colors & option to plot without average±SD
% 2014/1/8      Small bug fix of marker size for plots with colormap



% Check number of inputs
if nargin < 1
    error('plot_u_sd_jitterraw:notEnoughInputs', 'This function requires at least one input.');
end

if nargin < 2
    
    h = gcf;
end

d = data;
n_categories = size(d,2);
n_data = size(d,1);

%establish input argument parser

p = inputParser;

%Set default options

defaultOffset = 0.2;
defaultAverage = 'mean';
defaultSort = 'nosort';
defaultColheaders = cell(1,n_categories);
defaultFactor = 1;
defaultLeft_or_right = 'right';
defaultBarends = 'no';
defaultXsep = 1;
defaultColor_Opts = [1 0 0];
defaultColor_Map = 'no';
defaultData_marker = 'o';
defaultData_marker_size = 50;
defaultAve_marker = 'o';
defaultAve_marker_size = 50;
defaultDelete_Outliers = 'no';
defaultPlot_SDs = 'no';

std = zeros(n_data,n_categories);
std(:) = NaN;
defaultStd =std;

colors = cell(1,n_categories);
colors(:) = {'black'};
defaultColors = colors;
cmap=colormap;

defaultSave_direct = [];

%add required and optional inputs to parser

addRequired(p,'data',@isnumeric);
addOptional(p,'offset',defaultOffset,@isnumeric);
addOptional(p,'average',defaultAverage);
addOptional(p,'sort',defaultSort);
addOptional(p,'colheaders',defaultColheaders);
addOptional(p,'factor',defaultFactor,@isnumeric);
addOptional(p,'left_or_right',defaultLeft_or_right);
addOptional(p,'barends',defaultBarends);
addOptional(p,'Xsep',defaultXsep,@isnumeric);%not currently available as user-specific variable - does not plot pretty
addOptional(p,'std',defaultStd);
addOptional(p,'colors',defaultColors);
addOptional(p,'save_direct',defaultSave_direct);
addOptional(p,'color_opts',defaultColor_Opts);
addOptional(p,'color_map',defaultColor_Map);
addOptional(p,'data_marker',defaultData_marker);
addOptional(p,'data_marker_size',defaultData_marker_size);
addOptional(p,'ave_marker',defaultAve_marker);
addOptional(p,'ave_marker_size',defaultAve_marker_size);
addOptional(p,'delete_outliers',defaultDelete_Outliers);
addOptional(p,'plot_SDs',defaultPlot_SDs);

parse(p,data,varargin{:});

%redistribute parsed variables to original variable names

offset = p.Results.offset;
average = p.Results.average;
sort_or_nosort = p.Results.sort;
colheaders = p.Results.colheaders;
factor = p.Results.factor;
left_or_right = p.Results.left_or_right;
barends = p.Results.barends;
Xsep = p.Results.Xsep;
std = p.Results.std;
colors = p.Results.colors;
save_direct = p.Results.save_direct;
color_opts = p.Results.color_opts;
color_map = p.Results.color_map;
data_marker = p.Results.data_marker;
data_marker_size = p.Results.data_marker_size;
ave_marker = p.Results.ave_marker;
ave_marker_size = p.Results.ave_marker_size;
delete_outliers = p.Results.delete_outliers;
plot_SDs = p.Results.plot_SDs;

%set some other variables

if strcmp(left_or_right,'left') == 1
    offset = -offset;
end

if sum(cellfun(@isempty, colheaders)) == n_categories
    skip_colheaders = 1;
    if strcmp(sort_or_nosort,'sort') == 1
        fprintf('It is a bad idea to sort without providing colheaders.\n')
        fprintf('Good luck keeping track of your data!\n')
    end
else
    skip_colheaders = 0;
end

if ~isempty(save_direct)
    set(h,'Visible','off');
else
    set(h,'Visible','on');
end

if ~strcmp(delete_outliers,'no')
    for i = 1:n_categories
        d(:,i) = deleteoutliers(d(:,i),delete_outliers(1),delete_outliers(2));
    end
end

if ~strcmp(average,'none')
    if strcmp(average,'mean') == 1
        mean_d = nanmean(d);
    elseif strcmp(average,'median') == 1
        mean_d = nanmedian(d);
    end
    mean_original_d = mean_d;
    std_original_d = nanstd(d);
end

if strcmp(average,'none')
    offset = 0;
end

%make figure

%figure(h)
hold on

%add option to sort data by mean
%there must be an easier way to rearrange an array by a property of columns
if strcmp(sort_or_nosort,'sort') == 1
    
    [sorted_means,sort_idx] = sort(mean_d);    
    
    %now put original array in new order
    %if data in colheaders, sort that one too
    %there must be an easy way to do this
    sorted_d = d(:,sort_idx);
    d = sorted_d;
    if ~isempty(colors)
        sorted_colors = colors(sort_idx);
        colors = sorted_colors;
    end
    if ~isempty(std)
        sorted_std = std(:,sort_idx);
        std = sorted_std;
    end

    
    if skip_colheaders ~= 1
        sorted_colheaders = colheaders(sort_idx);
        colheaders = sorted_colheaders;
    end
    %recalculate mean for newly sorted data
    if strcmp(average,'mean') == 1
        mean_d = nanmean(d);
    elseif strcmp(average,'median') == 1
        mean_d = nanmedian(d);
    end
    
    std_d = nanstd(d);
end



%plot SD lines for control data (column 1)

if ~strcmp(plot_SDs,'no')
    SD = std_original_d(:,1);
    M = mean_original_d(:,1);
    line([0 n_categories+1],[M+plot_SDs*SD M+plot_SDs*SD],'LineWidth',0.5,'Color',[0.5 0.5 0.5])
    line([0 n_categories+1],[M-plot_SDs*SD M-plot_SDs*SD],'LineWidth',0.5,'Color',[0.5 0.5 0.5])
end

%for column in data
%plot (mean ±SD)
e = nanstd(d,1);
%define X axis positions
x = 1:1:n_categories;
x = x*Xsep;

%put column indices in each column
    
%plot mean and error bar
if ~strcmp(average,'none')
    if strcmp(barends,'no') == 1
        if color_opts(2) == 0
            if ave_marker == 'o'
                scatter(x+offset,mean_d,ave_marker_size,ave_marker,'k','filled');
            else
                scatter(x+offset,mean_d,ave_marker_size,ave_marker,'k');
            end
        else
            for i = 1:n_categories
                if ave_marker == 'o'
                    scatter(x(i)+offset,mean_d(i),ave_marker_size,ave_marker,'filled',colors{i});
                else
                    scatter(x(i)+offset,mean_d(i),ave_marker_size,ave_marker,colors{i});
                end
            end
        end
    else %draw with barends
        if color_opts(2) == 0
            if color_opts(3) == 0
                errorbar(mean_d,e,ave_marker,'MarkerFaceColor','k','MarkerEdgeColor','k','XData',x+offset);
            else
                for i = 1:n_categories
                    errorbar(mean_d(i),e(i),ave_marker,'MarkerFaceColor','k','MarkerEdgeColor','k','Color',colors{i},'XData',x(i)+offset);
                end
            end
        else
            if color_opts(3) == 1
                for i = 1:n_categories
                    errorbar(mean_d(i),e(i),ave_marker,'MarkerFaceColor',colors{i},'Color',colors{i},'XData',x(i)+offset);
                end
            else
                for i = 1:n_categories
                    errorbar(mean_d(i),e(i),ave_marker,'MarkerFaceColor',colors{i},'MarkerEdgeColor',colors{i},'Color','k','XData',x(i)+offset);
                end
            end
        end
    end
    
    %plot error lines
    for i = 1:n_categories
        if color_opts(3) == 0
            line([x(i)+offset x(i)+offset],[mean_d(i)-e(i) mean_d(i)+e(i)],'Color','k','LineWidth',.5);
        else
            line([x(i)+offset x(i)+offset],[mean_d(i)-e(i) mean_d(i)+e(i)],'Color',colors{i},'LineWidth',.5);
        end
    end
end

%plot raw data with jitter in x axis to left of each
if n_categories >1
    x = repmat(x,n_data,1);
    x = jitter(x,factor);
else %if have one column, need to add false second column to allow jitter to work properly, and then delete
    x=repmat(x,n_data,1);
    x(:,2) = 2;
    x = jitter(x,factor);
    x(:,2) = [];
end
x(isnan(d)) = NaN;

for i = 1:n_categories
    %if ismatrix(color_map)
    if size(x,1) == size(color_map,1) %WKM 12/20/16: Modified because 'no' 
                                      %(default for color_map) is technically a matrix
        try
            scatter(x(:,i)-offset,d(:,i),data_marker_size,color_map(:,i),'marker',data_marker)
        catch MEpropval
            length(x(:,i)-offset)
            length(d(:,i))
            length(color_map(:,i))
            rethrow(MEpropval)
        end
    else
        if color_opts(1) == 1
            scatter(x(:,i)-offset,d(:,i),data_marker_size,data_marker,'MarkerEdgeColor',colors{i})
        else
            scatter(x(:,i)-offset,d(:,i),data_marker_size,data_marker,'MarkerEdgeColor','k')
        end
    end
    
end

%plot error lines for each datum
if ~isempty(std)
    for i = 1:n_categories
        for j = 1:n_data
            %if ismatrix(color_map)
            if size(x,1) == size(color_map,1) %WKM 12/20/16: Modified because 'no' 
                                      %(default for color_map) is technically a matrix
                if ~isnan(color_map(j,i))
                    coloridx = ceil(max(color_map(j,i) * 64)); %get colormap value of datum
                    line([x(j,i)-offset x(j,i)-offset],[d(j,i)-std(j,i) d(j,i)+std(j,i)],'Color',cmap(coloridx,:),'LineWidth',.5)
                end
            else
                
                if color_opts(1) == 1
                    line([x(j,i)-offset x(j,i)-offset],[d(j,i)-std(j,i) d(j,i)+std(j,i)],'Color',colors{i},'LineWidth',.5)
                else
                    line([x(j,i)-offset x(j,i)-offset],[d(j,i)-std(j,i) d(j,i)+std(j,i)],'Color','k','LineWidth',.5)
                end
            end
        end
    end
end



set(gca,'XTick',[1:1:n_categories],'TickDir','Out','Xlim',[0 n_categories + 1])

if skip_colheaders == 0
    %add x axis labels
    set(gca,'XTickLabel',colheaders)
else
    set(gca,'XTickLabel',[])    
end
hold off
if ~isempty(save_direct);
    saveas(gcf,save_direct);
end

function y = jitter(x, factor, uniformOrGaussianFlag, smallOrRangeFlag, realOrImaginaryFlag)
% Adds a small amount of noise to an input vector, matrix or N-D array. The
% noise can be uniformly or normally distributed, and can have a magnitude
% based upon the range of values of X, or based upon the smallest
% difference between values of X (excluding 'fuzz').
% 
% NOTE: This function accepts complex values for the first input, X.  If
% any values of X have imaginary components (even zero-valued imaginary
% components), then by default the noise will be imaginary.  Otherwise, the
% default is for real noise.  You can choose between real and imaginary
% noise by setting the fifth input parameter (see below).
% 
% Y = JITTER(X) adds an amount of uniform noise to the input X, with a
% magnitude of one fifth of the smallest difference between X values
% (excluding 'fuzz'), i.e. the noise, n~U(-d/5, d/5), where d is the
% smallest difference between X values.
% 
% Y = JITTER(X, FACTOR) adds noise as above, but scaled by a factor
% of FACTOR, i.e. n~U(-FACTOR*d/5, FACTOR*d/5).
% 
% Y = JITTER(X, FACTOR, 1) adds noise as above, but normally distributed
% (white noise), i.e. n~N(0, FACTOR*d/5). JITTER(X, FACTOR, 0) works the
% same as JITTER(X, FACTOR). If the second parameter is left empty (for
% example JITTER(X, [], 1)), then a default scale factor of 1 is used.
% 
% Y = JITTER(X, FACTOR, [], 1) adds an amount of noise to X with a
% magnitude of one fiftieth of the range of X.  JITTER(X, FACTOR, [], 0)
% works the same as JITTER(X, FACTOR, []).  A value of 0 or 1 can be given as
% the third input to choose between uniform and normal noise (see above),
% i.e. n~U(-FACTOR*r/50, FACTOR*r/50) OR n~N(0, FACTOR*r/50), where r is
% the range of the values of X.  If the second parameter is left empty then
% a default scale factor of 1 is used.
% 
% Y = JITTER(X, FACTOR, [], [], 1) adds an amount of noise as above, but
% with imaginary noise.  The magnitude of the noise is the same as in the
% real case, but the phase angle is a uniform random variable, theta~U(0,
% 2*pi).  JITTER(X, FACTOR, [], [], 0) works the same as JITTER(X, FACTOR,
% [], []). A value of 0 or 1 can be given as the third input to choose
% between uniform and normal noise, and a value of 0 or 1 can be given as
% the fourth input to choose between using the smallest distance between
% values or the range for determining the magnitude of the noise.  If the
% second parameter is left empty then a default scale factor of 1 is used.
% 
% 
% EXAMPLE:  x = [1 -2 7; Inf 3.5 NaN; -Inf 0.001 3];
%           jitter(x)
% 
%           ans =
% 
%             0.9273   -2.0602    6.9569
%                Inf    3.4597       NaN
%               -Inf    0.0333    2.9130
% 
%           %Plot a noisy sine curve. 
%           x2 = sin(0:0.1:6);
%           plot(jitter(x2, [], 1, 1));  
% 
% 
% ACKNOWLEGEMENT: This function is based upon the R function of the same
% name, written by Werner Stahel and Martin Maechler, ETH Zurich.
% See http://stat.ethz.ch/R-manual/R-patched/library/base/html/jitter.html
% for details of the original.
% 
% 
%   Class support for input X:
%      float: double, single
% 
% 
%   See also RAND, RANDN.
% 
% 
% $ Author: Richie Cotton $     $ Date: 2006/03/21 $


% Check number of inputs
if nargin < 1
    error('jitter:notEnoughInputs', 'This function requires at least one input.');
end
    
% Set defaults where required
if nargin < 2 || isempty(factor)
    factor = 1;
end

if nargin < 3 || isempty(uniformOrGaussianFlag)
    uniformOrGaussianFlag = 0;
end

if nargin < 4 || isempty(smallOrRangeFlag)
    smallOrRangeFlag = 0;
end

if nargin < 5 || isempty(realOrImaginaryFlag)
    realOrImaginaryFlag = ~isreal(x);
end


% Find the range of X, ignoring infinite value and NaNs
xFinite = x(isfinite(x(:)));
xRange = max(xFinite) - min(xFinite);

if ~smallOrRangeFlag
    % Remove 'fuzz'
    dp = 3 - floor(log10(xRange));
    xFuzzRemoved = round(x * 10^dp) * 10^-dp;
    % Find smallest distance between values of X
    xUnique = unique(sort(xFuzzRemoved));
    xDifferences = diff(xUnique);
    if length(xDifferences)
        smallestDistance = min(xDifferences);
    elseif xUnique ~= 0 
        % In this case, all values are the same, so xUnique has length 1
        smallestDistance = 0.1 * xUnique;
    else
        % In this case, all values are 0
        smallestDistance = 0.1 * xRange;
    end
    scaleFactor = 0.2 * factor * smallestDistance;
else
    % Calc scale factor based upon range
    scaleFactor = 0.02 * factor * xRange;
end

% Add the noise
s = size(x);
if uniformOrGaussianFlag
    % Normal noise
    if realOrImaginaryFlag
        randomPhaseAngles = 2 * pi * rand(s);
        y = x + scaleFactor * randn(s) * exp(randomPhaseAngles * i);
    else
        y = x + scaleFactor * randn(s);
    end
else
    % Uniform noise
    if realOrImaginaryFlag
        randomPhaseAngles = 2 * pi * rand(s);
        y = x + scaleFactor * (2*rand(s)-1) * exp(randomPhaseAngles * i);
    else
        y = x + scaleFactor * (2*rand(s)-1);
    end
end


