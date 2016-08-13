function options = makeParameterStructure(options)

    if nargin < 1
        options = [];
    end

    
    %# of modes in clustering PCA (default = 50)
    template_pca_dimension  = 50;
    
    %first mode in the PCA to use in likelihood model (default = 1)
    first_mode = 1;
    
    %number of clusters in kmeans (default = 12)
    k = 12;
    
    %maximum # of peaks to use in clustering & GMM (default = 10000);
    maxNumPeaks = 10000;
        
    %# of replicates in kmeans (default = 5)
    kmeans_replicates = 5;
    
    %max # of iterations for kmeans (default = 10000)
    kmeans_maxIter = 10000;    
    
    %sampling frequency of original data (default = 1e4)
    fs = 1e4;
    
    %minimum # of peaks in a template (default = 50)
    min_template_size = 50;
    
    %maximum number of gaussians to use in noise fit (default = 2)
    maxNumGaussians_noise = 2;
    
    %number of replicates to use in GMM calculations (default = 3)
    replicates_GMM = 3;
    
    %smoothing length for noise filter in milliseconds (default = 4)
    smoothingLength_noise = 4;
    
    %minimum region length for signal determination in milliseconds
    %(default = 4)
    minRegionLength = 4;
    
    %maximum IPI to fit in milliseconds (default = 500)
    maxIPI = 500;
    
    %smoothing sigma for IPI kernel density estimation in milliseconds
    %(default = .05)
    IPI_sigma = .05;
    
    %number of FWHM widths if IPI distribution to set for diffThreshold
    %(default = 2)
    num_IPI_halfWidths = 2;
    
    %percentage of peak amplitudes below the nosie threshold to be called a
    %"noise" template (default = .9)
    amplitude_threshold = .9;
    
    
    if ~isfield(options,'template_pca_dimension') || isempty(options.template_pca_dimension)
        options.template_pca_dimension = template_pca_dimension;
    end
    
    if ~isfield(options,'first_mode') || isempty(options.first_mode)
        options.first_mode = first_mode;
    end
    
    if ~isfield(options,'k') || isempty(options.k)
        options.k = k;
    end
    
    if ~isfield(options,'maxNumPeaks') || isempty(options.maxNumPeaks)
        options.maxNumPeaks = maxNumPeaks;
    end
    
    if ~isfield(options,'kmeans_replicates') || isempty(options.kmeans_replicates)
        options.kmeans_replicates = kmeans_replicates;
    end
    
    if ~isfield(options,'kmeans_maxIter') || isempty(options.kmeans_maxIter)
        options.kmeans_maxIter = kmeans_maxIter;
    end
    
    if ~isfield(options,'fs') || isempty(options.fs)
        options.fs = fs;
    end
    
    if ~isfield(options,'min_template_size') || isempty(options.min_template_size)
        options.min_template_size = min_template_size;
    end
    
    if ~isfield(options,'maxNumGaussians_noise') || isempty(options.maxNumGaussians_noise)
        options.maxNumGaussians_noise = maxNumGaussians_noise;
    end
    
    if ~isfield(options,'replicates_GMM') || isempty(options.replicates_GMM)
        options.replicates_GMM = replicates_GMM;
    end
    
    if ~isfield(options,'smoothingLength_noise') || isempty(options.smoothingLength_noise)
        options.smoothingLength_noise = smoothingLength_noise;
    end
    
    if ~isfield(options,'minRegionLength') || isempty(options.minRegionLength)
        options.minRegionLength = minRegionLength;
    end
    
    if ~isfield(options,'maxIPI') || isempty(options.maxIPI)
        options.maxIPI = maxIPI;
    end
    
    if ~isfield(options,'IPI_sigma') || isempty(options.IPI_sigma)
        options.IPI_sigma = IPI_sigma;
    end
    
    if ~isfield(options,'num_IPI_halfWidths') || isempty(options.num_IPI_halfWidths)
        options.num_IPI_halfWidths = num_IPI_halfWidths;
    end
    
    if ~isfield(options,'amplitude_threshold') || isempty(options.amplitude_threshold)
        options.amplitude_threshold = amplitude_threshold;
    end
    
    
    