function options = makeParameterStructure(options)

    if nargin < 1
        options = [];
    end

    
    %run t-SNE embedding in createTemplates.m (default = false)
    run_tsne = false;
    
    %# of modes in clustering PCA (default = 50)
    template_pca_dimension  = 50;
    
    %first mode in the PCA to use in likelihood model (default = 1)
    first_mode = 1;
    
    %number of clusters in kmeans (default = 12)
    k = 12;
    
    %true if running human labeling of noise templates (default = false)
    humanLabel = false;
    
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
    
    %maximum carrier frequency in Hz (default = 1000)
    maxCarrierFrequency = 1000;
    
    %smoothing sigma for IPI kernel density estimation in milliseconds
    %(default = 1 ms)
    IPI_sigma = 1;
    
    %number of FWHM widths of IPI distribution to set for diffThreshold
    %(default = 2)
    num_IPI_halfWidths = 2;
    
    %percentage of peak amplitudes below the noise threshold to be called a
    %"noise" template (default = .75)
    amplitude_threshold = .75;
    
    %minimum noise threshold (set to be negative to ignore, default = -1)
    min_noise_threshold = -1;
    
    %length of data median filter in milliseconds (default = 1, not used if < 0)
    median_filter_length = 1;
    
    %minimum seperation between peaks in milliseconds (default = 1)
    min_seperation = 1;
    
    %noise posterior threshold (between 0 and 1, 1 most stringent, default = .5)
    noise_posterior_threshold = .5;
    
    %multiple of diff threshold, any small peaks within
    %diff_threshold_multiplier*diff_threshold are ignored (default = .75)
    diff_threshold_multiplier = .75;
    
    %high pass filter cut-off on data set in Hz (default = -1, < 0 to not filter)
    high_pass_filter_cutoff = -1;
    
    %butterworth high-pass filter cut-off order (default = 6)
    butterworth_order = 6;
    
    %Toggle for using likelihood to refine clusterings (default = false)
    refine_clusters = false;
    
    
    %%%%%%%% t-SNE options %%%%%%%%
    
    
    %2^H (H is the transition entropy)
    perplexity = 32;
    
    %relative convergence criterium for t-SNE
    relTol = 1e-4;
    
    %number of dimensions for use in t-SNE
    num_tsne_dim = 2;
    
    %binary search tolerance for finding pointwise transition region
    sigmaTolerance = 1e-5;
    
    %maximum number of non-zero neighbors in P
    maxNeighbors = 200;
    
    %initial momentum
    momentum = .5;
    
    %value to which momentum is changed
    final_momentum = 0.8;    
    
    %iteration at which momentum is changed
    mom_switch_iter = 250;      
    
    %iteration at which lying about P-values is stopped
    stop_lying_iter = 125;      
    
    %degree of P-value expansion at early iterations
    lie_multiplier = 4;
    
    %maximum number of iterations
    max_iter = 1000;  
    
    %initial learning rate
    epsilon = 500;  
    
    %minimum gain for delta-bar-delta
    min_gain = .01;   

    %readout variable for t-SNE
    tsne_readout = 10;
    
    %embedding batchsize
    embedding_batchSize = 20000;
    
    %maximum number of iterations for the Nelder-Mead algorithm
    maxOptimIter = 100;
    
    %number of points in the training set
    trainingSetSize = 35000;
    
    %local neighborhood definition in training set creation
    kdNeighbors = 5;
    
    %t-SNE training set stopping critereon
    training_relTol = 2e-3;
    
    %t-SNE training set perplexity
    training_perplexity = 20;
    
    %number of points to evaluate in each training set file
    training_numPoints = 10000;
    
    %minimum training set template length
    minTemplateLength = 1;
    
    
    
    if ~isfield(options,'high_pass_filter_cutoff') || isempty(options.high_pass_filter_cutoff)
        options.high_pass_filter_cutoff = high_pass_filter_cutoff;
    end
    
    if ~isfield(options,'butterworth_order') || isempty(options.butterworth_order)
        options.butterworth_order = butterworth_order;
    end
    
    if ~isfield(options,'diff_threshold_multiplier') || isempty(options.diff_threshold_multiplier)
        options.diff_threshold_multiplier = diff_threshold_multiplier;
    end
    
    if ~isfield(options,'noise_posterior_threshold') || isempty(options.noise_posterior_threshold)
        options.noise_posterior_threshold = noise_posterior_threshold;
    end
    
    if ~isfield(options,'min_seperation') || isempty(options.min_seperation)
        options.min_seperation = min_seperation;
    end
    
    if ~isfield(options,'median_filter_length') || isempty(options.median_filter_length)
        options.median_filter_length = median_filter_length;
    end
    
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
    
    if ~isfield(options,'humanLabel') || isempty(options.humanLabel)
        options.humanLabel = humanLabel;
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
    
    if ~isfield(options,'min_noise_threshold') || isempty(options.min_noise_threshold)
        options.min_noise_threshold = min_noise_threshold;
    end
    
    if ~isfield(options,'perplexity') || isempty(options.perplexity)
        options.perplexity = perplexity;
    end
    
    
    if ~isfield(options,'relTol') || isempty(options.relTol)
        options.relTol = relTol;
    end
    
    
    if ~isfield(options,'num_tsne_dim') || isempty(options.num_tsne_dim)
        options.num_tsne_dim = num_tsne_dim;
    end
    
    
    if ~isfield(options,'sigmaTolerance') || isempty(options.sigmaTolerance)
        options.sigmaTolerance = sigmaTolerance;
    end
    
    
    if ~isfield(options,'maxNeighbors') || isempty(options.maxNeighbors)
        options.maxNeighbors = maxNeighbors;
    end
    
    
    if ~isfield(options,'momentum') || isempty(options.momentum)
        options.momentum = momentum;
    end
    
    
    if ~isfield(options,'final_momentum') || isempty(options.final_momentum)
        options.final_momentum = final_momentum;
    end
    
    
    if ~isfield(options,'mom_switch_iter') || isempty(options.mom_switch_iter)
        options.mom_switch_iter = mom_switch_iter;
    end
    
    
    if ~isfield(options,'stop_lying_iter') || isempty(options.stop_lying_iter)
        options.stop_lying_iter = stop_lying_iter;
    end
    
    
    if ~isfield(options,'lie_multiplier') || isempty(options.lie_multiplier)
        options.lie_multiplier = lie_multiplier;
    end
    
    
    if ~isfield(options,'max_iter') || isempty(options.max_iter)
        options.max_iter = max_iter;
    end
    
    
    if ~isfield(options,'epsilon') || isempty(options.epsilon)
        options.epsilon = epsilon;
    end
    
    
    if ~isfield(options,'min_gain') || isempty(options.min_gain)
        options.min_gain = min_gain;
    end
    
    
    if ~isfield(options,'tsne_readout') || isempty(options.tsne_readout)
        options.tsne_readout = tsne_readout;
    end
    
    
    if ~isfield(options,'embedding_batchSize') || isempty(options.embedding_batchSize)
        options.embedding_batchSize = embedding_batchSize;
    end
    
    
    if ~isfield(options,'maxOptimIter') || isempty(options.maxOptimIter)
        options.maxOptimIter = maxOptimIter;
    end
    
    
    if ~isfield(options,'trainingSetSize') || isempty(options.trainingSetSize)
        options.trainingSetSize = trainingSetSize;
    end
    
    
    if ~isfield(options,'kdNeighbors') || isempty(options.kdNeighbors)
        options.kdNeighbors = kdNeighbors;
    end
    
    
    if ~isfield(options,'training_relTol') || isempty(options.training_relTol)
        options.training_relTol = training_relTol;
    end
    
    
    if ~isfield(options,'training_perplexity') || isempty(options.training_perplexity)
        options.training_perplexity = training_perplexity;
    end
    
    
    if ~isfield(options,'training_numPoints') || isempty(options.training_numPoints)
        options.training_numPoints = training_numPoints;
    end
    
    
    if ~isfield(options,'minTemplateLength') || isempty(options.minTemplateLength)
        options.minTemplateLength = minTemplateLength;
    end
    
    if ~isfield(options,'run_tsne') || isempty(options.run_tsne)
        options.run_tsne = run_tsne;
    end
    
    if ~isfield(options,'maxCarrierFrequency') || isempty(options.maxCarrierFrequency)
        options.maxCarrierFrequency = maxCarrierFrequency;
    end
    
    if ~isfield(options,'refine_clusters') || isempty(options.refine_clusters)
        options.refine_clusters = refine_clusters;
    end
    
    
    
    
    