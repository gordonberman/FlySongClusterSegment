function options = makeDefaultOptions(options)

    %highPassFilterFreq -> Frequency for use in high-pass butterworth
    %                       filter of data (default = 200 Hz)
    
    %smoothSigma -> standard dev. of gaussian with which to smooth data,
    %               use -1 to have no smoothing (default = 5)
    
    %sigmaThreshold -> # of standard deviations above zero required for
    %                  peaks (default = 2)

    %diffThreshold -> half-width of peak window (default = 150)
    
    %noiseLevel -> noise level (default = .005)
    
    %minNoiseLevel -> minimum value for noise threshold (default = 0)
    
    %template_pca_dimension -> # of modes in clustering PCA (default = 50)
    
    %first_mode -> first mode in the PCA to uselikelihood model (default = 5)
    
    %k -> number of clusters in kmeans (default = 12)
    
    %maxNumPeaks -> maximum # of peaks to use in clustering (default = 10000);
    
    %kmean_replicates -> # of replicates in kmeans (default = 5)
    
    %kmeans_maxIter -> max # of iterations for kmeans (default = 1000)
        
    %fs = sampling frequency of original data
    
    %baseline_quantile = quantile of noise data to take as baseline value.
    %                    If set less than or equal to zero, there will be no 
    %                    baseline (default = .9)                      
    
    %min_template_size = minimum # of peaks in a template (default = 50)
    
    %baseline_threshold = likelihood threshold for baseline (default = 0)
    
    %use_likelihood_threshold = logical to use baseline threshold or not (default = false)
    
    
    highPassFilterFreq_default = 100;
    %smoothSigma_default = 3; 
    sigmaThreshold_default = 2;
    diffThreshold_default = 150;
    noiseLevel_default = -1;
    first_mode_default = 2;
    template_pca_dimension_default = 50;
    k_default = 12;
    maxNumPeaks_default = 10000;
    kmeans_replicates_default = 5;
    kmeans_maxIter_default = 1000;
    fs_default = 1e4;
    baseline_quantile_default = .9;
    min_template_size_default = 50;
    baseline_threshold_default = 0;
    use_likelihood_threshold_default = false;
    minNoiseLevel_default = 0;
    
    
    
    if options.setAll
        
        %options.smoothSigma = smoothSigma_default;
        options.sigmaThreshold = sigmaThreshold_default;
        options.diffThreshold = diffThreshold_default;
        options.noiseLevel = noiseLevel_default;
        options.first_mode = first_mode_default;
        options.k = k_default;        
        options.maxNumPeaks = maxNumPeaks_default;
        options.kmeans_replicates = kmeans_replicates_default;
        options.template_pca_dimension = template_pca_dimension_default;
        options.kmeans_maxIter = kmeans_maxIter_default;
        options.fs = fs_default;
        options.baseline_quantile = baseline_quantile_default;
        options.min_template_size = min_template_size_default;
        options.baseline_threshold = baseline_threshold_default;
        options.use_likelihood_threshold = use_likelihood_threshold_default;
        options.minNoiseLevel = minNoiseLevel_default;
        options.highPassFilterFreq = highPassFilterFreq_default;
               
    else
        
        %         if ~isfield(options,'smoothSigma') || isempty(options.smoothSigma)
        %             options.smoothSigma = smoothSigma_default;
        %         end
        
        
        if ~isfield(options,'highPassFilterFreq') || isempty(options.highPassFilterFreq)
            options.highPassFilterFreq = highPassFilterFreq_default;
        end
        
        if ~isfield(options,'sigmaThreshold') || isempty(options.sigmaThreshold)
            options.sigmaThreshold = sigmaThreshold_default;
        end
        
        
        if ~isfield(options,'minNoiseLevel') || isempty(options.minNoiseLevel)
            options.minNoiseLevel = minNoiseLevel_default;
        end
        
        
        if ~isfield(options,'diffThreshold') || isempty(options.diffThreshold)
            options.diffThreshold = diffThreshold_default;
        end
        
        
        if ~isfield(options,'noiseLevel') || isempty(options.noiseLevel)
            options.noiseLevel = noiseLevel_default;
        end
        
                        
        if ~isfield(options,'first_mode') || isempty(options.first_mode)
            options.first_mode = first_mode_default;
        end
        
        
        if ~isfield(options,'k') || isempty(options.k)
            options.k = k_default;
        end
        
        
        if ~isfield(options,'maxNumPeaks') || isempty(options.maxNumPeaks)
            options.maxNumPeaks = maxNumPeaks_default;
        end
        
        
        if ~isfield(options,'kmeans_replicates') || isempty(options.kmeans_replicates)
            options.kmeans_replicates = kmeans_replicates_default;
        end
        
        
        if ~isfield(options,'template_pca_dimension') || isempty(options.template_pca_dimension)
            options.template_pca_dimension = template_pca_dimension_default;
        end
        
        
        if ~isfield(options,'kmeans_maxIter') || isempty(options.kmeans_maxIter)
            options.kmeans_maxIter = kmeans_maxIter_default;
        end
        
        
        if ~isfield(options,'fs') || isempty(options.fs)
            options.fs = fs_default;
        end
        
                
        if ~isfield(options,'baseline_quantile') || isempty(options.baseline_quantile)
            options.baseline_quantile = baseline_quantile_default;
        end
        
       
        if ~isfield(options,'min_template_size') || isempty(options.min_template_size)
            options.min_template_size = min_template_size_default;
        end
        
        if ~isfield(options,'baseline_threshold') || isempty(options.baseline_threshold)
            options.baseline_threshold = baseline_threshold_default;
        end
        
        
        if ~isfield(options,'use_likelihood_threshold') || isempty(options.use_likelihood_threshold)
            options.use_likelihood_threshold = use_likelihood_threshold_default;
        end
        
    end
    
    
    
    