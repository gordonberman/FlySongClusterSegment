function testMultipleNoisePostsSeparate(tempfile,valsignaldir,valnoisedir,...
    cattimes,outfile,np)
    %Run template assignment and false positive/negative estimation for a
    %range of noise_posterior_threshold values
    %outfile can be used as input for makePRPlot
    if nargin < 6 || isempty(np)
        np = (0:20)/20; %modify to allow user input
    end
    for n=1:length(np)
        SeparateSongsAssignTemplatesAndGetFPFN(tempfile,valsignaldir,...
            valnoisedir,cattimes,np(n),outfile,...
            'noise_posterior_threshold','valnp');
    end
end