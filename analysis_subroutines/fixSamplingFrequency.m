function newdata = fixSamplingFrequency(data,oldfs,newfs)
    %subsample data to equalize the sampling frequency
    fprintf(['Sub-sampling (initial sampling frequency %d Hz) to ' ...
             'match sampling frequency of templates (%d Hz)\n'], ...
             oldfs,newfs);
    [P,Q] = rat(newfs/oldfs);
    newdata = resample(data,P,Q);
end