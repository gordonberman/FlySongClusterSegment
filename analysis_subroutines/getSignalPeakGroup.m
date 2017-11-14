function signalPeakGroup = getSignalPeakGroup(signalPeakIdx, peakIdxGroup, ...
    isNoise)
    %Get templateGroupings for all signal peaks.
    %Modified from getSignalPeakGroup_unclustered. Now assumes peakIdxGroup
    %has length (# of groups) rather than (# of templates).
    signalPeakGroup = zeros(1,length(signalPeakIdx));
    for p = 1:length(peakIdxGroup)
        if isNoise(p) == 0
            spiind = find(ismember(signalPeakIdx,peakIdxGroup{p}));
            %signalPeakGroup(spiind) = templateGroupings(p);
            signalPeakGroup(spiind) = p;
        end
    end
end