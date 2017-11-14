function signalPeakIdx = getSignalPeakIdx(peakIdxGroup, isNoise)
    %Get peakIdx for all signal peaks.
    whichsignal = find(isNoise==0)';
    signalPeakIdx = {};
    if length(whichsignal) > length(peakIdxGroup)
        loopover = 1:(length(peakIdxGroup)-1);
    else
        loopover = whichsignal;
    end
    for i = 1:length(loopover)
        if isa(signalPeakIdx,'cell')
            signalPeakIdx = cell2mat(signalPeakIdx);
        end
        %signalPeakIdx = horzcat(signalPeakIdx,peakIdxGroup{whichsignal(i)});
        try
            signalPeakIdx = vertcat(signalPeakIdx,peakIdxGroup{loopover(i)});
        catch MEmatdim
            error(['Trying to get peakIdxGroup %.0f, but length is only' ...
                ' %.0f.\n Check whether you have templateGroupings (this ' ...
                'reduces the length of peakIdxGroup).\n'],...
                loopover(i),length(peakIdxGroup));
            rethrow(MEmatdim);
        end
    end
    signalPeakIdx = sort(signalPeakIdx)'; %to retain dimension from earlier version
end