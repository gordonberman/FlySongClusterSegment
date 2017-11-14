function [signalPeakIdx, signalFreqIdx] = getSignalFreqIdx(peakIdxGroup,...
    freqIdxGroup, isNoise)
    %Get peakIdx and freqIdx for all signal peaks.
    %Needed: way to test that dimensions of peakIdxGroup and freqIdxGroup
    %are the same
    whichsignal = find(isNoise==0)';
    signalPeakIdx = {};
    signalFreqIdx = {};
    if length(whichsignal) > length(peakIdxGroup)
        loopover = 1:(length(peakIdxGroup)-1);
    else
        loopover = whichsignal;
    end
    for i = 1:length(loopover)
        if isa(signalFreqIdx,'cell')
            signalFreqIdx = cell2mat(signalFreqIdx);
        end
        %signalPeakIdx = horzcat(signalPeakIdx,peakIdxGroup{whichsignal(i)});
        try
            signalFreqIdx = vertcat(signalFreqIdx,freqIdxGroup{loopover(i)});
        catch MEmatdim
            error(['Trying to get freqIdxGroup %.0f, but length is only' ...
                ' %.0f.\n Check whether you have templateGroupings (this ' ...
                'reduces the length of freqIdxGroup).\n'],...
                loopover(i),length(freqIdxGroup));
            rethrow(MEmatdim);
        end
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
    %sort freqs using same sort order as peaks
    [signalPeakIdx, sortidx] = sort(signalPeakIdx'); %to retain dimension from earlier version
    signalFreqIdx = signalFreqIdx(sortidx)';
end
