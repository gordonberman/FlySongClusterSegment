function [peakIdxGroup,toPlot] = editPeaks_wm112216(data,peakIdxGroup,toPlot)

    plotRange = 50000;
    edgeValue = 10000;
    maxNewPeaks = 1000000;
    clickThreshold = 50;
    
    if nargin < 3 || isempty(toPlot)
        toPlot = 1:length(peakIdxGroup);
    end
    

    [peaks,idx] = sort(cell2mat(peakIdxGroup(toPlot)));
    peakAssignments = zeros(size(peaks));
    count = 0;
    for i=1:length(toPlot)
        L = length(peakIdxGroup{toPlot(i)});
        peakAssignments((1:L) + count) = toPlot(i);
        count = count + L;
    end
    peakAssignments = peakAssignments(idx);
    
    
    
    newPeaks = zeros(maxNewPeaks,1);
    currentPeak = max([peaks(1) 5001]);
    xlimits = round(currentPeak + [-edgeValue, plotRange-edgeValue]);
    xlimits(1) = max(1,xlimits(1)); %added 11/22/16
    
    
    h = figure(9942325.);
    count = 1;
    while xlimits(1) < peaks(end)
        
        firstPeakInFrame = find(peaks > xlimits(1),1,'first');
        lastPeakInFrame = find(peaks < xlimits(2),1,'last');
        peaksInFrame = peaks(firstPeakInFrame:lastPeakInFrame);
        peakAssignmentsInFrame = peakAssignments(firstPeakInFrame:lastPeakInFrame);
        currentNewPeaks = zeros(1000,1);
        currentCount = 1;
        keep = true(size(peaksInFrame));
        
        figure(9942325.)
        clf

        makePeakPlot(data(xlimits(1):xlimits(2)),{peaksInFrame(keep)},[],xlimits(1):xlimits(2));
        xlim([xlimits(1) xlimits(2)])
        title('Click to add or remove peak, Return/Enter to move to next screen, E to end','fontweight','bold','fontsize',12)
        
        
        button = 1;
        while ~isempty(button) && button ~= 69 && button ~= 101
            
            
            [x,~,button] = ginput(1);
            while ~isempty(button) && ~(button == 1 || button == 69 || button == 101)
                [x,~,button] = ginput(1);
            end
                        
            if button == 1
                %if a location is clicked
                
                x = round(x);
                [minD,minDIdx] = min(abs(peaksInFrame - x));
                
                if minD <= clickThreshold
                    
                    %if close to a peak in the original data set
                    keep(minDIdx) = false;
                    
                else
                    
                    if currentCount > 1
                        
                        [minD,minDIdx] = min(abs(currentNewPeaks(1:currentCount-1) - x));
                        if minD <= clickThreshold
                            
                            %if close to one of the new peaks
                            if currentCount > 2
                                currentNewPeaks(1:currentCount-2) = currentNewPeaks(setdiff(1:currentCount-1,minDIdx));
                                currentCount = currentCount - 1;
                            else
                                currentNewPeaks(1) = 0;
                                currentCount = currentCount - 1;
                            end
                            
                        else
                            
                            %if not close to any previous peak & previous
                            %new peaks chosen
                            currentNewPeaks(currentCount) = x;
                            currentCount = currentCount + 1;
                            
                        end
                        
                    else
                        
                        %if not close to any previous peak & no previous
                        %new peaks chosen
                        currentNewPeaks(currentCount) = x;
                        currentCount = currentCount + 1;
                        
                    end
                    
                end
                
            else
                %if the data collection is ended in some form

                if isempty(button)

                    idx = find(peaks >= xlimits(2),1,'first');
                    currentPeak = peaks(idx);
                    xlimitsNext = round(currentPeak + [-edgeValue, plotRange-edgeValue]);
                    if xlimitsNext(2) > length(data)
                        xlimitsNext = xlimitsNext + length(data) - xlimitsNext(2);
                    end

                else

                    xlimitsNext = peaks(end) + [1 100];

                end

            end
            
            
            figure(9942325.)
            clf
            if currentCount > 1
                makePeakPlot(data(xlimits(1):xlimits(2)),{[peaksInFrame(keep);currentNewPeaks(1:currentCount-1)]},1,xlimits(1):xlimits(2));
            else
                makePeakPlot(data(xlimits(1):xlimits(2)),{peaksInFrame(keep)},1,xlimits(1):xlimits(2));
            end
            xlim([xlimits(1) xlimits(2)])
            title('Click to add or remove peak, Return/Enter to move to next screen, E to end','fontweight','bold','fontsize',12)
            
            
        end
        
        
        xlimits = xlimitsNext;
        
        if currentCount > 1
            newPeaks(count + (1:(currentCount-1))) = currentNewPeaks(1:(currentCount-1));
            count = count + currentCount - 1;
        end
        

        if sum(keep) < length(keep)
            idx = find(~keep);
            assignmentVals = peakAssignmentsInFrame(~keep);
            q = unique(assignmentVals);
            for ii=1:length(q)
                j = q(ii);
                idx2 = idx(assignmentVals == j);                
                peakIdxGroup{j} = setdiff(peakIdxGroup{j},peaksInFrame(idx2));
            end
        end
        
        
    end
    
    
    
    if count > 1
        peakIdxGroup = [peakIdxGroup; newPeaks(1:count-1)];
        toPlot = [toPlot;length(toPlot)+1];
    end
    
    for i=1:length(peakIdxGroup)
        peakIdxGroup{i} = sort(peakIdxGroup{i});
    end
    
    figure(9942325.)
    clf
    makePeakPlot(data,peakIdxGroup,toPlot);
    
