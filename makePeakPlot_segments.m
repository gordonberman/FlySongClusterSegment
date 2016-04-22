function makePeakPlot_segments(data,peakIdxGroup,toPlot,options)
    %makes peak plot (data, with lines, colored differently for different
    %templates)
    
    %Inputs:
    
    %data -> 1-D sounds trace
    %peakIdxGroup -> L x 1 cell array of peak values outputted from
    %                   createTemplates.m
    %toPlot -> 1-D array of the templates to plot (i.e. [1 3 4 5])

    if nargin < 4 || isempty(options);
        options.setAll = true;
    else
        options.setAll = false;
    end
    options = makeDefaultOptions(options);
    t = options.diffThreshold;
    
    %cs = 'r--k--g--m--c--r-.k-.g-.m-.c-.r-*k-*g-*m--c-*r-^k-^g-^m-^c-^';
    cs = [0 1 0;1 0 1;1 0 0;0 1 1;0 0 1;1 1 0;.5 0 1;0 .5 1;1 0 .5;0 1 .5;1 .5 0;.5 1 0];
    
    highVal = ceil(max(data)*2)/2;
    lowVal = floor(min(data)*2)/2;
    height = highVal - lowVal;
    
    L = length(toPlot);
    N = length(data);
    CCs = cell(L,1);
    for j=1:L
        
        isSignal = false(size(data));
        for i=1:length(peakIdxGroup{toPlot(j)})
            a = max([1, peakIdxGroup{toPlot(j)}(i)-t]);
            b = min([N, peakIdxGroup{toPlot(j)}(i)+t]);
            isSignal(a:b) = true;
        end
        
        CC = bwconncomp(isSignal);
        firstIdx = zeros(CC.NumObjects,1);
        lastIdx = zeros(CC.NumObjects,1);
        for i=1:CC.NumObjects
            firstIdx(i) = CC.PixelIdxList{i}(1);
            lastIdx(i) = CC.PixelIdxList{i}(end);
        end
        
        diffs = firstIdx(2:end) - lastIdx(1:end-1);
        idx = find(diffs <= 2*t + 1);
        for i=1:length(idx)
            isSignal(lastIdx(idx(i)):firstIdx(idx(i)+1)) = true;
        end
        
        CCs{j} = bwconncomp(isSignal);
        
    end

    hold on
    for i=1:L
        
        for j=1:CCs{i}.NumObjects
            x = CCs{i}.PixelIdxList{j}(1);
            y = CCs{i}.PixelIdxList{j}(end);
            rectangle('Position',[x lowVal y-x+1 height],'facecolor',cs(i,:),'edgecolor',cs(i,:));
        end
        
    end
    
    plot(data,'k-')
    
    
    
       %     plot(data,'b-');
    %     hold on
    %
    %     for i=1:L
    %         j = toPlot(i);
    %         for k=1:length(peakIdxGroup{j})
    %             plot([peakIdxGroup{j}(k) peakIdxGroup{j}(k)],[lowVal highVal],cs(3*i-2:3*i))
    %         end
    %     end
    %
    %     hold off
   