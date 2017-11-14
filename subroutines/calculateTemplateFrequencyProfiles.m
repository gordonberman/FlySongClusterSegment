function [frequencyTemplates,fs,meanFreqs,stdFreqs,f_carriers,D] = calculateTemplateFrequencyProfiles(templates,Fs,groups)
%   Inputs:
%       templates -> Lx1 cell array of templates
%       Fs -> sampling frequency (in Hz)
%       groups -> Lx1 array that groups templates into clusters (i.e.
%                   groups = [1 1 2 1 3]; would group templates 1, 2, and 4 
%                   into one group, 3 into one, and 5 into another). Leave
%                   empty to ignore
%
%
%   Outputs:
%       frequencyTemplates -> Lx1 cell array of frequency templates
%       fs -> sampled frequencies (in Hz)
%       meanFreqs -> mean frequency profiles for each template
%       stdFreqs -> standard deviation of frequency profiles for each template
%       D -> distance matrix between all pairs of templates 
%               (potentially useful for distinguishing templates into 
%               clusters of templates, but unused right now)




    L = length(templates);
    
    if nargin >= 3 && ~isempty(groups)
        idx = unique(groups);
        L = length(idx);
        newTemplates = cell(L,1);
        for i=1:L
            newTemplates{i} = cell2mat(templates(groups == idx(i)));
        end
        templates = newTemplates;
    end
        
    
    
    d = length(templates{1}(1,:));
    N = floor(d/2);
    fs = linspace(0,Fs/2,N+1);
    fs = fs(2:end);

    frequencyTemplates = cell(L,1);
    for i=1:L
        frequencyTemplates{i} = zeros(length(templates{i}(:,1)),N);
        for j=1:length(templates{i}(:,1))
            q = fft(templates{i}(j,:));
            q = q.*conj(q);
            frequencyTemplates{i}(j,:) = q(2:N+1);
        end
    end
    
    
    meanFreqs = zeros(L,N);
    stdFreqs = zeros(L,N);
    f_carriers = zeros(L,1);
    for i=1:L
        meanFreqs(i,:) = mean(frequencyTemplates{i});
        stdFreqs(i,:) = std(frequencyTemplates{i});
        q = fit(fs',meanFreqs(i,:)','spline');
        idx = argmax(meanFreqs(i,:));
        f_carriers(i) = fminsearch(@(x) -q(x),fs(idx));
    end
    
    
    
    
    if nargout > 5
        D = zeros(L);
        for i=1:L
            for j=(i+1):L
                d = (meanFreqs(i,:) - meanFreqs(j,:)).^2./(stdFreqs(i,:).^2 + stdFreqs(j,:).^2);
                D(i,j) = sqrt(sum(d));
                D(j,i) = D(i,j);
            end
        end
    end