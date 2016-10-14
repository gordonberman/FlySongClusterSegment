function [freqs,fs,meanFreqs,stdFreqs,f_carriers,D] = calculateTemplateFrequencyProfiles(templates,Fs)

    L = length(templates);
    d = length(templates{1}(1,:));
    N = floor(d/2);
    fs = linspace(0,Fs/2,N+1);
    fs = fs(2:end);

    freqs = cell(L,1);
    for i=1:L
        freqs{i} = zeros(length(templates{i}(:,1)),N);
        for j=1:length(templates{i}(:,1))
            q = fft(templates{i}(j,:));
            q = q.*conj(q);
            freqs{i}(j,:) = q(2:N+1);
        end
    end
    
    
    meanFreqs = zeros(L,N);
    stdFreqs = zeros(L,N);
    f_carriers = zeros(L,1);
    for i=1:L
        meanFreqs(i,:) = mean(freqs{i});
        stdFreqs(i,:) = std(freqs{i});
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