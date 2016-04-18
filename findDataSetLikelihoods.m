function likelihoods = findDataSetLikelihoods(peakGuesses,template,sigmas)

    N = length(peakGuesses(:,1));
    L = length(sigmas);    
    
    likelihoods = zeros(N,1);
    for i=1:N 
        likelihoods(i) = -.5*sum((peakGuesses(i,:) - template(1:L)).^2./sigmas.^2 + log(2*pi.*sigmas.^2));
    end
    
    