function [amp,W] = fastWavelet_morlet_convolution_parallel(x,f,omega0,dt)


    N = length(x);
    L = length(f);
    amp = zeros(L,N);
    if mod(N,2) == 1
        x(end+1) = 0;
        N = N + 1;
        test = true;
    else
        test = false;
    end
    
    
    s = size(x);
    if s(2) == 1
        x = x';
    end
    
    x = [zeros(1,N/2) x zeros(1,N/2)];
    M = N;
    N = length(x);
    
    scales = (omega0 + sqrt(2+omega0^2))./(4*pi.*f);
    Omegavals = 2*pi*(-N/2:N/2-1)./(N*dt);
    
    xHat = fft(x);
    xHat = fftshift(xHat);
    
    if test
        idx = (M/2+1):(M/2+M-1);
    else
        idx = (M/2+1):(M/2+M);
    end
    
    if nargout == 2
        W = zeros(size(amp));
        test2 = true;
    else
        test2 = false;
    end
    
    for i=1:L
        
        m = morletConjFT(-Omegavals*scales(i),omega0);
        q  = ifft(m.*xHat)*sqrt(scales(i));
        
        q = q(idx);
        
        amp(i,:) = abs(q)*pi^-.25*exp(.25*(omega0-sqrt(omega0^2+2))^2)/sqrt(2*scales(i));
       
        if test2
            W(i,:) = q;
        end
    end
    
    