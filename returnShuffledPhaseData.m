function out = returnShuffledPhaseData(data)

    f = fft(data);
    a = real(f);
    b = imag(f);
    
    m = abs(f);
    p = atan2(b,a);
    
    p2 = p(randperm(length(p)));
    f2 = m.*cos(p2) + 1i*sin(p2);
    
    out = real(ifft(f2));