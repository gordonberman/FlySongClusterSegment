function out = morletConjFT(w,omega0)

    out = pi^(-1/4).*exp(-.5.*(w-omega0).^2);