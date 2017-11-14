function Y = normalizeHist(X,Y)

    Y = Y ./ (sum(Y)*(X(2)-X(1)));