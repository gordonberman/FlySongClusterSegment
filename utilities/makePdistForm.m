function DD = makePdistForm(D)

    L = length(D(:,1));
    DD = zeros(1,L*(L-1)/2);
    count = 1;
    for i=1:(L-1)
        
        q = D((i+1):end,i);
        DD(count:count+length(q)-1) = q';
        count = count + length(q);
        
    end