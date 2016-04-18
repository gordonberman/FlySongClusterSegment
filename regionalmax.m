function vals = regionalmax(x)

    idx = (x(2:end-1) > x(3:end)) & (x(2:end-1) > x(1:end-2));
    
    firstValue = x(1) > x(2);
    lastValue = x(end) > x(end-1);

    s = size(idx);
    if s(2) > s(1)
        vals = [firstValue idx lastValue];
    else
        vals = [firstValue; idx; lastValue];
    end
    