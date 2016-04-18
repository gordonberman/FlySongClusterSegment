function vals = regionalmin(x)
    
    x = -x;

    idx = (x(2:end-1) > x(3:end)) & (x(2:end-1) > x(1:end-2));
    
    firstValue = x(1) > x(2);
    lastValue = x(end) > x(end-1);
    
    vals = [firstValue idx lastValue];
    