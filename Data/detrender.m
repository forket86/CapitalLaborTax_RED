function ydt = detrender(y,method)

switch method
    case 'off'
        ydt = y;
    case 'hp100'
        ydt = y - hpfilter(y,'Smoothing',100);
    case 'lin'
        T = length(y);
        X = [ones(T,1),(1:T)'];
        beta = X\y;
        ydt = y - X*beta;
    otherwise

end

