function r = rebuildVector(grid, vector, symbol, minValue, maxValue)

goodValuesIndeces=find(vector~=symbol);
r=vector;

if length(goodValuesIndeces)==0
    return;
elseif length(goodValuesIndeces)==1
    for iter=1:length(vector)
        if vector(iter)==symbol
            if iter<goodValuesIndeces
                r(iter)=minValue;
            else
                r(iter)=maxValue;
            end
        end
    end
    return;
end

for iter=1:length(goodValuesIndeces)
    gridTemp(iter)=grid(goodValuesIndeces(iter));
    vectorTemp(iter)=vector(goodValuesIndeces(iter));
end


for iter=1:length(vector)
    if vector(iter)==symbol
        r(iter)=interpAle(gridTemp,vectorTemp,grid(iter));
    end
end

end

