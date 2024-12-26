function value = poly(order,x,phi)

[row col]=size(x);
dimension=col;

if dimension==1
    if order==0
    	value=[ones(row,1)]*phi';
    elseif order==1
    	value=[ones(row,1) x(:,1)]*phi';
    elseif order==2
        value=[ones(row,1) x(:,1) x(:,1).^2]*phi';
    elseif order==3
        value=[ones(row,1) x(:,1) x(:,1).^2 x(:,1).^3]*phi';
    end
end

if dimension==2
    if order==0
    	value=[ones(row,1)]*phi';
    elseif order==1
    	value=[ones(row,1) x(:,1) x(:,2)]*phi';
    elseif order==2
        value=[ones(row,1) x(:,1) x(:,1).^2 x(:,2) x(:,2).^2]*phi';
    elseif order==3
        value=[ones(row,1) x(:,1) x(:,1).^2 x(:,1).^3 x(:,2) x(:,2).^2 x(:,2).^3]*phi';
    end
end

if dimension==3
    if order==0
    	value=[ones(row,1)]*phi';
    elseif order==1
    	value=[ones(row,1) x(:,1) x(:,2) x(:,3)]*phi';
    elseif order==2
        value=[ones(row,1) x(:,1) x(:,1).^2 x(:,2) x(:,2).^2 x(:,3) x(:,3).^2]*phi';
    elseif order==3
        value=[ones(row,1) x(:,1) x(:,1).^2 x(:,1).^3 x(:,2) x(:,2).^2 x(:,2).^3 x(:,3) x(:,3).^2 x(:,3).^3]*phi';
    end
end

end