function y0 = interpAle(X,Y,x0)

for i=1:length(x0)
    m=1;
    q=1;
    
    if x0(i)<=X(1)
        m=(Y(2)-Y(1))/(X(2)-X(1));
        q=Y(1)-m*X(1);
    elseif x0(i)>X(1) && x0(i)<X(end)
        [val J]=min(abs(X-x0(i)));
        if X(J)>x0(i)
            J=J-1;
        end
        m=(Y(J+1)-Y(J))/(X(J+1)-X(J));
        q=Y(J)-m*X(J);
    elseif x0(i)>=X(end)
        m=(Y(end)-Y(end-1))/(X(end)-X(end-1));
        q=Y(end)-m*X(end);
    end
    
    y0(i)=m*x0(i)+q;
end

end

