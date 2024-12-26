function z0 = interpAle2(X,Y,Z,x0,y0)

for i=1:length(x0)
    Ytemp = interpAleMatrix(X,Z,x0(i));
    z0(i) = interpAle(Y,Ytemp,y0(i));
end

end

