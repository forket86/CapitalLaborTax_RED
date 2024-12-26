function value = chebyAle(order,x,w,lowX,highX,lowT,highT,cross)

[row col]=size(x);
dimension=col;
if strcmp(class(x),'gpuArray')==1
    value=gpuArray(zeros(row,1));
else
    value=zeros(row,1);
end

if dimension==1
    x(:,1)=chebyTransformAle(lowX(1),highX(1),x(:,1),false);
    for i=0:order
        value(:)=value(:)+w(i+1,1)*chebyTAle(i,x(:,1));
    end
end

if dimension==2 && cross==false
    x(:,1)=chebyTransformAle(lowX(1),highX(1),x(:,1),false);
    x(:,2)=chebyTransformAle(lowX(2),highX(2),x(:,2),false);
    for i=0:order
        value(:)=value(:)+w(i+1,1)*chebyTAle(i,x(:,1))+w(i+1,2)*chebyTAle(i,x(:,2));
    end
end

if dimension==2 && cross==true
    x(:,1)=chebyTransformAle(lowX(1),highX(1),x(:,1),false);
    x(:,2)=chebyTransformAle(lowX(2),highX(2),x(:,2),false);
    for i=0:order
        value(:)=value(:)+w(i+1,1)*chebyTAle(i,x(:,1))+w(i+1,2)*chebyTAle(i,x(:,2))+w(i+1,3)*chebyTAle(i,x(:,1).*x(:,2));
    end
end

if dimension==3 && cross==false
    x(:,1)=chebyTransformAle(lowX(1),highX(1),x(:,1),false);
    x(:,2)=chebyTransformAle(lowX(2),highX(2),x(:,2),false);
    x(:,3)=chebyTransformAle(lowX(3),highX(3),x(:,3),false);
    for i=0:order
        value(:)=value(:)+w(i+1,1)*chebyTAle(i,x(:,1))+w(i+1,2)*chebyTAle(i,x(:,2))+w(i+1,3)*chebyTAle(i,x(:,3));
    end
end

if dimension==3 && cross==true
    x(:,1)=chebyTransformAle(lowX(1),highX(1),x(:,1),false);
    x(:,2)=chebyTransformAle(lowX(2),highX(2),x(:,2),false);
    x(:,3)=chebyTransformAle(lowX(3),highX(3),x(:,3),false);
    for i=0:order
        value(:)=value(:)+w(i+1,1)*chebyTAle(i,x(:,1))+w(i+1,2)*chebyTAle(i,x(:,2))+w(i+1,3)*chebyTAle(i,x(:,3))...
            +w(i+1,4)*chebyTAle(i,x(:,1).*x(:,2))+w(i+1,5)*chebyTAle(i,x(:,2).*x(:,3));
    end
end

if dimension==4 && cross==true
    x(:,1)=chebyTransformAle(lowX(1),highX(1),x(:,1),false);
    x(:,2)=chebyTransformAle(lowX(2),highX(2),x(:,2),false);
    x(:,3)=chebyTransformAle(lowX(3),highX(3),x(:,3),false);
    x(:,4)=chebyTransformAle(lowX(4),highX(4),x(:,4),false);
    for i=0:order
        value(:)=value(:)+w(i+1,1)*chebyTAle(i,x(:,1))+w(i+1,2)*chebyTAle(i,x(:,2))+w(i+1,3)*chebyTAle(i,x(:,3))+w(i+1,4)*chebyTAle(i,x(:,4))...
                         +w(i+1,5)*chebyTAle(i,x(:,1).*x(:,2))+w(i+1,6)*chebyTAle(i,x(:,2).*x(:,3))+w(i+1,7)*chebyTAle(i,x(:,3).*x(:,4))...
                         +w(i+1,8)*chebyTAle(i,x(:,1).*x(:,3))+w(i+1,9)*chebyTAle(i,x(:,1).*x(:,4))+w(i+1,10)*chebyTAle(i,x(:,2).*x(:,4))...
                         +w(i+1,11)*chebyTAle(i,x(:,1).*x(:,2).*x(:,3))+w(i+1,12)*chebyTAle(i,x(:,1).*x(:,2).*x(:,4))+w(i+1,13)*chebyTAle(i,x(:,2).*x(:,3).*x(:,4))...
                         +w(i+1,14)*chebyTAle(i,x(:,1).*x(:,2).*x(:,3).*x(:,4));
    end
end

if dimension==5 && cross==true
    x(:,1)=chebyTransformAle(lowX(1),highX(1),x(:,1),false);
    x(:,2)=chebyTransformAle(lowX(2),highX(2),x(:,2),false);
    x(:,3)=chebyTransformAle(lowX(3),highX(3),x(:,3),false);
    x(:,4)=chebyTransformAle(lowX(4),highX(4),x(:,4),false);
    x(:,5)=chebyTransformAle(lowX(5),highX(5),x(:,5),false);
    for i=0:order
        value(:)=value(:)+w(i+1,1)*chebyTAle(i,x(:,1))+w(i+1,2)*chebyTAle(i,x(:,2))+w(i+1,3)*chebyTAle(i,x(:,3))+w(i+1,4)*chebyTAle(i,x(:,4))...
                         +w(i+1,5)*chebyTAle(i,x(:,1).*x(:,2))+w(i+1,6)*chebyTAle(i,x(:,2).*x(:,3))+w(i+1,7)*chebyTAle(i,x(:,3).*x(:,4))...
                         +w(i+1,8)*chebyTAle(i,x(:,1).*x(:,3))+w(i+1,9)*chebyTAle(i,x(:,1).*x(:,4))+w(i+1,10)*chebyTAle(i,x(:,2).*x(:,4))...
                         +w(i+1,11)*chebyTAle(i,x(:,1).*x(:,2).*x(:,3))+w(i+1,12)*chebyTAle(i,x(:,1).*x(:,2).*x(:,4))+w(i+1,13)*chebyTAle(i,x(:,2).*x(:,3).*x(:,4))...
                         +w(i+1,14)*chebyTAle(i,x(:,1).*x(:,2).*x(:,3).*x(:,4));
                         %+w(i+1,15)*chebyTAle(i,x(:,5))+...
                         %+w(i+1,16)*chebyTAle(i,x(:,1).*x(:,5))+w(i+1,17)*chebyTAle(i,x(:,2).*x(:,5))+w(i+1,18)*chebyTAle(i,x(:,3).*x(:,5))+w(i+1,19)*chebyTAle(i,x(:,4).*x(:,5))...
                         %+w(i+1,20)*chebyTAle(i,x(:,1).*x(:,2).*x(:,5))+w(i+1,21)*chebyTAle(i,x(:,2).*x(:,3).*x(:,5))+w(i+1,22)*chebyTAle(i,x(:,3).*x(:,4).*x(:,5))...
                         %+w(i+1,23)*chebyTAle(i,x(:,1).*x(:,3).*x(:,5))+w(i+1,24)*chebyTAle(i,x(:,1).*x(:,4).*x(:,5))+w(i+1,25)*chebyTAle(i,x(:,2).*x(:,4).*x(:,5))...
                         %+w(i+1,26)*chebyTAle(i,x(:,1).*x(:,2).*x(:,3).*x(:,4).*x(:,5));
                     %+w(i+1,26)*chebyTAle(i,x(:,1).*x(:,2).*x(:,3).*x(:,5))+w(i+1,27)*chebyTAle(i,x(:,1).*x(:,2).*x(:,4).*x(:,5))+w(i+1,28)*chebyTAle(i,x(:,1).*x(:,3).*x(:,4).*x(:,5))+w(i+1,29)*chebyTAle(i,x(:,2).*x(:,3).*x(:,4).*x(:,5))...
                         
                         
    end
end

value(:)=chebyTransformAle(lowT,highT,value(:),true);

end

