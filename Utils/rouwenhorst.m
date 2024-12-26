function [grid_rouv, P_rouv] = rouwenhorst(mu,rho,sigma_eps,numb_std,n)
v = numb_std*sigma_eps/(sqrt(1-rho^2));
nu = mu/(1-rho);
grid_rouv=linspace(nu-v,nu+v,n)';

% transition matrix  % need to check how to choose these p and q
qp = 1+rho;
p=0.5*qp;
q=p;
P= [p 1-p;1-q q];
P_rouv=P;
for i=3:length(grid_rouv)
    P_rouv = p*[P_rouv zeros(length(P_rouv),1);zeros(1,length(P_rouv)) 0]+...
        (1-p)*[zeros(length(P_rouv),1) P_rouv; 0 zeros(1,length(P_rouv))]+...
        (1-q)*[zeros(1,length(P_rouv),1) 0;P_rouv zeros(length(P_rouv),1)]+...
        q*[0 zeros(1,length(P_rouv)); zeros(length(P_rouv),1) P_rouv ];
    
    P_rouv(2:end-1,:) = 0.5*P_rouv(2:end-1,:);
end

grid_rouv = exp(grid_rouv);

end