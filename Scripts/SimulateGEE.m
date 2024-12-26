rng(1)

[g_grd, Pg] = rouwen(gnum,log(g),rho,sigmasq);
g_grd=exp(g_grd);

%Generate shock sequences of TFP
chain = zeros(1,T);
chain(1)=1;
g_shocks=zeros(1,length(chain));
g_shocks_indeces=zeros(1,length(chain));
g_shocks(1)= g_grd(chain(1));
g_shocks_indeces(1)= chain(1);
for i=2:T
    this_step_distribution = Pg(chain(i-1),:);
    cumulative_distribution = cumsum(this_step_distribution);
    r = rand();
    chain(i) = find(cumulative_distribution>r,1);
    g_shocks_indeces(i)=chain(i);
    g_shocks(i)= g_grd(chain(i));
end

fval=NaN(5,T);
tic
for t=1:T
    
    state.k=k(t);
    state.tauKbar=tauKbar(t);
    state.g=g_shocks(t);
    state.g_iter=chain(t);
    state.tauLbar=tauLbar(t);
    
    X(:,t)=[k(t);tauKbar(t);g_shocks(t);tauLbar(t)];
    
    %% Solver
    
    %tauK0=tauKinit(t);
    %tauL0=tauLinit(t);
    %sigmanext0=sigmainit(t+1);
    %tauKbarnext0=tauKbarinit(t+1);
    %tauLbarnext0=tauLbarinit(t+1);
    
    %if iter>1 && t>1
    tauK0=tauK(t);
    tauL0=tauL(t);
    sigmanext0=sigma(t+1);
    tauKbarnext0=tauKbar(t+1);
    tauLbarnext0=tauLbar(t+1);
    %end
    
    
    [x,fval(:,t)]=fsolve(@(x) FOCs5(par,utility,state,x(1),x(2),x(3),x(4),x(5),Approximator,g_grd,Pg,tauKmax,penalty),[tauL0,tauK0,sigmanext0,tauLbarnext0,tauKbarnext0],opt);
    tauL(t)     = x(1);
    tauK(t)     = x(2);
    sigma(t+1)  = x(3);
    tauLbar(t+1)= x(4);
    tauKbar(t+1)= x(5);
    
    [labor(t),c(t),k(t+1),labor_tauL(t),c_tauL(t),labor_tauK(t),c_tauK(t),labor_k(t),c_k(t)]=...
        getLaborCKnext(k(t),tauK(t),tauL(t),Z,g_shocks(t),alpha,delta,gamma,eta,D);
    
    lambda(t)=(du(c(t))-sigma(t+1).*ddu(c(t)).*k(t+1)-chi_tauL(tauL(t),tauLbar(t))./c_tauL(t)+(-dv(labor(t))).*labor_tauL(t)./c_tauL(t))./(1-(1-alpha).*Z.*k(t).^(alpha).*labor(t).^(-alpha).*labor_tauL(t)./c_tauL(t));
    nu(t) = du(c(t))-sigma(t+1).*ddu(c(t)).*k(t+1)-lambda(t);
    mu(t)=(nu(t).*c_tauL(t)-chi_tauL(tauL(t),tauLbar(t)))./labor_tauL(t);
    %check=dv(labor(t))-lambda(t)*(1-alpha).*Z.*k(t).^(alpha).*labor(t).^(-alpha)-mu(t);
    
    %r(t)=alpha*z*k(t).^(alpha-1).*labor(t).^(1-alpha);
    %wage(t)=(1-alpha).*Z.*k(t)^alpha.*labor(t).^(-alpha);
    %check_GovBudget(t)= tauK(t)*(r(t)-delta).*k(t)+tauL(t).*wage(t).*labor(t)-g_shocks(t);
    
    Output1(t)=(du(c(t)).*(c(t)+k(t+1))-dv(labor(t)).*labor(t));
    Output2(t)=lambda(t).*(Z.*alpha.*k(t).^(alpha-1).*labor(t).^(1-alpha)+1-delta);
    Output3(t)=(mu(t).*labor_k(t)-nu(t).*c_k(t));
    Output4(t)=tauL(t);
    Output5(t)=tauK(t);
    
end