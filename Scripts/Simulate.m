[g_grd, Pg] = rouwen(gnum,log(g),rho,sigmasq);
g_grd=exp(g_grd);
rng(1)
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

%% Simulate existing solution again
tauKbar_LTC(1)=tauK_LTC(1);
for t=1:T
    tauL_LTC(t)=tauLtilde(k_LTC(t),tauK_LTC(t),g_shocks(t),sigma_LTC(t),tauLbar_LTC(t),phitauL_LTC);
    sigma_LTC(t+1)=sigmatilde(k_LTC(t),tauK_LTC(t),g_shocks(t),sigma_LTC(t),tauLbar_LTC(t),phisigma_LTC);
    [labor_LTC(t),c_LTC(t),k_LTC(t+1),labor_LTC_tauL(t),c_LTC_tauL(t),labor_LTC_tauK(t),c_LTC_tauK(t),labor_LTC_k(t),c_LTC_k(t)]=...
        getLaborCKnext(k_LTC(t),tauK_LTC(t),tauL_LTC(t),Z,g_shocks(t),alpha,delta,gamma,eta,D);
    
    lambda_LTC(t)=(du(c_LTC(t))-(sigma_LTC(t+1)-sigma_LTC(t)).*ddu(c_LTC(t)).*k_LTC(t+1)+sigma_LTC(t).*(du(c_LTC(t))+ddu(c_LTC(t)).*c_LTC(t))-chi_tauL(tauL_LTC(t),tauLbar_LTC(t))./c_LTC_tauL(t)+((-dv(labor_LTC(t))+sigma_LTC(t).*(-dv(labor_LTC(t))-ddv(labor_LTC(t)).*labor_LTC(t)))).*labor_LTC_tauL(t)./c_LTC_tauL(t))./(1-(1-alpha).*Z.*k_LTC(t).^(alpha).*labor_LTC(t).^(-alpha).*labor_LTC_tauL(t)./c_LTC_tauL(t));
    nu_LTC(t) = du(c_LTC(t))-(sigma_LTC(t+1)-sigma_LTC(t)).*ddu(c_LTC(t)).*k_LTC(t+1)+sigma_LTC(t).*(du(c_LTC(t))+ddu(c_LTC(t)).*c_LTC(t))-lambda_LTC(t);
    mu_LTC(t)=(nu_LTC(t).*c_LTC_tauL(t)-chi_tauL(tauL_LTC(t),tauLbar_LTC(t)))./labor_LTC_tauL(t);
    
    tauLbar_LTC(t+1)=tauLbartilde(k_LTC(t),tauK_LTC(t),g_shocks(t),sigma_LTC(t),tauLbar_LTC(t),phitauLbar_LTC);
    tauK_LTC(t+1)=tauKtilde(k_LTC(t),tauK_LTC(t),g_shocks(t),sigma_LTC(t),tauLbar_LTC(t),phitauK_LTC);
    tauKbar_LTC(t+1)=tauK_LTC(t+1);
    
    ERHS=0;
    EOutput1(t)=0;
    EOutput2(t)=0;
    EOutput3(t)=0;
    EOutput4(t)=0;
    EOutput5(t)=0;
    knext=k_LTC(t+1);
    tauKnext=tauK_LTC(t+1);
    sigmanext=sigma_LTC(t+1);
    tauLbarnext=tauLbar_LTC(t+1);
    z=Z;
    for g_iter_next=1:length(grids.g_grd)
        g_next=grids.g_grd(g_iter_next);
        tauLnext=tauLtilde(k_LTC(t+1),tauK_LTC(t+1),g_next,sigma_LTC(t+1),tauLbar_LTC(t+1),phitauL_LTC);
        [labornext,cnext,knextnext]=getLaborCKnext(k_LTC(t+1),tauK_LTC(t+1),tauLnext,Z,g_next,alpha,delta,gamma,eta,D);
        
        RHS=du(cnext).*(1+(alpha*grids.z_grd.*k_LTC(t+1).^(alpha-1).*labornext.^(1-alpha)-delta).*(1-tauK_LTC(t+1)));
        ERHS=ERHS+Pg(chain(t),g_iter_next).*RHS;
        
        [lnext,cnext,knextnext,lnext_tauL,cnext_tauL,lnext_tauK,cnext_tauK,lnext_k,cnext_k] = getLaborCKnext(knext,tauKnext,tauLnext,z,g_next,alpha,delta,gamma,eta,D);
        sigmanextnext=sigmatilde(knext,tauKnext,g_next,sigmanext,tauLbar_LTC(t+1),phisigma_LTC);
        lambdanext=(du(cnext)-(sigmanextnext-sigmanext).*ddu(cnext).*knextnext+sigmanext.*(du(cnext)+ddu(cnext).*cnext)-chi_tauL(tauLnext,tauLbarnext)./cnext_tauL+((-dv(lnext)+sigmanext.*(-dv(lnext)+-ddv(lnext).*lnext))).*lnext_tauL./cnext_tauL)./(1-(1-alpha).*z.*knext.^(alpha).*lnext.^(-alpha).*lnext_tauL./cnext_tauL);
        nunext = du(cnext)-(sigmanextnext-sigmanext).*ddu(cnext).*knextnext+sigmanext.*(du(cnext)+ddu(cnext).*cnext)-lambdanext;
        % Eq.17
        munext=(nunext.*cnext_tauL-chi_tauL(tauLnext,tauLbarnext))./lnext_tauL;
        
        EOutput1(t)=EOutput1(t)+(du(cnext).*(cnext+knextnext)-dv(labornext).*labornext).*Pg(chain(t),g_iter_next);
        EOutput2(t)=EOutput2(t)+lambdanext.*(Z.*alpha.*knext.^(alpha-1).*labornext.^(1-alpha)+1-delta).*Pg(chain(t),g_iter_next);
        EOutput3(t)=EOutput3(t)+(munext.*lnext_k-nunext.*cnext_k).*Pg(chain(t),g_iter_next);
        EOutput4(t)=EOutput4(t)+tauLnext.*Pg(chain(t),g_iter_next);
        EOutput5(t)=EOutput5(t)+tauK_LTC(t+1).*Pg(chain(t),g_iter_next);
    end
    
    err(t)=du(c_LTC(t))-beta*ERHS;
    
    X(:,t)=[k_LTC(t);tauKbar_LTC(t);g_shocks(t);sigma_LTC(t);tauLbar_LTC(t)];
    Output1(t)=(du(c_LTC(t)).*(c_LTC(t)+knext)-dv(labor_LTC(t)).*labor_LTC(t));
    Output2(t)=lambda_LTC(t).*(Z.*alpha.*k_LTC(t).^(alpha-1).*labor_LTC(t).^(1-alpha)+1-delta);
    Output3(t)=(mu_LTC(t).*labor_LTC_k(t)-nu_LTC(t).*c_LTC_k(t));
    Output4(t)=tauL_LTC(t);
    Output5(t)=tauK_LTC(t);
    
    errOutput1(t)=du(c_LTC(t)).*k_LTC(t+1)-beta*EOutput1(t);
    errOutput2(t)=lambda_LTC(t)-beta*EOutput2(t)+du(c_LTC(t)).*(sigmanext-sigma_LTC(t))+beta*EOutput3(t);
    errOutput4(t)=tauLbar_LTC(t+1)-EOutput4(t);
    errOutput5(t)=tauKbar_LTC(t+1)-EOutput5(t);
    
end

clear k tauKbar tauLbar sigma tauL tauK labor lambda mu nu c labor_tauL c_tauL labor_tauK c_tauK labor_k c_k
k(1)=k_LTC(1);
tauKbar(1)=tauKbar_LTC(1);
tauLbar(1)=tauLbar_LTC(1);
sigma(1)=sigma_LTC(1);

fval=NaN(5,T);
tic
for t=1:T
    
    state.k=k(t);
    state.tauKbar=tauKbar(t);
    state.g=g_shocks(t);
    state.g_iter=chain(t);
    state.sigma=sigma(t);
    state.tauLbar=tauLbar(t);
    
    X(:,t)=[k(t);tauKbar(t);g_shocks(t);sigma(t);tauLbar(t)];
    
    %% Solver
    tauK0=tauKbar(t);
    tauL0=tauLFC;
    sigmanext0=sigmaFC;
    tauLbarnext0=tauLbar(t);
    tauKbarnext0=tauKbar(t);
    if t>1
        tauK0=tauK(t-1);
        tauL0=tauL(t-1);
        sigmanext0=sigma(t-1);
    end
    
    %
    init=[tauL_LTC(t),tauK_LTC(t),sigma_LTC(t+1),tauLbar_LTC(t+1),tauKbar_LTC(t+1)];
    
    [x,fval(:,t)]=fsolve(@(x) FOCs2(par,utility,state,x(1),x(2),x(3),x(4),x(5),Approximator,g_grd,Pg),init,opt);
    %[x,fval(:,t)]=fminunc(@(x) sum(FOCs2(par,utility,state,x(1),x(2),x(3),x(4),x(5),Approximator,g_grd,Pg).^2),init,options2);
    
    tauL(t)     = x(1);
    tauK(t)     = x(2);
    sigma(t+1)  = x(3);
    tauLbar(t+1)= x(4);
    tauKbar(t+1)= x(5);
    
    %FOCs2(par,utility,state,x(1),x(2),x(3),x(4),x(5),netParams,g_grd,Pg)
    
    [labor(t),c(t),k(t+1),labor_tauL(t),c_tauL(t),labor_tauK(t),c_tauK(t),labor_k(t),c_k(t)]=...
        getLaborCKnext(k(t),tauK(t),tauL(t),Z,g_shocks(t),alpha,delta,gamma,eta,D);
    
    lambda(t)=(du(c(t))-(sigma(t+1)-sigma(t)).*ddu(c(t)).*k(t+1)+sigma(t).*(du(c(t))+ddu(c(t)).*c(t))-chi_tauL(tauL(t),tauLbar(t))./c_tauL(t)+((-dv(labor(t))+sigma(t).*(-dv(labor(t))-ddv(labor(t)).*labor(t)))).*labor_tauL(t)./c_tauL(t))./(1-(1-alpha).*Z.*k(t).^(alpha).*labor(t).^(-alpha).*labor_tauL(t)./c_tauL(t));
    nu(t) = du(c(t))-(sigma(t+1)-sigma(t)).*ddu(c(t)).*k(t+1)+sigma(t).*(du(c(t))+ddu(c(t)).*c(t))-lambda(t);
    mu(t)=(nu(t).*c_tauL(t)-chi_tauL(tauL(t),tauLbar(t)))./labor_tauL(t);
    
end