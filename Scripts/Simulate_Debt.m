[g_grd, Pg] = rouwen(gnum,log(g),rho,sigmasq);
g_grd=exp(g_grd);

%Generate shock sequences
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
%k(1)=k_LTC(1);
%tauKbar(1)=tauKbar_LTC(1);
%mu(1)=sigma_LTC(1);
%tauLbar(1)=tauLbar_LTC(1);

Xguess=[k(1:1000);tauKbar(1:1000);g_shocks(1:1000);mu(1:1000);tauLbar(1:1000);nu(1:1000);b(1:1000)];
Yguess=[k(2:1001);mu(2:1001);nu(2:1001);tauKbar(2:1001);tauLbar(2:1001); b(2:1001)];

netGuess=feedforwardnet(7);
netGuess.trainParam.showWindow = 0; 
netGuess=train(netGuess,Xguess,Yguess);

fval=zeros(6,T);
tic
for t=1:T
    
    state.k=k(t);
    state.tauKbar=tauKbar(t);
    state.g=g_shocks(t);
    state.g_iter=chain(t);
    state.mu=mu(t);
    state.tauLbar=tauLbar(t);
    state.b=b(t);
    state.nu=nu(t);
    
    X(:,t)=[k(t);tauKbar(t);g_shocks(t);mu(t);tauLbar(t);nu(t);b(t)];
    
    %% Solver
    if t<=1000
        init=[k(t+1),mu(t+1),nu(t+1),tauKbar(t+1),tauLbar(t+1), b(t+1)];
    else
        init=netGuess(X(:,t));%[k(t),mu(t),nu(t),tauKbar(t),tauLbar(t), b(t)];
    end
    
    phi=1500;
    
    L=-H;
    [x,fval(:,t)]=fsolve(@(x) FOCs2Debt(par,utility,state,x(1),x(2),x(3),x(4),x(5),x(6),netParams,g_grd,Pg,phi,L,H),init,opt);
    x = getX(par,utility,state,x(1),x(2),x(3),x(4),x(5),x(6),netParams,g_grd,Pg);
    c(t)=x(1);
    l(t)=x(2);
    k(t+1)=x(3);
    tauK(t)=x(4);
    tauL(t)=x(5);
    lambda(t)=x(6);
    mu(t+1)=x(7);
    nu(t+1)=x(8);
    xi(t)=x(9);
    b(t+1)=x(10);
    tauKbar(t+1)=x(11);
    tauLbar(t+1)=x(12);
    
    F(t)=Z*k(t)^alpha*l(t)^(1-alpha);
    Fk(t)=Z*alpha*k(t)^(alpha-1)*l(t)^(1-alpha);
    Fkk(t)=Z*alpha*(alpha-1)*k(t)^(alpha-2)*l(t)^(1-alpha);
    Fl(t)=Z*k(t)^alpha*(1-alpha)*l(t)^(-alpha);
    Fll(t)=Z*k(t)^alpha*(1-alpha)*(-alpha)*l(t)^(-alpha-1);
    Fkl(t)=Z*alpha*k(t)^(alpha-1)*(1-alpha)*l(t)^(-alpha);
    
    Eoutput=zeros(9,1);
    for g_iter_next=1:length(g_grd)
        Xtemp=[k(t+1);tauKbar(t+1);g_grd(g_iter_next);mu(t+1);tauLbar(t+1);nu(t+1);b(t+1)];
        Prediction=SimNetwork(Xtemp,netParams.ymaxin,netParams.yminin,netParams.xmaxin,netParams.xminin,netParams.ymaxout,netParams.yminout,netParams.xmaxout,netParams.xminout,netParams.IW,netParams.B1,netParams.LW,netParams.b2);
        Eoutput=Eoutput+Pg(state.g_iter,g_iter_next)*Prediction;
    end
    
    capitalTaxIncome(t)=tauK(t).*(Fk(t)-delta).*k(t);
    laborTaxIncome(t)=tauL(t).*Fl(t).*l(t);
    borrowedAmount(t)=beta*b(t+1).*Eoutput(6)./du(c(t))-b(t);
    
    Output1(t)=du(c(t)).*(1+(Fk(t)-delta)*(1-tauK(t)));
    Output2(t)=lambda(t).*(1+Fk(t)-delta);
    Output3(t)=du(c(t)).*Fkk(t).*(1-tauK(t));
    Output4(t)=nu(t+1).*du(c(t)).*(tauK(t).*(Fkk(t).*k(t)+Fk(t)-delta)+tauL(t).*Fkl(t).*l(t));
    Output5(t)=xi(t).*Fkl(t).*du(c(t)).*(1-tauL(t));
    Output6(t)=du(c(t));
    Output7(t)=nu(t+1).*du(c(t));
    Output8(t)=tauL(t);
    Output9(t)=tauK(t);
    
end