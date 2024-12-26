clear all;
clc;
close all;

addpath('../Utils/');

%% Data
warning off
load('./INIT/035_Farhi_Extreme')
warning on

rng(1)

%% Change cost
chi0=0.01;
psi0=10.0;
par.chi0=chi0;
par.psi0=psi0;

chi=@(tauL,tauLbar) (chi0/2)*(tauL-tauLbar).^2;    utility.chi=chi;
chi_tauL=@(tauL,tauLbar) chi0*(tauL-tauLbar);      utility.chi_tauL=chi_tauL;
chi_tauLbar=@(tauL,tauLbar) -chi0*(tauL-tauLbar);  utility.chi_tauLbar=chi_tauLbar;
psi=@(tauK,tauKbar) (psi0/2)*(tauK-tauKbar).^2;    utility.psi=psi;
psi_tauK=@(tauK,tauKbar) psi0*(tauK-tauKbar);      utility.psi_tauK=psi_tauK;
psi_tauKbar=@(tauK,tauKbar) -psi0*(tauK-tauKbar);  utility.psi_tauKbar=psi_tauKbar;

%% Let's Dance
netParams.ymaxin=net.inputs{1}.processSettings{1}.ymax;
netParams.yminin=net.inputs{1}.processSettings{1}.ymin;
netParams.xmaxin=net.inputs{1}.processSettings{1}.xmax;
netParams.xminin=net.inputs{1}.processSettings{1}.xmin;

netParams.ymaxout=net.outputs{2}.processSettings{1}.ymax;
netParams.yminout=net.outputs{2}.processSettings{1}.ymin;
netParams.xmaxout=net.outputs{2}.processSettings{1}.xmax;
netParams.xminout=net.outputs{2}.processSettings{1}.xmin;

netParams.IW = net.IW{1,1};
netParams.B1 = net.b{1};
netParams.LW = net.LW{2,1};
netParams.b2 = net.b{2};

k_prev=k;
c_prev=c;
labor_prev=l;
b_prev=b;
tauK_prev=tauK;
tauL_prev=tauL;

H=0.040;
for iter=1:1000
    
    iter
    
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
        if t<T
            init=[k(t+1),mu(t+1),nu(t+1),tauKbar(t+1),tauLbar(t+1), b(t+1)];
        else
            init=[k(t),mu(t),nu(t),tauKbar(t),tauLbar(t), b(t)];
        end
        
        phi=50;
        
        L=-H;
        [x,fval(:,t)]=fsolve(@(x) FOCs2(par,utility,state,x(1),x(2),x(3),x(4),x(5),x(6),netParams,g_grd,Pg,phi,L,H),init,opt);
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
    toc
    
    FOC_mean_residual=mean(fval,2)
    
    Poutput=SimNetwork(X(:,start:end),netParams.ymaxin,netParams.yminin,netParams.xmaxin,netParams.xminin,netParams.ymaxout,netParams.yminout,netParams.xmaxout,netParams.xminout,netParams.IW,netParams.B1,netParams.LW,netParams.b2);
    
    error1(iter)=mean((Output1(1,start:end)'-Poutput(1,:)').^2);
    error2(iter)=mean((Output2(1,start:end)'-Poutput(2,:)').^2);
    error3(iter)=mean((Output3(1,start:end)'-Poutput(3,:)').^2);
    error4(iter)=mean((Output4(1,start:end)'-Poutput(4,:)').^2);
    error5(iter)=mean((Output5(1,start:end)'-Poutput(5,:)').^2);
    error6(iter)=mean((Output6(1,start:end)'-Poutput(6,:)').^2);
    error7(iter)=mean((Output7(1,start:end)'-Poutput(7,:)').^2);
    error8(iter)=mean((Output8(1,start:end)'-Poutput(8,:)').^2);
    error9(iter)=mean((Output9(1,start:end)'-Poutput(9,:)').^2);
    
    errK(iter)=max(abs(k_prev(1:T)-k(1:T)));
    errC(iter)=max(abs(c_prev(:)-c(:)));
    errLabor(iter)=max(abs(labor_prev(:)-l(:)));
    errTauK(iter)=max(abs(tauK_prev(:)-tauK(:)));
    errTauB(iter)=max(abs(b_prev(:)-b(:)));
    errTauL(iter)=max(abs(tauL_prev(:)-tauL(:)));
    %errMainLoop(iter)=max([errK(iter) errC(iter) errLabor(iter) error1(iter) error2(iter) error3(iter) error4(iter) error5(iter) error6(iter) error7(iter) error8(iter) error9(iter)]);
    errMainLoop(iter)=max([error1(iter) error2(iter) error3(iter) error4(iter) error5(iter) error6(iter) error7(iter) error8(iter) error9(iter)]);
    errMainLoop(iter)
    %min(b)
    
    if errMainLoop(iter)<1e-4
       save(['./INIT/035_Farhi_Extreme_temp'])
       break;
    end
    
    k_prev=k;
    c_prev=c;
    labor_prev=l;
    b_prev=b;
    tauK_prev=tauK;
    tauL_prev=tauL;
    
    for a_iter=1:30
        [net,tr] = adapt(net,X(:,start:end),[Output1(:,start:end);Output2(:,start:end);Output3(:,start:end);Output4(:,start:end);Output5(:,start:end);Output6(:,start:end);Output7(:,start:end);Output8(:,start:end);Output9(:,start:end)]);
    end
    
    netParams.ymaxin=net.inputs{1}.processSettings{1}.ymax;
    netParams.yminin=net.inputs{1}.processSettings{1}.ymin;
    netParams.xmaxin=net.inputs{1}.processSettings{1}.xmax;
    netParams.xminin=net.inputs{1}.processSettings{1}.xmin;
    
    netParams.ymaxout=net.outputs{2}.processSettings{1}.ymax;
    netParams.yminout=net.outputs{2}.processSettings{1}.ymin;
    netParams.xmaxout=net.outputs{2}.processSettings{1}.xmax;
    netParams.xminout=net.outputs{2}.processSettings{1}.xmin;
    
    netParams.IW = net.IW{1,1};
    netParams.B1 = net.b{1};
    netParams.LW = net.LW{2,1};
    netParams.b2 = net.b{2};
    
    Poutput=SimNetwork(X,netParams.ymaxin,netParams.yminin,netParams.xmaxin,netParams.xminin,netParams.ymaxout,netParams.yminout,netParams.xmaxout,netParams.xminout,netParams.IW,netParams.B1,netParams.LW,netParams.b2);
    
end