clear all;
clc;
close all;

addpath('../Utils/');

%% Data
load(['./INIT/LowPers39_39'])

rng(1)

cutInit=50;

neural=true;
T=1000;
adaptNum=6;
smoothadapt=true;
options2 = optimoptions('fminunc','Display','off');
chi0_grd=39;
for chi0_iter=1:length(chi0_grd)
    chi0=chi0_grd(chi0_iter);
    display(['************** chi_0=' num2str(chi0) ' ********************'])

    %% Parameters
    psi0=39;
    %chi0=0.01;
    par.chi0=chi0;
    par.psi0=psi0;

    chi=@(tauL,tauLbar) (chi0/2)*(tauL-tauLbar).^2;    utility.chi=chi;
    chi_tauL=@(tauL,tauLbar) chi0*(tauL-tauLbar);      utility.chi_tauL=chi_tauL;
    chi_tauLbar=@(tauL,tauLbar) -chi0*(tauL-tauLbar);  utility.chi_tauLbar=chi_tauLbar;
    psi=@(tauK,tauKbar) (psi0/2)*(tauK-tauKbar).^2;    utility.psi=psi;
    psi_tauK=@(tauK,tauKbar) psi0*(tauK-tauKbar);      utility.psi_tauK=psi_tauK;
    psi_tauKbar=@(tauK,tauKbar) -psi0*(tauK-tauKbar);  utility.psi_tauKbar=psi_tauKbar;


    %% Time simulation
    opt = optimoptions('fsolve','Display','off','Algorithm','Levenberg-Marquardt','MaxIterations',10000,'MaxFunctionEvaluation',10000);%,'FunctionTolerance',1e-9,'OptimalityTolerance',1e-14);

    Approximator.neural=neural;
    if neural==false
        const=1;
        Approximator.const=const;
        Approximator.poly=@(x,beta) exp(beta(1)+(beta(2:end)')*log(1+x))-const;

        Approximator.beta1=regress(log(const+Output1(:,start:T))',[ones(length(Output1(:,start:T)),1),log(const+X(:,start:T)')]);
        Approximator.beta2=regress(log(const+Output2(:,start:T))',[ones(length(Output2(:,start:T)),1),log(const+X(:,start:T)')]);
        Approximator.beta3=regress(log(const+Output3(:,start:T))',[ones(length(Output3(:,start:T)),1),log(const+X(:,start:T)')]);
        Approximator.beta4=regress(log(const+Output4(:,start:T))',[ones(length(Output4(:,start:T)),1),log(const+X(:,start:T)')]);
        Approximator.beta5=regress(log(const+Output5(:,start:T))',[ones(length(Output5(:,start:T)),1),log(const+X(:,start:T)')]);

    else


        Approximator.ymaxin=net.inputs{1}.processSettings{1}.ymax;
        Approximator.yminin=net.inputs{1}.processSettings{1}.ymin;
        Approximator.xmaxin=net.inputs{1}.processSettings{1}.xmax;
        Approximator.xminin=net.inputs{1}.processSettings{1}.xmin;

        Approximator.ymaxout=net.outputs{2}.processSettings{1}.ymax;
        Approximator.yminout=net.outputs{2}.processSettings{1}.ymin;
        Approximator.xmaxout=net.outputs{2}.processSettings{1}.xmax;
        Approximator.xminout=net.outputs{2}.processSettings{1}.xmin;

        Approximator.IW = net.IW{1,1};
        Approximator.B1 = net.b{1};
        Approximator.LW = net.LW{2,1};
        Approximator.b2 = net.b{2};
    end
    for iter=1:100000

        iter

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

            init=[tauL_LTC(t),tauK_LTC(t),sigma_LTC(t+1),tauLbar_LTC(t+1),tauKbar_LTC(t+1)];

            [x,fval(:,t)]=fsolve(@(x) FOCs2(par,utility,state,x(1),x(2),x(3),x(4),x(5),Approximator,g_grd,Pg),init,opt);

            tauL(t)     = x(1);
            tauK(t)     = x(2);
            sigma(t+1)  = x(3);
            tauLbar(t+1)= x(4);
            tauKbar(t+1)= x(5);


            [labor(t),c(t),k(t+1),labor_tauL(t),c_tauL(t),labor_tauK(t),c_tauK(t),labor_k(t),c_k(t)]=...
                getLaborCKnext(k(t),tauK(t),tauL(t),Z,g_shocks(t),alpha,delta,gamma,eta,D);

            lambda(t)=(du(c(t))-(sigma(t+1)-sigma(t)).*ddu(c(t)).*k(t+1)+sigma(t).*(du(c(t))+ddu(c(t)).*c(t))-chi_tauL(tauL(t),tauLbar(t))./c_tauL(t)+((-dv(labor(t))+sigma(t).*(-dv(labor(t))-ddv(labor(t)).*labor(t)))).*labor_tauL(t)./c_tauL(t))./(1-(1-alpha).*Z.*k(t).^(alpha).*labor(t).^(-alpha).*labor_tauL(t)./c_tauL(t));
            nu(t) = du(c(t))-(sigma(t+1)-sigma(t)).*ddu(c(t)).*k(t+1)+sigma(t).*(du(c(t))+ddu(c(t)).*c(t))-lambda(t);
            mu(t)=(nu(t).*c_tauL(t)-chi_tauL(tauL(t),tauLbar(t)))./labor_tauL(t);

            Output1(t)=(du(c(t)).*(c(t)+k(t+1))-dv(labor(t)).*labor(t));
            Output2(t)=lambda(t).*(Z.*alpha.*k(t).^(alpha-1).*labor(t).^(1-alpha)+1-delta);
            Output3(t)=(mu(t).*labor_k(t)-nu(t).*c_k(t));
            Output4(t)=tauL(t);
            Output5(t)=tauK(t);

        end
        toc

        FOC_mean_residual=mean(fval,2)


        Poutput=SimNetwork(X(:,start:T),Approximator.ymaxin,Approximator.yminin,Approximator.xmaxin,Approximator.xminin,Approximator.ymaxout,Approximator.yminout,Approximator.xmaxout,Approximator.xminout,Approximator.IW,Approximator.B1,Approximator.LW,Approximator.b2);

        error1(iter)=mean((Output1(1,start:T)'-Poutput(1,:)').^2);
        error2(iter)=mean((Output2(1,start:T)'-Poutput(2,:)').^2);
        error3(iter)=mean((Output3(1,start:T)'-Poutput(3,:)').^2);
        error4(iter)=mean((Output4(1,start:T)'-Poutput(4,:)').^2);
        error5(iter)=mean((Output5(1,start:T)'-Poutput(5,:)').^2);

        errK(iter)=max(abs(k_prev(:)-k(:)));
        errC(iter)=max(abs(c_prev(:)-c(:)));
        errLabor(iter)=max(abs(labor_prev(:)-labor(:)));
        errTauK(iter)=max(abs(tauK_prev(:)-tauK(:)));
        errTauL(iter)=max(abs(tauL_prev(:)-tauL(:)));
        errMainLoop(iter)=max([errK(iter) errC(iter) errLabor(iter) errTauK(iter) errTauL(iter) error1(iter) error2(iter) error3(iter) error4(iter) error5(iter)]);
        errMainLoop(iter)
        if errMainLoop(iter)<1e-3
            save(['./INIT/LowPers39_39_temp'])
            break;
        end

        k_prev=k;
        c_prev=c;
        labor_prev=labor;
        tauK_prev=tauK;
        tauL_prev=tauL;

        if neural==false
            Approximator.beta1=regress(log(const+Output1(:,start:T))',[ones(length(Output1(:,start:T)),1),log(const+X(:,start:T)')]);
            Approximator.beta2=regress(log(const+Output2(:,start:T))',[ones(length(Output2(:,start:T)),1),log(const+X(:,start:T)')]);
            Approximator.beta3=regress(log(const+Output3(:,start:T))',[ones(length(Output3(:,start:T)),1),log(const+X(:,start:T)')]);
            Approximator.beta4=regress(log(const+Output4(:,start:T))',[ones(length(Output4(:,start:T)),1),log(const+X(:,start:T)')]);
            Approximator.beta5=regress(log(const+Output5(:,start:T))',[ones(length(Output5(:,start:T)),1),log(const+X(:,start:T)')]);
        else


            for a_iter=1:1
                [net,tr] = train(net,X(:,start:T),[Output1(:,start:T);Output2(:,start:T);Output3(:,start:T);Output4(:,start:T);Output5(:,start:T)]);
            end

            Approximator.ymaxin=net.inputs{1}.processSettings{1}.ymax;
            Approximator.yminin=net.inputs{1}.processSettings{1}.ymin;
            Approximator.xmaxin=net.inputs{1}.processSettings{1}.xmax;
            Approximator.xminin=net.inputs{1}.processSettings{1}.xmin;

            Approximator.ymaxout=net.outputs{2}.processSettings{1}.ymax;
            Approximator.yminout=net.outputs{2}.processSettings{1}.ymin;
            Approximator.xmaxout=net.outputs{2}.processSettings{1}.xmax;
            Approximator.xminout=net.outputs{2}.processSettings{1}.xmin;

            Approximator.IW = net.IW{1,1};
            Approximator.B1 = net.b{1};
            Approximator.LW = net.LW{2,1};
            Approximator.b2 = net.b{2};

            Poutput=SimNetwork(X(:,start:T),Approximator.ymaxin,Approximator.yminin,Approximator.xmaxin,Approximator.xminin,Approximator.ymaxout,Approximator.yminout,Approximator.xmaxout,Approximator.xminout,Approximator.IW,Approximator.B1,Approximator.LW,Approximator.b2);
        end

    end
end