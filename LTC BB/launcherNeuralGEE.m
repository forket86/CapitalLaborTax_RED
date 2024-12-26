clear all;
clc;
close all;

addpath('../Utils/');

%% Choose a model to solve.
% 1. Optimal time consistent policy with limited state contigent labor and capital tax with symmetric cost=2
% 2. Optimal time consistent policy with limited state contigent labor and
% capital tax with symmetric cost=39
MODEL=1;
warning off
if MODEL==1
    load(['../FC BB/INIT/CurveGEE/Both/Curve2'])
elseif MODEL==2
    load(['../FC BB/INIT/CurveGEE/Both/Curve39'])
end
warning on

k_FC=k;
c_FC=c;
labor_FC=labor;
tauL_FC=tauL;
tauK_FC=tauK;

rng(1)
chi0_grd=chi0;
for chi0_iter=1:length(chi0_grd)
    chi0=chi0_grd(chi0_iter);
    psi0=chi0;
    display(['************** chi_0=' num2str(chi0) ' ********************'])

    %% Parameters
    par.chi0=chi0;
    par.psi0=psi0;

    chi=@(tauL,tauLbar) (chi0/2)*(tauL-tauLbar).^2;    utility.chi=chi;
    chi_tauL=@(tauL,tauLbar) chi0*(tauL-tauLbar);      utility.chi_tauL=chi_tauL;
    chi_tauLbar=@(tauL,tauLbar) -chi0*(tauL-tauLbar);  utility.chi_tauLbar=chi_tauLbar;
    psi=@(tauK,tauKbar) (psi0/2)*(tauK-tauKbar).^2;    utility.psi=psi;
    psi_tauK=@(tauK,tauKbar) psi0*(tauK-tauKbar);      utility.psi_tauK=psi_tauK;
    psi_tauKbar=@(tauK,tauKbar) -psi0*(tauK-tauKbar);  utility.psi_tauKbar=psi_tauKbar;


    %% Time simulation
    opt = optimoptions('fsolve','Display','none','Algorithm','Levenberg-Marquardt','MaxIterations',10000,'MaxFunctionEvaluation',10000,'FunctionTolerance',1e-9,'OptimalityTolerance',1e-14,'StepTolerance',1e-9);

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

    for iter=1:10000
        iter

        k(1)=k_FC(1);
        tauKbar(1)=tauK(1);
        tauLbar(1)=tauL(1);
        sigma(1)=sigma(1);

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
            tauK0=tauK(t);
            tauL0=tauL(t);
            sigmanext0=sigma(t+1);
            tauKbarnext0=tauKbar(t+1);
            tauLbarnext0=tauLbar(t+1);

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

            Output1(t)=(du(c(t)).*(c(t)+k(t+1))-dv(labor(t)).*labor(t));
            Output2(t)=lambda(t).*(Z.*alpha.*k(t).^(alpha-1).*labor(t).^(1-alpha)+1-delta);
            Output3(t)=(mu(t).*labor_k(t)-nu(t).*c_k(t));
            Output4(t)=tauL(t);
            Output5(t)=tauK(t);

        end
        toc

        FOC_residual=mean(fval,2)

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
            if MODEL==1
                save(['../FC BB/INIT/CurveGEE/Both/Curve2_temp'])
            elseif MODEL==2
                save(['../FC BB/INIT/CurveGEE/Both/Curve39_temp'])
            end
            break;
        end

        k_prev=k;
        c_prev=c;
        labor_prev=labor;
        tauK_prev=tauK;
        tauL_prev=tauL;

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