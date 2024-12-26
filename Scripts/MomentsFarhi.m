%% Generate Markov chain
infinityHorizon=10000;
T=infinityHorizon;
rng(1)

Simulate_Debt

labor=l;
tauLbar_temp = tauLbar(1:end-1);
tauKbar_temp = tauKbar(1:end-1);
U=u(c)-v(labor)-chi(tauL,tauLbar_temp)-psi(tauK,tauKbar_temp);
UtauL=u(c)-v(labor)-chi(tauL,tauLbar_temp);
UtauK=u(c)-v(labor)-psi(tauK,tauKbar_temp);
UHH=u(c)-v(labor);
CUTEND_WELFARE=500;
V_Fahri=zeros(T-CUTEND_WELFARE,1);
VtauL_Fahri=zeros(T-CUTEND_WELFARE,1);
VtauK_Fahri=zeros(T-CUTEND_WELFARE,1);
VHH_Fahri=zeros(T-CUTEND_WELFARE,1);
for t=1:T-CUTEND_WELFARE
    gridC=0:length(U(t:t+CUTEND_WELFARE))-1;
    V_Fahri(t)=beta.^gridC*U(t:t+CUTEND_WELFARE)';
    VtauL_Fahri(t)=beta.^gridC*UtauL(t:t+CUTEND_WELFARE)';
    VtauK_Fahri(t)=beta.^gridC*UtauK(t:t+CUTEND_WELFARE)';
    VHH_Fahri(t)=beta.^gridC*UHH(t:t+CUTEND_WELFARE)';
end
T0=50;Tend=T;

Fahrisamecost1=[
mean(tauK(T0:Tend));
mean(tauL(T0:Tend));
std(log(tauK(T0:Tend)+1));
std(log(tauL(T0:Tend)+1));
corr(log(tauK(T0+1:Tend)+1)',log(tauK(T0:Tend-1)+1)');
corr(log(tauL(T0+1:Tend)+1)',log(tauL(T0:Tend-1)+1)');
corr(log(tauK(T0:Tend)+1)',log(g_shocks(T0:Tend))');
corr(log(tauL(T0:Tend)+1)',log(g_shocks(T0:Tend))');];

Y=Z*k(1:Tend).^alpha.*labor(1:Tend).^(1-alpha);
I=k(2:Tend)-(1-delta)*k(1:Tend-1);
Fahrisamecost2=[
mean(Y);
mean(c(T0:Tend));
mean(labor(T0:Tend));
mean(I);
std(Y);
std(log(c(T0:Tend)));
std(log(labor(T0:Tend)));
std(log(I));
corr(log(Y(T0+1:Tend))',log(Y(T0:Tend-1))');
corr(log(c(T0+1:Tend))',log(c(T0:Tend-1))');
corr(log(labor(T0+1:Tend))',log(labor(T0:Tend-1))');
corr(log(I(T0+1:end))',log(I(T0:end-1))');
corr(log(Y(T0:Tend))',log(g_shocks(T0:Tend))');
corr(log(c(T0:Tend))',log(g_shocks(T0:Tend))');
corr(log(labor(T0:Tend))',log(g_shocks(T0:Tend))');
corr(log(I(T0:end))',log(g_shocks(T0:Tend-1))');
];

Fahrisamecost3=[
mean(tauKbar(T0:Tend));
mean(tauLbar(T0:Tend));
std(log(tauKbar(T0:Tend)+1));
std(log(tauLbar(T0:Tend)+1));
mean(tauL(T0:Tend)-tauLbar(T0:Tend));
mean(tauK(T0:Tend)-tauKbar(T0:Tend));
std(log(abs(tauL(T0:Tend)-tauLbar(T0:Tend))));
std(log(abs(tauK(T0:Tend)-tauKbar(T0:Tend))));];

