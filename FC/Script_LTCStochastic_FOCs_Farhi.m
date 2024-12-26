explore=0.1;
LOW_C=cFC*0.01;HIGH_C=cFC*1.5;
cFCunit=chebyTransformAle(LOW_C,HIGH_C,cFC,0);
LOW_K=kFC*(1-explore);HIGH_K=kFC*(1+explore);
LOW_TAUK=-0.15;HIGH_TAUK=0.15;
tauKFCunit=chebyTransformAle(LOW_TAUK,HIGH_TAUK,tauKFC,0);
LOW_TAUL=tauLFC*(1-explore); HIGH_TAUL=tauLFC*(1+explore);
tauLFCunit=chebyTransformAle(LOW_TAUL,HIGH_TAUL,tauLFC,0);
LOW_SIGMA=sigmaFC*(1-0.15);HIGH_SIGMA=sigmaFC*(1+0.1);
sigmaFCunit=chebyTransformAle(LOW_SIGMA,HIGH_SIGMA,sigmaFC,0);

% Grids structure
k_grd=linspace(LOW_K,HIGH_K,num_nodes)';           grids.k_grd=k_grd;
tauK_grd=linspace(LOW_TAUK,HIGH_TAUK,num_nodes)';  grids.tauK_grd=tauK_grd;
tauL_grd=linspace(LOW_TAUL,HIGH_TAUL,num_nodes)';  grids.tauL_grd=tauL_grd;
sigma_grd=linspace(LOW_SIGMA,HIGH_SIGMA,10);       grids.sigma_grd=sigma_grd;
grids.g_grd=g_grd; grids.z_grd=z_grd;
[grids.k_s,grids.tauK_s,grids.giter_s,grids.sigma_s]=ndgrid(k_grd,tauK_grd,1:length(g_grd),sigma_grd);

sigmatilde   =@(k,tauK,g,sigma,phisigma)    chebyAle(order,[k,tauK,g,sigma],phisigma,   [LOW_K,LOW_TAUK,g_grd(1),LOW_SIGMA],[HIGH_K,HIGH_TAUK,g_grd(end),HIGH_SIGMA],LOW_SIGMA,HIGH_SIGMA,cross);
tauKtilde    =@(k,tauK,g,sigma,phitauK)     chebyAle(order,[k,tauK,g,sigma],phitauK,    [LOW_K,LOW_TAUK,g_grd(1),LOW_SIGMA],[HIGH_K,HIGH_TAUK,g_grd(end),HIGH_SIGMA],LOW_TAUK,HIGH_TAUK,cross);
tauLtilde    =@(k,tauK,g,sigma,phitauL)     chebyAle(order,[k,tauK,g,sigma],phitauL,    [LOW_K,LOW_TAUK,g_grd(1),LOW_SIGMA],[HIGH_K,HIGH_TAUK,g_grd(end),HIGH_SIGMA],LOW_TAUL,HIGH_TAUL,cross);

numCrossTerms=14;
phisigma_LTC_Fahri=zeros(order+1,numCrossTerms);
phisigma_LTC_Fahri(1,:)=sigmaFCunit/numCrossTerms;

phitauK_LTC_Fahri=zeros(order+1,numCrossTerms);
phitauK_LTC_Fahri(1,:)=tauKFCunit/numCrossTerms;

phitauL_LTC_Fahri=zeros(order+1,numCrossTerms);
phitauL_LTC_Fahri(1,:)=tauLFCunit/numCrossTerms;

[r,c]=size(phisigma_LTC_Fahri);
display(['Num. of sigma   params: ' num2str(r*c)])
[r,c]=size(phitauK_LTC_Fahri);
display(['Num. of phitauK params: ' num2str(r*c)])
[r,c]=size(phitauL_LTC_Fahri);
display(['Num. of phitauL params: ' num2str(r*c)])

load('./INIT/Fahri/phisigma_LTC_Fahri')
load('./INIT/Fahri/phitauK_LTC_Fahri')
load('./INIT/Fahri/phitauL_LTC_Fahri')
clear phi
phi(:,:,1)=phisigma_LTC_Fahri;
phi(:,:,2)=phitauK_LTC_Fahri;
phi(:,:,3)=phitauL_LTC_Fahri;
obj=@(phi) LTC_FOCsVectorized_Farhi(par,utility,sigmatilde,tauKtilde,tauLtilde,phi(:,:,1),phi(:,:,2),phi(:,:,3),grids,Pg);

options=optimoptions(@fsolve,'Algorithm','Levenberg-Marquardt','MaxFunEvals',10000,'MaxIterations',10000,'Display','iter','StepTolerance', 1.0000e-010);
[phi,fval]=fsolve(obj,phi,options);

phisigma_LTC_Fahri=phi(:,:,1);
phitauK_LTC_Fahri=phi(:,:,2);
phitauL_LTC_Fahri=phi(:,:,3);

display(['MSE on sigma: ' num2str(sum(fval(:,1).^2)/(length(grids.k_s)))])
display(['MSE on phitauK: ' num2str(sum(fval(:,2).^2)/(length(grids.k_s)))])
display(['MSE on phitauL: ' num2str(sum(fval(:,3).^2)/(length(grids.k_s)))])