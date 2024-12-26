LOW_C=cFC*0.1;HIGH_C=cFC*1.9;
cFCunit=chebyTransformAle(LOW_C,HIGH_C,cFC,0);
LOW_K=kFC*(1-0.1);HIGH_K=kFC*(1+0.1);
LOW_TAUK=-0.4;HIGH_TAUK=0.4;
LOW_TAUL=tauLFC*(1-0.4); HIGH_TAUL=tauLFC*(1+0.4);
LOW_SIGMA=sigmaFC*(1-0.2);HIGH_SIGMA=sigmaFC*(1+0.2);
sigmaFCunit=chebyTransformAle(LOW_SIGMA,HIGH_SIGMA,sigmaFC,0);

% Grids structure
k_grd=linspace(LOW_K,HIGH_K,num_nodes)';           grids.k_grd=k_grd;
tauK_grd=linspace(LOW_TAUK,HIGH_TAUK,num_nodes)';  grids.tauK_grd=tauK_grd;
tauL_grd=linspace(LOW_TAUL,HIGH_TAUL,num_nodes)';  grids.tauL_grd=tauL_grd;
sigma_grd=linspace(LOW_SIGMA,HIGH_SIGMA,10);       grids.sigma_grd=sigma_grd;
grids.g_grd=g_grd; grids.z_grd=z_grd;
[grids.k_s,grids.sigma_s,grids.giter_s]=ndgrid(k_grd,sigma_grd,1:length(g_grd));

sigmatilde=@(k,sigma,g,phisigma) chebyAle(order,[k,sigma,g],phisigma,[LOW_K,LOW_SIGMA,g_grd(1)],[HIGH_K HIGH_SIGMA g_grd(end)],LOW_SIGMA,HIGH_SIGMA,cross);
ctilde    =@(k,sigma,g,phic)     chebyAle(order,[k,sigma,g],phic,    [LOW_K,LOW_SIGMA,g_grd(1)],[HIGH_K HIGH_SIGMA g_grd(end)],LOW_C,HIGH_C,cross);

phisigma_FC=zeros(order+1,3+cross*2);
phisigma_FC(1,1)=sigmaFCunit/5;phisigma_FC(1,2)=sigmaFCunit/5;phisigma_FC(1,3)=sigmaFCunit/5;phisigma_FC(1,4)=sigmaFCunit/5;phisigma_FC(1,5)=sigmaFCunit/5;

phic_FC=zeros(order+1,3+cross*2);
phic_FC(1,1)=cFCunit/5;phic_FC(1,2)=cFCunit/5;phic_FC(1,3)=cFCunit/5;phic_FC(1,4)=cFCunit/5;phic_FC(1,5)=cFCunit/5;

[r,c]=size(phisigma_FC);
display(['Num. of sigma params: ' num2str(r*c)])
[r,c]=size(phic_FC);
display(['Num. of c params: ' num2str(r*c)])

load('./INIT/FC/phic_FC')
load('./INIT/FC/phisigma_FC')

phi(:,:,1)=phisigma_FC;
phi(:,:,2)=phic_FC;
obj=@(phi) FC_FOCs(par,utility,sigmatilde,ctilde,phi(:,:,1),phi(:,:,2),grids,Pg,options);
[phi,fval]=fsolve(obj,phi,options);
%fval=obj(phi);
phisigma_FC=phi(:,:,1);
phic_FC=phi(:,:,2);

display(['MSE on sigma: ' num2str(sum(fval(:,1).^2)/(num_nodes*length(g_grd)))])
display(['MSE on c: ' num2str(sum(fval(:,2).^2)/(num_nodes*length(g_grd)))])