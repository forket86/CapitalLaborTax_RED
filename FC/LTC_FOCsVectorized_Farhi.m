function residuals = LTC_FOCsVectorized_Fahri(par,utility,sigmatilde,tauKtilde,tauLtilde,phisigma,phitauK,phitauL,grids,Pg)

beta  = par.beta;
delta = par.delta;
alpha = par.alpha;
gamma = par.gamma;
eta   = par.eta;
D     = par.D;

du    = utility.du;
ddu   = utility.ddu;
dv    = @(l) -utility.dv(l);
ddv   = @(l) -utility.ddv(l);

k_grd=grids.k_grd;
sigma_grd=grids.sigma_grd;
z_grd=grids.z_grd;
z=z_grd;
z_next=z;
g_grd=grids.g_grd;
tauK_grd=grids.tauK_grd;

K=length(k_grd);

k=grids.k_s(:);
ONES=ones(length(k),1);
tauK=grids.tauK_s(:);
g_iter=grids.giter_s(:);
g=g_grd(g_iter);
sigma=grids.sigma_s(:);

tauL=tauLtilde(k,tauK,g,sigma,phitauL);
sigmanext=sigmatilde(k,tauK,g,sigma,phisigma);

[l,c,knext,ltauL,ctauL] = getLaborCKnext(k,tauK,tauL,z,g,alpha,delta,gamma,eta,D);

lambda=(du(c)-(sigmanext-sigma).*ddu(c).*knext+sigma.*(du(c)+ddu(c).*c)+((dv(l)+sigma.*(dv(l)+ddv(l).*l))).*ltauL./ctauL)./(1-(1-alpha).*z.*k.^(alpha).*l.^(-alpha).*ltauL./ctauL);
%nu = du(c)-(sigmanext-sigma).*ddu(c).*knext+sigma.*(du(c)+ddu(c).*c)-lambda;
%mu=(nu.*ctauL-chi(tauL,tauLbar))./ltauL;

%% BEGIN Check partial derivatives with numerical ones
% [~,~,~,labor_tauL,c_tauL,labor_tauK,c_tauK,labor_k,c_k] = getLaborCKnext(k,tauK,tauL,z,g,alpha,delta,gamma,eta,D);
% Delta=1e-4;
% lU_k    = getLaborCKnext(k+Delta,tauK,tauL,z,g,alpha,delta,gamma,eta,D);  lD_k      = getLaborCKnext(k-Delta,tauK,tauL,z,g,alpha,delta,gamma,eta,D);
% lU_tauK = getLaborCKnext(k,tauK+Delta,tauL,z,g,alpha,delta,gamma,eta,D);  lD_tauK   = getLaborCKnext(k,tauK-Delta,tauL,z,g,alpha,delta,gamma,eta,D);
% lU_tauL = getLaborCKnext(k,tauK,tauL+Delta,z,g,alpha,delta,gamma,eta,D);  lD_tauL   = getLaborCKnext(k,tauK,tauL-Delta,z,g,alpha,delta,gamma,eta,D);
% labor_k2=(lU_k-lD_k)/(2*Delta);          err1=max(abs(labor_k2-labor_k));
% labor_tauK2=(lU_tauK-lD_tauK)/(2*Delta); err2=max(abs(labor_tauK2-labor_tauK));
% labor_tauL2=(lU_tauL-lD_tauL)/(2*Delta); err3=max(abs(labor_tauL2-labor_tauL));
% 
% 
% [~,cU_k]  = getLaborCKnext(k+Delta,tauK,tauL,z,g,alpha,delta,gamma,eta,D);    [~,cD_k]  = getLaborCKnext(k-Delta,tauK,tauL,z,g,alpha,delta,gamma,eta,D);
% [~,cU_tauK]  = getLaborCKnext(k,tauK+Delta,tauL,z,g,alpha,delta,gamma,eta,D); [~,cD_tauK]  = getLaborCKnext(k,tauK-Delta,tauL,z,g,alpha,delta,gamma,eta,D);
% [~,cU_tauL]  = getLaborCKnext(k,tauK,tauL+Delta,z,g,alpha,delta,gamma,eta,D); [~,cD_tauL]  = getLaborCKnext(k,tauK,tauL-Delta,z,g,alpha,delta,gamma,eta,D);
% c_k2=(cU_k-cD_k)/(2*Delta);              err4=max(abs(c_k2-c_k));
% c_tauK2=(cU_tauK-cD_tauK)/(2*Delta);     err5=max(abs(c_tauK2-c_tauK));
% c_tauL2=(cU_tauL-cD_tauL)/(2*Delta);     err6=max(abs(c_tauL2-c_tauL));
% 
% err=[err1,err2,err3,err4,err5,err6];

%% END Check partial derivatives with numerical ones

%% BEGIN Countercheck with Stockman
% % Eq. 14 David Stockman (with nu)
% check1 = du(c)-(sigmanext-sigma).*ddu(c).*knext+sigma.*(du(c)+ddu(c).*c)-nu-lambda;
% 
% % Eq. 15 David Stockman (with mu)
% check2=-(dv(l)+lambda.*(1-alpha).*z.*k.^(alpha).*l.^(-alpha)+sigma.*(dv(l)+ddv(l).*l))-mu;
% check3=nu.*ctauL-chi(tauL,tauLbar)-mu.*ltauL;
%% END Countercheck with Stockman

tauKnext=tauKtilde(k,tauK,g,sigma,phitauK);
Eeuler1=0;
Eeuler2A=0;
Eeuler2B=0;
Euler3 = 0;
for g_iter_next=1:length(g_grd)
    g_next=g_grd(g_iter_next).*ONES;
    tauLnext=tauLtilde(knext,tauKnext,g_next,sigmanext,phitauL);
    sigmanextnext=sigmatilde(knext,tauKnext,g_next,sigmanext,phisigma);
    
    [lnext,cnext,knextnext,lnext_tauL,cnext_tauL,lnext_tauK,cnext_tauK,lnext_k,cnext_k] = getLaborCKnext(knext,tauKnext,tauLnext,z,g_next,alpha,delta,gamma,eta,D);
    
    lambdanext=(du(cnext)-(sigmanextnext-sigmanext).*ddu(cnext).*knextnext+sigmanext.*(du(cnext)+ddu(cnext).*cnext)+((dv(lnext)+sigmanext.*(dv(lnext)+ddv(lnext).*lnext))).*lnext_tauL./cnext_tauL)./(1-(1-alpha).*z.*knext.^(alpha).*lnext.^(-alpha).*lnext_tauL./cnext_tauL);
    nunext = du(cnext)-(sigmanextnext-sigmanext).*ddu(cnext).*knextnext+sigmanext.*(du(cnext)+ddu(cnext).*cnext)-lambdanext;
    % Eq.17
    munext=(nunext.*cnext_tauL)./lnext_tauL;
    
    %% BEGIN Countercheck with Stockman
%     % Eq. 14 David Stockman (with nu)
%     check1next = du(cnext)-(sigmanextnext-sigmanext).*ddu(cnext).*knextnext+sigmanext.*(du(cnext)+ddu(cnext).*cnext)-nunext-lambdanext;
%     
%     % Eq. 15 David Stockman (with mu)
%     check2next=-(dv(lnext)+lambdanext.*(1-alpha).*z.*knext.^(alpha).*lnext.^(-alpha)+sigmanext.*(dv(lnext)+ddv(lnext).*lnext))-munext;
%     check3next=nunext.*cnext_tauL-chi(tauLnext,tauLbarnext)-munext.*lnext_tauL;
    %% END Countercheck with Stockman
    
    Eeuler1=Eeuler1+(du(cnext).*(cnext+knextnext)+dv(lnext).*lnext).*Pg(g_iter,g_iter_next);
    % Expectations in Eq.15
    Eeuler2A=Eeuler2A+lambdanext.*(z_next.*alpha.*knext.^(alpha-1).*lnext.^(1-alpha)+1-delta).*Pg(g_iter,g_iter_next);
    Eeuler2B=Eeuler2B+(munext.*lnext_k-nunext.*cnext_k).*Pg(g_iter,g_iter_next);
    % Expectation in Eq.16
    Euler3=Euler3+(munext.*lnext_tauK-nunext.*cnext_tauK).*Pg(g_iter,g_iter_next);
end

%Implementability
residuals(:,:,:,:,1)=reshape(du(c).*knext-beta*Eeuler1,K,length(tauK_grd),length(g_grd),length(sigma_grd));
%Eq. 15
residuals(:,:,:,:,2)=reshape(lambda-beta*Eeuler2A+du(c).*(sigmanext-sigma)+beta*Eeuler2B,K,length(tauK_grd),length(g_grd),length(sigma_grd));
%Eq. 16
residuals(:,:,:,:,3)=reshape(Euler3,K,length(tauK_grd),length(g_grd),length(sigma_grd));