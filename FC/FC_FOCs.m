function residuals = FC_FOCs(par,utility,sigmatilde,ctilde,phisigma,phic,grids,Pg,options)

beta  = par.beta;
delta = par.delta;
alpha = par.alpha;
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
                     
K=length(k_grd);
ONE=ones(K,1);

residuals=zeros(K,length(sigma_grd),length(g_grd));
for g_iter=1:length(g_grd)
    for sigma_iter=1:length(sigma_grd)

        g=g_grd(g_iter).*ONE;
        sigma=sigma_grd(sigma_iter).*ONE;
        
        sigmanext=sigmatilde(k_grd,sigma,g,phisigma);
        c=ctilde(k_grd,sigma,g,phic);
        knextF=@(l) z.*k_grd.^alpha.*l.^(1-alpha)+(1-delta)*k_grd-g-c;
        %Ulc=0 preference are separable
        % Equation 14 David Stockman
        lambdaF = @(l) du(c)+...
            -(sigmanext-sigma).*ddu(c).*knextF(l)+...
            +sigma.*(du(c)+ddu(c).*c);
        
        % Equation 15 David Stockman
        eq15 = @(l) dv(l)+lambdaF(l).*(1-alpha).*z.*k_grd.^(alpha).*l.^(-alpha)+sigma.*(dv(l)+ddv(l).*l);
        [l,fval]=fsolve(eq15,1.0.*ONE,options);
        
        knext=knextF(l);
        lambda=lambdaF(l);
        
        Eeuler1=0;
        Eeuler2=0;
        for g_iter_next=1:length(g_grd)
            g_next=g_grd(g_iter_next).*ONE;
            sigmanextnext=sigmatilde(knext,sigmanext,g_next,phisigma);
            cnext=ctilde(knext,sigmanext,g_next,phic);
            
            knextnextF=@(lnext) z.*knext.^alpha.*lnext.^(1-alpha)+(1-delta)*knext-g_next-cnext;
            lambdanextF = @(lnext) du(cnext)+...
                -(sigmanextnext-sigmanext).*ddu(cnext).*knextnextF(lnext)+...
                +sigma.*(du(cnext)+ddu(cnext).*cnext);
            
            % Equation 15 David Stockman
            eq15next = @(lnext) dv(lnext)+lambdanextF(lnext)*(1-alpha).*z_next.*knext.^(alpha).*lnext.^(-alpha)+...
                +sigmanext.*(dv(lnext)+ddv(lnext).*lnext);
            [lnext,fval]=fsolve(eq15next,1.0.*ONE,options);
            
            lambdanext=lambdanextF(lnext);
            knextnext=knextnextF(lnext);
            Eeuler1=Eeuler1+(du(cnext).*(cnext+knextnext)+dv(lnext).*lnext)*Pg(g_iter,g_iter_next);
            Eeuler2=Eeuler2+lambdanext.*(z_next.*alpha.*knext.^(alpha-1).*lnext.^(1-alpha)+1-delta)*Pg(g_iter,g_iter_next);
        end
        
        %Private sector FOC
        residuals(:,sigma_iter,g_iter,1)=du(c).*knext-beta*Eeuler1;
        %Ramsey inter-temporal
        residuals(:,sigma_iter,g_iter,2)=lambda-beta*Eeuler2+du(c).*(sigmanext-sigma);
    end
end