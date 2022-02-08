function y = CouplingIntegral2d_2(wval,et,ez,mut,muz,wm,vm,vn,a,b,d,Mhxm,Mhxn,Mhym,lamx,lamy)

global eps0 mu0;

et=et*eps0;
ez=ez*eps0;
mut=mut*mu0;
muz=muz*mu0;

kym=(wm.*pi)/b;

Clamx= ((1-(-1).^vm.*exp(1j.*lamx.*a)).*(1-(-1).^vn.*exp(-1j.*lamx.*a)))./...
    ((lamx.^2-(vm.*pi./a).^2).*(lamx.^2-(vn.*pi./a).^2));
lamyterm=( (1-exp(1j.*lamy.*b) ).* (1-exp(-1j.*lamy.*b) ) )./(lamy.^2.*(lamy.^2-kym.^2));
kt=sqrt(w^2*et*mut);
ktz=sqrt(w^2*et*muz);
kzt=sqrt(w^2*ez*mut);
lamr=sqrt(lamx.^2 + lamy.^2);
lamzth=sqrt(kt.^2 - (mut/muz).*lamr.^2);
lamzpsi=sqrt(kt.^2 - (et/ez).*lamr.^2);

% tez term of gfxx
ghh_xx_te_stable=(1j.*lamzth.*lamx.^2)./(lamr.^2.*w.*mut);
ghh_xx_te_unstable=1./sin(lamzth.*d);
ghh_xx_te_unstable(isnan(ghh_xx_te_unstable))=1j;

% tmz term of gfxx
ghh_xx_tm_stable=(1j.*w.*et.*lamy.^2)./(lamr.^2.*lamzpsi);
ghh_xx_tm_unstable=1./sin(lamzpsi.*d);
ghh_xx_tm_unstable(isnan(ghh_xx_tm_unstable))=1j;

% total gfxx term
ghh_xx=(ghh_xx_te_stable.*ghh_xx_te_unstable) + (ghh_xx_tm_stable.*ghh_xx_tm_unstable);

% multiply by leading coefficients
y1=( (Mhxm.*Mhxn.*vm.*vn.*lamy.^2)./ (a^2) ).*(ghh_xx);


%%%% left off here...
% tez term of gfyx
ghh_yx_te_stable=(1j.*lamzth.*lamx.^2)./(lamr.^2.*w.*mut);
ghh_yx_te_unstable=cos(lamzth.*d)./sin(lamzth.*d);
ghh_yx_te_unstable(isnan(ghh_xx_te_unstable))=1j;

% tmz term of gfyx
ghh_yx_tm_stable=(1j.*w.*et.*lamy.^2)./(lamr.^2.*lamzpsi);
ghh_yx_tm_unstable=cos(lamzpsi.*d)./sin(lamzpsi.*d);
ghh_yx_tm_unstable(isnan(ghh_xx_tm_unstable))=1j;

% multiply by leading coefficients
y2=( (Mhym.*Mhxn.*wm.*vn.*lamx.*lamy)./ (a*b) ).*(ghh_yx);

y=Clamx.*(y1+y2).*lamyterm;
end
