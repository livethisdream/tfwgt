function y = SelfIntegral2d_1(w,et,ez,mut,muz,vm,vn,a,b,d,Mhxm,Mhxn,lamx,lamy)

global eps0 mu0;

et=et*eps0;
ez=ez*eps0;
mut=mut*mu0;
muz=muz*mu0;

Clamx= ((1-(-1).^vm.*exp(1j.*lamx.*a)).*(1-(-1).^vn.*exp(-1j.*lamx.*a)))./...
    ((lamx.^2-(vm.*pi./a).^2).*(lamx.^2-(vn.*pi./a).^2));
kt=sqrt(w^2*et*mut);
ktz=sqrt(w^2*et*muz);
kzt=sqrt(w^2*ez*mut);
lamr=sqrt(lamx.^2 + lamy.^2);

% % tez term
lamzth=sqrt(kt.^2 - (mut/muz).*lamr.^2);
ghh_xx_te_stable=(1j.*lamzth.*lamx.^2)./(lamr.^2.*w.*mut);
ghh_xx_te_unstable=cos(lamzth.*d)./sin(lamzth.*d);
ghh_xx_te_unstable(isnan(ghh_xx_te_unstable))=1j;

% tmz term of gf
lamzpsi=sqrt(kt.^2 - (et/ez).*lamr.^2);
ghh_xx_tm_stable=(1j.*w.*et.*lamy.^2)./(lamr.^2.*lamzpsi);
ghh_xx_tm_unstable=cos(lamzpsi.*d)./sin(lamzpsi.*d);
ghh_xx_tm_unstable(isnan(ghh_xx_tm_unstable))=1j;

ghh_xx=(ghh_xx_te_stable.*ghh_xx_te_unstable) + (ghh_xx_tm_stable.*ghh_xx_tm_unstable);

lamyterm=( (1-exp(1j.*lamy.*b) ).* (1-exp(-1j.*lamy.*b) ) )./(lamy.^2);


y=( (Mhxm.*Mhxn.*vm.*vn)./ (a^2) ).*(Clamx.*ghh_xx.*lamyterm);

end
