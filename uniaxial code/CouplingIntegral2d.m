function y = CouplingIntegral2d(w,et,ez,mut,muz,vm,vn,wm,wn,a,b,d,Mhxm,Mhxn,Mhym,Mhyn,lamx,lamy)

global eps0 mu0;

et=et*eps0;
ez=ez*eps0;
mut=mut*mu0;
muz=muz*mu0;

kym=(wm.*pi)/b;
kyn=(wn.*pi)/b;


Clamx= ((1-(-1).^vm.*exp(1j.*lamx.*a)).*(1-(-1).^vn.*exp(-1j.*lamx.*a)))./...
    ((lamx.^2-(vm.*pi./a).^2).*(lamx.^2-(vn.*pi./a).^2));
lamyterm=( (1-((-1).^(wm)).*exp(1j.*lamy.*b) ).* (1-((-1).^(wn)).*exp(-1j.*lamy.*b) ) )./...
    ((lamy.^2-kyn.^2).*(lamy.^2-kym.^2));
kt=sqrt(w^2*et*mut);
ktz=sqrt(w^2*et*muz);
kzt=sqrt(w^2*ez*mut);
lamr=sqrt(lamx.^2 + lamy.^2);
% lamzth=sqrt(kt.^2 - (mut/muz).*lamr.^2);
% lamzpsi=sqrt(kt.^2 - (et/ez).*lamr.^2);
lamzth=kt.*sqrt(1 - (lamr.^2./ktz.^2));
lamzpsi=kt.*sqrt(1 - (lamr.^2./kzt.^2));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% G_HH_xx
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

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
term1=( (Mhxm.*Mhxn.*vm.*vn.*lamy.^2)./ (a^2) ).*(ghh_xx);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% G_HH_xy
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% tez term of gfxy
ghh_xy_te_stable=(1j.*lamzth.*lamx.*lamy)./(lamr.^2.*w.*mut);
ghh_xy_te_unstable=1./sin(lamzth.*d);
ghh_xy_te_unstable(isnan(ghh_xy_te_unstable))=1j;

% tmz term of gfxy
ghh_xy_tm_stable=-(1j.*w.*et.*lamx.*lamy)./(lamr.^2.*lamzpsi);
ghh_xy_tm_unstable=1./sin(lamzpsi.*d);
ghh_xy_tm_unstable(isnan(ghh_xy_tm_unstable))=1j;

% total gfxy term
ghh_xy=(ghh_xy_te_stable.*ghh_xy_te_unstable) + (ghh_xy_tm_stable.*ghh_xy_tm_unstable);

% multiply by leading coefficients
term2=( (Mhxm.*Mhyn.*vm.*wn.*lamx.*lamy)./ (a*b) ).*(ghh_xy);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% G_HH_yx
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% tez term of gfyx
ghh_yx_te_stable=(1j.*lamzth.*lamx.*lamy)./(lamr.^2.*w.*mut);
ghh_yx_te_unstable=1./sin(lamzth.*d);
ghh_yx_te_unstable(isnan(ghh_yx_te_unstable))=1j;

% tmz term of gfyx
ghh_yx_tm_stable=-(1j.*w.*et.*lamx.*lamy)./(lamr.^2.*lamzpsi);
ghh_yx_tm_unstable=1./sin(lamzpsi.*d);
ghh_yx_tm_unstable(isnan(ghh_yx_tm_unstable))=1j;

% total gfyx term
ghh_yx=(ghh_yx_te_stable.*ghh_yx_te_unstable) + (ghh_yx_tm_stable.*ghh_yx_tm_unstable);

% multiply by leading coefficients
term3=( (Mhym.*Mhxn.*wm.*vn.*lamx.*lamy)./ (a*b) ).*(ghh_yx);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% G_HH_yy
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% tez term of gfyy
ghh_yy_te_stable=(1j.*lamzth.*lamy.^2)./(lamr.^2.*w.*mut);
ghh_yy_te_unstable=1./sin(lamzth.*d);
ghh_yy_te_unstable(isnan(ghh_yy_te_unstable))=1j;

% tmz term of gfyy
ghh_yy_tm_stable=(1j.*w.*et.*lamx.^2)./(lamr.^2.*lamzpsi);
ghh_yy_tm_unstable=1./sin(lamzpsi.*d);
ghh_yy_tm_unstable(isnan(ghh_yy_tm_unstable))=1j;

% total gfxx term
ghh_yy=(ghh_yy_te_stable.*ghh_yy_te_unstable) + (ghh_yy_tm_stable.*ghh_yy_tm_unstable);

% multiply by leading coefficients
term4=( (Mhym.*Mhyn.*wm.*wn.*lamx.^2)./ (b.^2) ).*(ghh_yy);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% Put it all together
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

y=Clamx.*(term1+term2+term3+term4).*lamyterm;
end
