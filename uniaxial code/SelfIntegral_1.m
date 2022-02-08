function y = SelfIntegral_1(w,et,ez,mut,muz,vm,vn,a,b,d,Mhxm,Mhxn,lamx)

global eps0 mu0;
Clamx= ((1-(-1).^vm.*exp(1j.*lamx.*a)).*(1-(-1).^vn.*exp(-1j.*lamx.*a)))./...
    ((lamx.^2-(vm.*pi./a).^2).*(lamx.^2-(vn.*pi./a).^2));

% % tez term
% the sum term
lmax=100;
kt=sqrt(w^2*eps0*et*mu0*mut);
ktz=sqrt(w^2*eps0*et*mu0*muz);
kzt=sqrt(w^2*eps0*ez*mu0*mut);
[l,lamxl] = ndgrid(0:lmax,lamx);
lamylth=sqrt(ktz^2)*sqrt(1-((pi^2*l.^2)/(d^2*kt^2))-((lamxl.^2)/(ktz^2)));
sumterm1= -((4*pi*mu0*muz.*lamx.^2.*Clamx)./(w*mu0^2*mut^2*d)).*...
    sum( ( ( (pi*l/d).^2 ).*(1-exp(-1j*lamylth*b) ) )./( lamylth.^3.*(lamylth.^2 + lamxl.^2) ),1);

% the other term
lamztha=kt*sqrt( 1 - ( lamx/ktz ).^2 );
term2stable=((1j*2*pi*b.*lamztha.*Clamx)./(w*mut*mu0));
term2unstable=((cos(lamztha*d))./(sin(lamztha*d)));
% if any(isnan(term2unstable))==1
%     display 'found nan'
% end
term2unstable(isnan(term2unstable))=1j;

% the whole tez term
omega11tez=(term2stable.*term2unstable + sumterm1);

% % tmz term
lamylpsi=sqrt(kzt^2)*sqrt(1-((pi^2*l.^2)/(d^2*kt^2))-((lamxl.^2)/(kzt^2)));
omega11tmz=-((4.*pi.*w.*ez.*eps0.*Clamx)./(d)).*sum( (1-exp(-1j.*lamylpsi.*b))./...
    (lamylpsi.*(lamylpsi.^2 + lamxl.^2).*(1+1.*(l==0))),1);

y=( (Mhxm.*Mhxn.*vm.*vn)./ (a^2) ).*(omega11tez+omega11tmz);

end
