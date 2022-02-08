function y = CouplingIntegral_3(w,et,ez,mut,muz,wn,vm,vn,a,b,d,Mhxm,Mhxn,Mhyn,lamx)

global eps0 mu0;

% waveguide parameters
kyn=wn*pi/b;
kxm=vm*pi/a;
kxn=vn*pi/a;

Clamx= ((1-(-1).^vm.*exp(1j.*lamx.*a)).*(1-(-1).^vn.*exp(-1j.*lamx.*a)))./...
    ((lamx.^2-(vm.*pi./a).^2).*(lamx.^2-(vn.*pi./a).^2));

% % A_lamy term
% the tez sum term
lmax=100;
kt=sqrt(w^2*eps0*et*mu0*mut);
ktz=sqrt(w^2*eps0*et*mu0*muz);
kzt=sqrt(w^2*eps0*ez*mu0*mut);
[l,lamxl] = ndgrid(0:lmax,lamx);
lamylth=sqrt(ktz^2)*sqrt(1-((pi^2*l.^2)/(d^2*kt^2))-((lamxl.^2)/(ktz^2)));
A11sumterm1= ((mu0*muz.*lamx.^2.*Clamx)./(w*mu0^2*mut^2)).*...
    sum( ( ( (-1).^l ).*( (pi*l/d).^2 ).*(1-exp(-1j*lamylth*b) ) )./...
        ( lamylth.*(lamylth.^2 + lamxl.^2).*(lamylth.^2 - kyn^2) ),1);


% % tmz sum term
lamylpsi=sqrt(kzt^2)*sqrt(1-((pi^2*l.^2)/(d^2*kt^2))-((lamxl.^2)/(kzt^2)));
A11sumterm2=(w.*ez.*eps0.*Clamx).*sum( ( lamylpsi.*(1-exp(-1j.*lamylpsi.*b) ) )./...
    ( (lamylpsi.^2 + lamxl.^2).*(lamylpsi.^2 - kyn^2).*( ( (-1).^l ) + 1.*(l==0) ) ),1);

% put it all together
A11lamy= -( (4*pi*Mhxm*Mhxn*vm*vn  )./ (d* a^2 ) ).*(A11sumterm1+A11sumterm2);

% % B_lamy term

% the tez sum term
B11sumterm1= -((Clamx.*mu0*muz)./(w*mu0^2*mut^2)).*...
    sum( ( ( (-1).^l ).*( (pi*l/d).^2 ).*(1-exp(-1j*lamylth*b) ) )./...
    ( lamylth.*(lamylth.^2 + lamxl.^2).*(lamylth.^2 - kyn^2) ),1);

% % tmz sum term
B11sumterm2=(w.*ez.*eps0.*Clamx).*sum( (1-exp(-1j.*lamylpsi.*b) )./...
    ( lamylpsi.*(lamylpsi.^2 + lamxl.^2).*(lamylpsi.^2 - kyn^2).*( ( (-1).^l ) + 1.*(l==0) ) ),1);

% put it all together
B11lamy= ( (4*pi*Mhxm*Mhyn*vm*wn.*lamx.^2)./ (a*b*d) ).*(B11sumterm1+B11sumterm2);

y=A11lamy+B11lamy;

end
