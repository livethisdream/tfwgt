function y = CouplingIntegral_1(w,et,ez,mut,muz,vm,vn,a,b,d,Mhxm,Mhxn,lamx)

global eps0 mu0;

Clamx= ((1-(-1).^vm.*exp(1i.*lamx.*a)).*(1-(-1).^vn.*exp(-1i.*lamx.*a)))./...
    ((lamx.^2-(vm.*pi./a).^2).*(lamx.^2-(vn.*pi./a).^2));

% % tez term
% the sum term
lmax=100;
kt=sqrt(w^2*eps0*et*mu0*mut);
ktz=sqrt(w^2*eps0*et*mu0*muz);
kzt=sqrt(w^2*eps0*ez*mu0*mut);
[l,lamxl] = ndgrid(0:lmax,lamx);
% if imag(muz/mut) < 0
%     lamylth = sqrt( (mu0*muz)/(mu0*mut) )*sqrt( kt^2-(pi*l/d).^2-( ( (mu0*mut)/(mu0*muz) )*lamxl.^2));    
% else
%     lamylth = sqrt( ( (mu0*muz)/(mu0*mut) )*(kt^2-(pi*l/d).^2)-lamxl.^2);
% end
lamylth=sqrt(ktz^2)*sqrt(1-((pi^2*l.^2)/(d^2*kt^2))-((lamxl.^2)/(ktz^2)));
sumterm1= -((4*pi*mu0*muz.*lamx.^2.*Clamx)./(w*mu0^2*mut^2*d)).*...
    sum( ( ( (-1).^l).*( (pi*l/d).^2 ).*(1-exp(-1i*lamylth*b) ) )./...
    ( lamylth.^3.*(lamylth.^2 + lamxl.^2) ),1);

% the other term
% lamztha=sqrt( w^2*eps0*et*mu0*mut - ((mu0*mut)/(mu0*muz))*lamx.^2 );
lamztha=kt*sqrt( 1 - (lamx/ktz).^2 );
term2stable=((1j.*2*pi*b.*lamztha.*Clamx)./(w*mut*mu0));
term2unstable=(1)./(sin(lamztha*d));
% if any(isnan(term2unstable))==1
%     display 'found nan'
% end
term2unstable(isnan(term2unstable))=1j;

% the whole tez term
omega12tez=(term2stable.*term2unstable + sumterm1);

% % tmz term
% if imag(ez/et) < 0
%     lamylpsi = sqrt( (eps0*ez)/(eps0*et) )*sqrt( kt^2-(pi*l/d).^2-( ( (eps0*et)/(eps0*ez) )*lamxl.^2));
% else
%     lamylpsi = sqrt( ( (eps0*ez)/(eps0*et) )*(kt^2-(pi*l/d).^2)-lamxl.^2);
% end
lamylpsi=sqrt(kzt^2)*sqrt(1-((pi^2*l.^2)/(d^2*kt^2))-((lamxl.^2)/(kzt^2)));
omega12tmz=-((4*pi.*w.*ez.*eps0.*Clamx)./(d)).*sum((1-exp(-1j.*lamylpsi.*b))./...
    (lamylpsi.*(lamylpsi.^2 + lamxl.^2).*( ( (-1).^l) +1.*(l==0))),1);

y=( (Mhxm.*Mhxn.*vm.*vn)./ (a^2) ).*(omega12tez+omega12tmz);

end

