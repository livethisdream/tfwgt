function y = SelfIntegral(w,et,ez,mut,muz,vm,vn,a,b,d,lamx)

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
% if imag(muz/mut) < 0
%     lamylth = sqrt( ( (mu0*muz)/(mu0*mut) ) )*sqrt( kt^2-(pi*l/d).^2-( ( (mu0*mut)/(mu0*muz) )*lamxl.^2));    
% else
%     lamylth = sqrt( ( (mu0*muz)/(mu0*mut) )*(kt^2-(pi*l/d).^2)-lamxl.^2);
% end
lamylth=sqrt(ktz^2)*sqrt(1-((pi^2*l.^2)/(d^2*kt^2))-((lamxl.^2)/(ktz^2)));
sumterm1= ((Clamx.*lamx.^2*1j*4*pi*mu0*muz)./(w*mu0^2*mut^2*d)).*...
    sum( ( ( (pi*l/d).^2 ).*(1-exp(-1j*lamylth*b) ) )./( lamylth.^3.*(lamylth.^2 + lamxl.^2) ),1);

% the other term
% lamztha2=sqrt(kt^2 - ((mu0*mut)/(mu0*muz))*lamx.^2 );
lamztha=kt*sqrt( 1 - ( lamx/ktz ).^2 );
term2stable=((2*pi*b.*lamztha.*Clamx)./(w*mut*mu0));
term2unstable=((cos(lamztha*d))./(sin(lamztha*d)));
term2unstable(isnan(term2unstable))=1j;

% the whole tez term
omega11tez=(term2stable.*term2unstable + sumterm1);

% % tmz term
% lamylpsi= sqrt( (ez./et).*(kt^2-(pi*l./d).^2)-lamxl.^2);
% if imag(ez/et) < 0
%     lamylpsi = sqrt( ( (eps0*ez)/(eps0*et) ) )*sqrt( kt^2-(pi*l/d).^2-( ( (eps0*et)/(eps0*ez) )*lamxl.^2));
% else
%     lamylpsi = sqrt( ( (eps0*ez)/(eps0*et) )*(kt^2-(pi*l/d).^2)-lamxl.^2);
% end
lamylpsi=sqrt(kzt^2)*sqrt(1-((pi^2*l.^2)/(d^2*kt^2))-((lamxl.^2)/(kzt^2)));
omega11tmz=((1j.*w.*ez.*eps0.*4.*pi.*Clamx)./(d)).*sum( (1-exp(-1j.*lamylpsi.*b))./...
    (lamylpsi.*(lamylpsi.^2 + lamxl.^2).*(1+1.*(l==0))),1);

y=omega11tez+omega11tmz;

end
