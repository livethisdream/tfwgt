function y = CouplingIntegral_5(w,et,ez,mut,muz,wm,wn,vm,vn,a,b,d,Mhxm,Mhxn,Mhym,Mhyn,lamx)

global eps0 mu0;

% some special terms we'll need throughout
kym=wm*pi/b;
kyn=wn*pi/b;
kxm=vm*pi/a;
kxn=vn*pi/a;
kt=sqrt(w^2*eps0*et*mu0*mut);
ktz=sqrt(w^2*eps0*et*mu0*muz);
kzt=sqrt(w^2*eps0*ez*mu0*mut);

Clamx= ((1-(-1).^vm.*exp(1j.*lamx.*a)).*(1-(-1).^vn.*exp(-1j.*lamx.*a)))./...
    ((lamx.^2-(vm.*pi./a).^2).*(lamx.^2-(vn.*pi./a).^2));
lmax=100;
[l,lamxl] = ndgrid(0:lmax,lamx);
lamylth=sqrt(ktz^2)*sqrt(1-((pi^2*l.^2)/(d^2*kt^2))-((lamxl.^2)/(ktz^2)));
lamylpsi=sqrt(kzt^2)*sqrt(1-((pi^2*l.^2)/(d^2*kt^2))-((lamxl.^2)/(kzt^2)));


% % % A11lamy term
% % tez term
% the sum term
Aterm1= ((4*pi*mu0*muz.*lamx.^2.*Clamx)./(w*mu0^2*mut^2*d)).*...
    sum(...
        ( ( (-1).^(l) ).*( (pi*l/d).^2 ).*(lamylth).*(1-exp(-1j*lamylth*b) ) )./...
        ( (lamxl.^2 + lamylth.^2).*(lamylth.^2 - kym.^2).*(lamylth.^2 - kyn.^2) )...
        ,1);

% % tmz term
% the sum term 
Aterm2=( (4.*pi.*w.*ez.*eps0.*Clamx)./(d) ).*...
    sum(...
    ( (lamylpsi.^3).*( 1-exp(-1j.*lamylpsi.*b) ) )./...
    ( (lamylpsi.^2 + lamxl.^2).*(lamylpsi.^2 - kym.^2).*(lamylpsi.^2 - kyn.^2).*( ( (-1).^l) + 1.*(l==0) ) )...
        ,1);


% put it all together
Alamy=-( (Mhxm.*Mhxn.*vm.*vn)./ (a^2) ).*(Aterm1+Aterm2);

% % % B11lamy term
% % tez term
% the sum term
Bterm1= -( (4*pi*mu0*muz.*lamx.^2.*Clamx)./(w*mu0^2*mut^2*d) ).*...
    sum( ...
	   (  ( (-1).^(l) ).*( (pi*l/d).^2 ).*(lamylth).*( 1-exp(-1j*lamylth*b) )  )./...
       ( (lamxl.^2 + lamylth.^2).*(lamylth.^2 - kym.^2).*(lamylth.^2 - kyn.^2) )...
      ,1);


% % tmz term
% the sum term 
Bterm2= (  (4.*pi.*w.*ez.*eps0.*lamx.^2.*Clamx)./(d) ).*...
		sum( ...
			(  (lamylpsi).*( 1-exp(-1j.*lamylpsi.*b) )  )./...
    		(  (lamylpsi.^2 + lamxl.^2).*(lamylpsi.^2 - kym.^2).*(lamylpsi.^2 - kyn.^2).*( (-1).^(l) + 1.*(l==0) )  )...
            ,1);


% put it all together
Blamy=( (Mhxm.*Mhyn.*wn.*vm)./ (a*b) ).*(Bterm1+Bterm2);

% % % C11lamy term
% % tez term
% the sum term
Cterm1= -( (4*pi*mu0*muz.*lamx.^2.*Clamx)./(w*mu0^2*mut^2*d) ).*...
    sum( ...
	(  ( (-1).^(l) ).*( (pi*l/d).^2 ).*(lamylth).*( 1-exp(-1j*lamylth*b) )  )./...
      ( (lamxl.^2 + lamylth.^2).*(lamylth.^2 - kym.^2).*(lamylth.^2 - kyn.^2) )...
       ,1);

% % tmz term
% the sum term 
Cterm2= (  (4.*pi.*w.*ez.*eps0.*lamx.^2.*Clamx)./(d) ).*...
		sum( ...
			(  (lamylpsi).*( 1-exp(-1j.*lamylpsi.*b) )  )./...
    		(  (lamylpsi.^2 + lamxl.^2).*(lamylpsi.^2 - kym.^2).*(lamylpsi.^2 - kyn.^2).*( (-1).^(l) + 1.*(l==0) )  ),1 ...
		    );


% put it all together
Clamy=( (Mhym.*Mhxn.*wm.*vn)./ (a*b) ).*(Cterm1+Cterm2);


% % % D11lamy term
% % tez term
% the sum term
Dterm1= ( (4*pi*mu0*muz.*lamx.^2.*Clamx)./(w*mu0^2*mut^2*d) ).*...
    sum( ...
	(  ( (-1).^(l) ).*( (pi*l/d).^2 ).*(lamylth).*( 1-exp(-1j*lamylth*b) )  )./...
      ( (lamxl.^2 + lamylth.^2).*(lamylth.^2 - kym.^2).*(lamylth.^2 - kyn.^2) )...
      ,1);

% % tmz term
% the sum term 
Dterm2=  (  (4.*pi.*w.*ez.*eps0.*lamx.^4.*Clamx)./(d) ).*...
		sum( ...
			(  ( 1-exp(-1j.*lamylpsi.*b) )  )./...
    		(  (lamylpsi).*(lamylpsi.^2 + lamxl.^2).*(lamylpsi.^2 - kym.^2).*(lamylpsi.^2 - kyn.^2).*( (-1).^(l) + 1.*(l==0) )  )...
            ,1);

% put it all together
Dlamy=-( (Mhym.*Mhyn.*wm.*wn)./ (b^2) ).*(Dterm1+Dterm2);

y=Alamy+Blamy+Clamy+Dlamy;

end








