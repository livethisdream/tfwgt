function y = SelfIntegral_4(w,et,ez,mut,muz,wm,vm,vn,a,b,d,Mhxm,Mhxn,Mhym,Mhyn,lamx)

global eps0 mu0;

% some special terms we'll need throughout
wa=wm;
kya=wm*pi/b;
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
lamztha=kt*sqrt(  1 - ( ( lamx.^2 + kya.^2)./ktz^2 )  );
lamzpsia=kt*sqrt(  1 - ( ( lamx.^2 + kya.^2)./kzt^2 )  );


% % % A11lamy term
% % tez terms
% the sum term
Aterm1= -((4*pi*mu0*muz.*lamx.^2.*Clamx)./(w*mu0^2*mut^2*d)).*...
    sum( ( ( (pi*l/d).^2 ).*(lamylth).*(1-exp(-1j*lamylth*b) ) )./...
    ( (lamxl.^2 + lamylth.^2).*(lamylth.^2 - kya.^2).^2 ),1);

% the second term
Aterm2stable=( (1j*pi*b.*lamztha.*lamx.^2.*Clamx)./( w*mut*mu0.*(lamx.^2 + kya.^2) ) );
Aterm2unstable=(  ( cos(lamztha*d) )./( sin(lamztha*d) )  );
% if any(isnan(Aterm2unstable))==1
%     display 'found nan'
% end
Aterm2unstable(isnan(Aterm2unstable))=1j;
Aterm2=Aterm2stable.*Aterm2unstable;

% % tmz terms
% the sum term 
Aterm3=-( (4.*pi.*w.*ez.*eps0.*Clamx)./(d) ).*sum( ( (lamylpsi.^3).*( 1-exp(-1j.*lamylpsi.*b) ) )./...
    ( (lamylpsi.^2 + lamxl.^2).*(lamylpsi.^2 - kya.^2).^2.*( 1+1.*(l==0) ) ),1);

% the fourth term
Aterm4stable=(  (1j*pi*b.*w.*et.*eps0.*kya.^2.*Clamx)./( (lamx.^2 + kya.^2) )  );
Aterm4unstable=(  (cos(lamzpsia*d) )./( lamzpsia.*sin(lamzpsia*d) )  );
% if any(isnan(Aterm4unstable))==1
%     display 'found nan'
% end
Aterm4unstable(isnan(Aterm4unstable))=1j;
Aterm4=Aterm4stable.*Aterm4unstable;

% put it all together
Alamy=( (Mhxm.*Mhxn.*vm.*vn)./ (a^2) ).*(Aterm1+Aterm2+Aterm3+Aterm4);

% % % B11lamy term
% % tez terms
% the sum term
Bterm1= -( (4*pi*mu0*muz.*lamx.^2.*Clamx)./(w*mu0^2*mut^2*d) ).*...
    sum( ...
	(  ( (pi*l/d).^2 ).*(lamylth).*( 1-exp(-1j*lamylth*b) )  )./...
      ( (lamxl.^2 + lamylth.^2).*(lamylth.^2 - kya.^2).^2 )...
      ,1);

% the second term
Bterm2stable=(  (1j*pi*b.*lamztha.*lamx.^2.*Clamx)./( w*mut*mu0.*(lamx.^2 + kya.^2) )  );
Bterm2unstable=( cos(lamztha*d) )./( sin(lamztha*d) ) ;
% if any(isnan(Bterm2unstable))==1
%     display 'found nan'
% end
Bterm2unstable(isnan(Bterm2unstable))=1j;
Bterm2=Bterm2stable.*Bterm2unstable;

% % tmz terms
% the sum term 
Bterm3= (  (4.*pi.*w.*ez.*eps0.*lamx.^2.*Clamx)./(d) ).*...
		sum( ...
			(  (lamylpsi).*( 1-exp(-1j.*lamylpsi.*b) )  )./...
    			(  (lamylpsi.^2 + lamxl.^2).*(lamylpsi.^2 - kya.^2).^2.*( 1+1.*(l==0) )  )...
            ,1);

% the fourth term
Bterm4stable= -( 1j*pi*b.*w.*et.*eps0.*lamx.^2.*Clamx)./( (lamx.^2 + kya.^2) ) ;
Bterm4unstable= ( cos(lamzpsia*d) )./( lamzpsia.*sin(lamzpsia*d) );
% if any(isnan(Bterm4unstable))==1
%     display 'found nan'
% end
Bterm4unstable(isnan(Bterm4unstable))=1j;
Bterm4=Bterm4stable.*Bterm4unstable;

% put it all together
Blamy=( (Mhxm.*Mhyn.*wa.*vm)./ (a*b) ).*(Bterm1+Bterm2+Bterm3+Bterm4);

% % % C11lamy term
% % tez terms
% the sum term
Cterm1= -( (4*pi*mu0*muz.*lamx.^2.*Clamx)./(w*mu0^2*mut^2*d) ).*...
    sum( ...
	(  ( (pi*l/d).^2 ).*(lamylth).*( 1-exp(-1j*lamylth*b) )  )./...
      ( (lamxl.^2 + lamylth.^2).*(lamylth.^2 - kya.^2).^2 )...
       ,1);

% the second term
Cterm2stable=(  (1j*pi*b.*lamztha.*lamx.^2.*Clamx)./( w*mut*mu0.*(lamx.^2 + kya.^2) )  );
Cterm2unstable=( cos(lamztha*d) )./( sin(lamztha*d) ) ;
% if any(isnan(Cterm2unstable))==1
%     display 'found nan'
% end
Cterm2unstable(isnan(Cterm2unstable))=1j;
Cterm2=Cterm2stable.*Cterm2unstable;

% % tmz terms
% the sum term 
Cterm3= (  (4.*pi.*w.*ez.*eps0.*lamx.^2.*Clamx)./(d) ).*...
		sum( ...
			(  (lamylpsi).*( 1-exp(-1j.*lamylpsi.*b) )  )./...
    			(  (lamylpsi.^2 + lamxl.^2).*(lamylpsi.^2 - kya.^2).^2.*( 1+1.*(l==0) )  ),1 ...
		    );

% the fourth term
Cterm4stable= -( 1j*pi*b.*w.*et.*eps0.*lamx.^2.*Clamx)./( (lamx.^2 + kya.^2) ) ;
Cterm4unstable= ( cos(lamzpsia*d) )./( lamzpsia.*sin(lamzpsia*d) );
% if any(isnan(Cterm4unstable))==1
%     display 'found nan'
% end
Cterm4unstable(isnan(Cterm4unstable))=1j;
Cterm4=Cterm4stable.*Cterm4unstable;

% put it all together
Clamy=( (Mhym.*Mhxn.*wa.*vn)./ (a*b) ).*(Cterm1+Cterm2+Cterm3+Cterm4);


% % % D11lamy term
% % tez terms
% the sum term
Dterm1= -( (4*pi*mu0*muz.*lamx.^2.*Clamx)./(w*mu0^2*mut^2*d) ).*...
    sum( ...
	(  ( (pi*l/d).^2 ).*(lamylth).*( 1-exp(-1j*lamylth*b) )  )./...
      ( (lamxl.^2 + lamylth.^2).*(lamylth.^2 - kya.^2).^2 )...
      ,1);

% the second term
Dterm2stable=(  (1j*pi*b.*lamztha.*lamx.^2.*Clamx)./( w*mut*mu0.*(lamx.^2 + kya.^2) )  );
Dterm2unstable=( cos(lamztha*d) )./( sin(lamztha*d) ) ;
% if any(isnan(Dterm2unstable))==1
%     display 'found nan'
% end
Dterm2unstable(isnan(Dterm2unstable))=1j;
Dterm2=Dterm2stable.*Dterm2unstable;

% % tmz terms
% the sum term 
Dterm3= - (  (4.*pi.*w.*ez.*eps0.*lamx.^4.*Clamx)./(d) ).*...
		sum( ...
			(  ( 1-exp(-1j.*lamylpsi.*b) )  )./...
    			(  (lamylpsi).*(lamylpsi.^2 + lamxl.^2).*(lamylpsi.^2 - kya.^2).^2.*( 1+1.*(l==0) )  )...
            ,1);

% the fourth term
Dterm4stable= ( 1j*pi*b.*w.*et.*eps0.*lamx.^4.*Clamx)./( kya.^2.*(lamx.^2 + kya.^2) ) ;
Dterm4unstable= ( cos(lamzpsia*d) )./( lamzpsia.*sin(lamzpsia*d) );
% if any(isnan(Dterm4unstable))==1
%     display 'found nan'
% end
Dterm4unstable(isnan(Dterm4unstable))=1j;
Dterm4=Dterm4stable.*Dterm4unstable;

% put it all together
Dlamy=( (Mhym.*Mhyn.*wa^2)./ (b^2) ).*(Dterm1+Dterm2+Dterm3+Dterm4);

y=Alamy+Blamy+Clamy+Dlamy;

end








