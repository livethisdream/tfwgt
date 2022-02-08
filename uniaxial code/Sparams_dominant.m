function Sthy=Sparams_dominant(X,wval)
% export the theoretical scattering coefficients with the real part on top
% and the imaginary part below
% INPUTS:
% 1) the vector X contains the input arguments in the order:
% - real(et), imag(et), real(ez), imag(ez)
% - real(mut), imag(mut), real(muz), imag(muz)
% 2) wval is the angular frequency (single value)
% % Calculation of A coefficients and scattering parameters

% global myrez myiez myrmuz myimuz
global myway widx vm vn a b dmat m n eps0 mu0 deltamn Zntez Zmtez Zntmz Zmtmz ;
global Mhxmtez Mhxntez Mhymtez Mhyntez Mhxmtmz Mhxntmz Mhymtmz Mhyntmz;

% % pre-allocate for speed
Sthy=0;

% real and imaginary parts need to be considered separately
ret=X(1);
iet=X(2);
rez=X(3);
iez=X(4);
rmut=X(5);
imut=X(6);
rmuz=X(7);
imuz=X(8);

% % for a known ez and muz
% ret=X(1);
% iet=X(2);
% rez=myrez(widx);
% iez=myiez(widx);
% rmut=X(3);
% imut=X(4);
% rmuz=myrmuz(widx);
% imuz=myimuz(widx);


% put the real and imaginary parts together
et=ret+1j*iet;
ez=rez+1j*iez;
mut=rmut+1j*imut;
muz=rmuz+1j*imuz;

% % display some useful info
% display(['Current Values are:'])
% display(['widx = ' num2str(widx)])
% display(['et = ' num2str(et)]) 
% display(['ez = ' num2str(ez)])
% display(['mut = ' num2str(mut)])
% display(['muz = ' num2str(muz)])

for didx=1:length(dmat) % if we have two thicknesses
    d=dmat(didx);

    if myway==1
      % my spectral integrals
        intval11=quadgk(@(lamx) SelfIntegral(wval,et,ez,mut,muz,vm,vn,a,b,d,lamx),0,50e3,'abstol',1e-10,'reltol',1e-10);
        intval12=quadgk(@(lamx) CouplingIntegral(wval,et,ez,mut,muz,vm,vn,a,b,d,lamx),0,50e3,'abstol',1e-10,'reltol',1e-10);
    else    
      % maj hyde's spectral integrals
        intval11= -quadgk(@(lamx) Self_Integral(lamx,{eps0*et,eps0*ez,mu0*mut,mu0*muz,wval,pi/a,a,b,d}),0,50e3,'abstol',1e-10,'reltol',1e-10);         
        intval12= -quadgk(@(lamx) Coupling_Integral(lamx,{eps0*et,eps0*ez,mu0*mut,mu0*muz,wval,pi/a,a,b,d}),0,50e3,'abstol',1e-10,'reltol',1e-10);
    end

  % the A & B coefficients
    A11=deltamn./(Zmtez(widx).*Zntez(widx)) - ((1j.*Zntez(widx).*Mhxmtez(widx).*Mhxntez(widx).*vm.*vn)./...
        (4.*a.^2)).* 2.*intval11;
    A12 =  ((1j.*Zntez(widx).*Mhxmtez(widx).*Mhxntez(widx).*vm.*vn)./...
        (4.*a.^2)).* 2.*intval12;
    A21=A12;
    A22=A11;
    B1=(2./(Zmtez(widx).^2)).*(m==1);
    B2=0;

  % solve for the C matrix
    Amat=[A11 A12; A21 A22];
    Bmat=[B1 ; B2];
    Cmat=Amat\Bmat;

  % scattering parameters
    S11thy=Cmat(1)-1;
    S21thy=Cmat(2);
    S22thy=S11thy;
    S12thy=S21thy;
    Sthy=cat(1,Sthy,real(S11thy),imag(S11thy),real(S21thy),imag(S21thy),real(S12thy),imag(S12thy),...
          real(S22thy),imag(S22thy)); % real, then imag of each S parameter
end % end of dmat loop
Sthy=Sthy(2:end); % strip off the leading zero
end % end of function

    