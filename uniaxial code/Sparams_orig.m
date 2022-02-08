function Sthy=Sparams(X,wval)
% export the theoretical scattering coefficients with the real part on top
% and the imaginary part below
% INPUTS:
% 1) the vector X contains the input arguments in the order:
% - real(et), imag(et), real(ez), imag(ez)
% - real(mut), imag(mut), real(muz), imag(muz)
% 2) wval is the angular frequency (single value)
% % Calculation of A coefficients and scattering parameters

global widx v wi a b dmat eps0 mu0 numModes includeModes solveCase;

k0=sqrt(wval.^2.*eps0.*mu0); % gonna need this later

% pre-allocate for speed
Sthy=0;
A11=zeros(numModes,numModes);
A12=A11;
Bmat=zeros(2*numModes,1);
% B2=zeros(numModes,1);

% cases:
%  1 - isotropic, dielectric, non-magnetic (et=ez=er) & (mut=muz=mu0)
%  2 - isotropic, non-dielectric, magnetic (et=ez=eps0) & (mut=muz=mur)
%  3 - isotropic, dielectric, magnetic (et=ez=er) & (mut=muz=mur)
%  4 - uniaxial, dielectric, non-magnetic (et,ez) & (mut=muz=mu0)
%  5 - uniaxial, non-dielectric, magnetic (et=ez=eps0) & (mut,muz)
%  6 - uniaxial, dielectric, magnetic (et,ez) & (mut,muz) 
switch solveCase
    case 1
        ret=X(1);
        iet=X(2);
        rez=ret;
        iez=iet;
        rmut=1;
        imut=0;
        rmuz=rmut;
        imuz=imut;
    
    case  2 
        ret=1;
        iet=0;
        rez=ret;
        iez=iet;
        rmut=X(1);
        imut=X(2);
        rmuz=rmut;
        imuz=imut;
    
    case 3
        ret=X(1);
        iet=X(2);
        rez=ret;
        iez=iet;
        rmut=X(3);
        imut=X(4);
        rmuz=rmut;
        imuz=imut;
    
    case 4
        ret=X(1);
        iet=X(2);
        rez=X(3);
        iez=X(4);
        rmut=1;
        imut=0;
        rmuz=rmut;
        imuz=imut;
        
    case 5
        ret=1;
        iet=0;
        rez=ret;
        iez=iet;
        rmut=X(1);
        imut=X(2);
        rmuz=X(3);
        imuz=X(4);
    case 6 
        ret=X(1);
        iet=X(2);
        rez=X(3);
        iez=X(4);
        rmut=X(5);
        imut=X(6);
        rmuz=X(7);
        imuz=X(8);
end

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
    midx=0;
    for m=includeModes % source modes
        midx=midx+1;
        nidx=0;
        vm=v(m); 
        wm=wi(m);
        if wm==0
            delta0wm = 1;
        else
            delta0wm = sqrt(2);
        end
        
      % waveguide parameters for m
        kxm=vm.*pi/a;
        kym=wm.*pi/b;
        kcm=sqrt(kxm.^2+kym.^2);
        kzm=sqrt(k0.^2-kcm.^2);
                 
      % m index normalization coefficients
        if any(m==[4 7 9 15 19]) % tmz
            Zm=kzm./(wval.*eps0);  % TMZ m index
            Mhxm=(sqrt(2).*kym.*delta0wm)./(Zm.*kcm*sqrt(a*b));  % TMZ
            Mhym=-(sqrt(2).*kxm.*delta0wm)./(Zm.*kcm*sqrt(a*b));  % TMZ
        else % it's tez
            Zm=wval.*mu0./kzm;  % TEZ m index
            Mhxm=(sqrt(2).*kxm.*delta0wm)./(Zm.*kcm*sqrt(a*b));  % TEZ
            Mhym=(sqrt(2).*kym.*delta0wm)./(Zm.*kcm*sqrt(a*b));  % TEZ
        end
        
        for n=includeModes  % observation modes
            nidx=nidx+1;
         % the source terms for the MFIE's  
            if m==1 && n==1 % only the dominant source is excited
                Bmat(1)=(2./(Zm.^2)); 
            end
            vn=v(n); 
            wn=wi(n);
            if wn==0
                delta0wn = 1;
            else
                delta0wn = sqrt(2);
            end
            
            deltamn=1.*(m==n); 
            
            % waveguide parameters for n
            kxn=vn.*pi/a;
            kyn=wn.*pi/b; 
            kcn=sqrt(kxn.^2+kyn.^2);
            kzn=sqrt(k0.^2-kcn.^2);
            
            % n index normalization coefficients
            if any(n==[4 7 9 15 19]) 
              Zn=kzn./(wval.*eps0);  % TMZ n index  
              Mhxn=(sqrt(2).*kyn.*delta0wn)./(Zn.*kcn*sqrt(a*b));  % TMZ
              Mhyn=-(sqrt(2).*kxn.*delta0wn)./(Zn.*kcn*sqrt(a*b));  % TMZ
            else
              Zn=wval.*mu0./kzn;  %  TEZ n index         
              Mhxn=(sqrt(2).*kxn.*delta0wn)./(Zn.*kcn*sqrt(a*b));  % TEZ
              Mhyn=(sqrt(2).*kyn.*delta0wn)./(Zn.*kcn*sqrt(a*b));  % TEZ
            end
                      
          % % which set of lamy solutions we use depends on wm and wn 
          if wm == 0 && wn == 0  % Case 1
               % quadgk numerically solves the lamx integrals
%                 disp(['v_m = ' num2str(vm)  ' and v_n = ' num2str(vn)])
%                 disp(['w_m = ' num2str(wm)  ' and w_n = ' num2str(wn)])
%                 disp('Using Case 1')
                intval11=quadgk(@(lamx) SelfIntegral_1(wval,et,ez,mut,muz,vm,vn,a,b,d,Mhxm,Mhxn,lamx),0,50e3,'abstol',1e-10,'reltol',1e-10);
                intval12=quadgk(@(lamx) CouplingIntegral_1(wval,et,ez,mut,muz,vm,vn,a,b,d,Mhxm,Mhxn,lamx),0,50e3,'abstol',1e-10,'reltol',1e-10);
           elseif wm ~=0 && wn==0 % Case 2
%                 disp(['v_m = ' num2str(vm)  ' and v_n = ' num2str(vn)])
%                 disp(['w_m = ' num2str(wm)  ' and w_n = ' num2str(wn)])
%                 disp('Using Case 2')
                intval11=quadgk(@(lamx) SelfIntegral_2(wval,et,ez,mut,muz,wm,vm,vn,a,b,d,Mhxm,Mhxn,Mhym,lamx),0,50e3,'abstol',1e-10,'reltol',1e-10);
                intval12=quadgk(@(lamx) CouplingIntegral_2(wval,et,ez,mut,muz,wm,vm,vn,a,b,d,Mhxm,Mhxn,Mhym,lamx),0,50e3,'abstol',1e-10,'reltol',1e-10);
           elseif wm==0 && wn~=0  % Case 3 
%                 disp(['v_m = ' num2str(vm)  ' and v_n = ' num2str(vn)])
%                 disp(['w_m = ' num2str(wm)  ' and w_n = ' num2str(wn)])
%                 disp('Using Case 3')
                intval11=quadgk(@(lamx) SelfIntegral_3(wval,et,ez,mut,muz,wn,vm,vn,a,b,d,Mhxm,Mhxn,Mhyn,lamx),0,50e3,'abstol',1e-10,'reltol',1e-10);
                intval12=quadgk(@(lamx) CouplingIntegral_3(wval,et,ez,mut,muz,wn,vm,vn,a,b,d,Mhxm,Mhxn,Mhyn,lamx),0,50e3,'abstol',1e-10,'reltol',1e-10);    
           elseif wm~=0 && wn~=0 && wm==wn % Case 4
%                disp(['v_m = ' num2str(vm)  ' and v_n = ' num2str(vn)])
%                disp(['w_m = ' num2str(wm)  ' and w_n = ' num2str(wn)])
%                disp('Using Case 4')
               intval11=quadgk(@(lamx) SelfIntegral_4(wval,et,ez,mut,muz,wm,vm,vn,a,b,d,Mhxm,Mhxn,Mhym,Mhyn,lamx),0,50e3,'abstol',1e-10,'reltol',1e-10); 
               intval12=quadgk(@(lamx) CouplingIntegral_4(wval,et,ez,mut,muz,wm,vm,vn,a,b,d,Mhxm,Mhxn,Mhym,Mhyn,lamx),0,50e3,'abstol',1e-10,'reltol',1e-10); 
           elseif wm~=0 && wn~=0 && wm~=wn % Case 5
%                disp(['v_m = ' num2str(vm)  ' and v_n = ' num2str(vn)])
%                disp(['w_m = ' num2str(wm)  ' and w_n = ' num2str(wn)])
%                disp('Using Case 5')
               intval11=quadgk(@(lamx) SelfIntegral_5(wval,et,ez,mut,muz,wm,wn,vm,vn,a,b,d,Mhxm,Mhxn,Mhym,Mhyn,lamx),0,50e3,'abstol',1e-10,'reltol',1e-10); 
               intval12=quadgk(@(lamx) CouplingIntegral_5(wval,et,ez,mut,muz,wm,wn,vm,vn,a,b,d,Mhxm,Mhxn,Mhym,Mhyn,lamx),0,50e3,'abstol',1e-10,'reltol',1e-10); 
          end  
        % the A & B coefficients 
          A11(midx,nidx)=deltamn./(Zm.*Zn) - ( (Zn)./(4).* 2.*intval11 );
          A12(midx,nidx) =  ((Zn)./(4)).* 2.*intval12;
          A21=A12;
          A22=A11;
        end % end of n loop
    end % end of m loop
    
  % solve for the C matrix
    Amat=[A11 A12; A21 A22];
%     Bmat=[B1 ; B2];
    Cmat=Amat\Bmat;

  % scattering parameters
    S11thy=Cmat(1)-1;
    S21thy=Cmat(numModes+1);
    S22thy=S11thy;
    S12thy=S21thy;
    Sthy=cat(1,Sthy,real(S11thy),imag(S11thy),real(S21thy),imag(S21thy),real(S12thy),imag(S12thy),...
          real(S22thy),imag(S22thy)); % real, then imag of each S parameter
end % end of dmat loop
Sthy=Sthy(2:end); % strip off the leading zero


end % end of function

    