%%% dominant mode uniaxial permittivity extraction
% things to do to add more modes:
% - the vectors that depend on m and n will have to be indexed, not equated
% profile on;
clear all;
clc;
myway=1;
plots=1;
% close all;

%% configuration/constants
% input file - the frequency ranges should be equal, or there's a problem
% f=fileformat3('meas1','fgm125_1.txt');
% f=fileformat3('meas2','fgm125_2.txt'); 
f=fileformat3('meas1','plexiglass.txt');

% some test values
myet=2.7997 - 0.016298i;
myez=myet;
mymu=1;
mymut=1;
mymuz=1;

% define constants
global a b dmat vm vn m n eps0 mu0 deltamn Zntez Zmtez Zntmz Zmtmz Mhxmtez;
global Mhxntez Mhymtez Mhyntez Mhxmtmz Mhxntmz Mhymtmz Mhyntmz;
m=[1];
n=[1];
v=[1 3 1 1 5 3 3 5 5 7 7 7 9 1 1 3 3 9 9 5]; % the first 20 modes
wi=[0 0 2 2 0 2 2 2 2 0 2 2 0 4 4 4 4 2 2 4];
vm=v(m);  % for dominant mode, will need to vectorize later on 
vn=v(n);  % for dominant mode, will need to vectorize later on 
wm=wi(m);  % for dominant mode, will need to vectorize later on 
wn=wi(n);  % for dominant mode, will need to vectorize later on 
c=3.0e8;           
eps0=8.854e-12;
mu0=pi*4e-7;
a=0.9*2.54/100;  % inches to m
b=0.4*2.54/100;  % inches to m
dmat=0.260.*2.54/100; % inches to m

% waveguide parameters
kxm=vm*pi/a;
kym=wm*pi/b;
kxn=vn*pi/a;
kyn=wn*pi/b;
wf=2*pi*f;  % for all freq's
% wf=2*pi*f(1); % for 1 freq
% wf=2*pi*f(1:2); % for 2 freqs
% wf=wf'; % because i like column vectors
k0=sqrt(wf.^2.*eps0.*mu0);
kcm=sqrt(kxm^2+kym^2);
kcn=sqrt(kxn^2+kyn^2);
kzm=sqrt(k0.^2-kcm.^2);
kzn=sqrt(k0.^2-kcn.^2);

if wn==0
    delta0wn = 1;
else
    delta0wn = sqrt(2);
end

if wm==0
    delta0wm = 1;
else
    delta0wm = sqrt(2);
end

Zmtez=wf.*mu0./kzm;  % TEZ m index
Zntez=wf.*mu0./kzn;  %  TEZ n index

Zmtmz=kzm./(wf.*eps0);  % TMZ m index 
Zntmz=kzn./(wf.*eps0);  % TMZ n index

% m index normalization coefficients
Mhxmtez=(sqrt(2).*kxm.*delta0wm)./(Zmtez.*kcm*sqrt(a*b));  % TEZ
Mhymtez=(sqrt(2).*kym.*delta0wm)./(Zmtez.*kcm*sqrt(a*b));  % TEZ

Mhxmtmz=(sqrt(2).*kym.*delta0wm)./(Zmtez.*kcm*sqrt(a*b));  % TMZ
Mhymtmz=-(sqrt(2).*kxm.*delta0wm)./(Zmtez.*kcm*sqrt(a*b));  % TMZ

% n index normalization coefficients
Mhxntez=(sqrt(2).*kxn.*delta0wn)./(Zntez.*kcn*sqrt(a*b));  % TEZ
Mhyntez=(sqrt(2).*kyn.*delta0wn)./(Zntez.*kcn*sqrt(a*b));  % TEZ

Mhxntmz=(sqrt(2).*kyn.*delta0wn)./(Zntez.*kcn*sqrt(a*b));  % TMZ
Mhyntmz=-(sqrt(2).*kxn.*delta0wn)./(Zntez.*kcn*sqrt(a*b));  % TMZ

deltamn=1.*(m==n); 

%% Calculation of A coefficients and scattering parameters
% what i'm about to do will only work for the dominant mode - will have to
% rework a new for loop for multiple values of m and n

% pre-allocate for speed
S11thy=zeros(length(wf),1);
S21thy=zeros(length(wf),1);

% call the functions and do the integrations in one loop

if myway==1
%     % my way
%     for widx=1:length(wf)
%         wval=wf(widx);
%     Sthy=Sparams(wf,myet,myez,mymut,mymuz,vm,vn,a,b,dmat);
%     Sthy=Sparams(wf,myet,myez,mymut,mymuz);
          Sthy=Sparams({myet,myez,mymut,mymuz},wf);
%         intval11=quadgk(@(lamx) SelfIntegral(wval,myet,myet,1,1,vm,vn,a,b,d,lamx),0,inf);
%         intval12=quadgk(@(lamx) CouplingIntegral(wval,myet,myet,1,1,vm,vn,a,b,d,lamx),0,inf);
% 
%         %the A & B coefficients
%         A11=deltamn./(Zmtez(widx).*Zntez(widx)) - ((1j.*Zntez(widx).*Mhxmtez(widx).*Mhxntez(widx).*vm.*vn)./...
%             (4.*a.^2)).* 2.*intval11;
%         A12 =  ((1j.*Zntez(widx).*Mhxmtez(widx).*Mhxntez(widx).*vm.*vn)./...
%             (4.*a.^2)).* 2.*intval12;
%         A21=A12;
%         A22=A11;
%         B1=(2./(Zmtez(widx).^2)).*(m==1);
%         B2=0;
% 
%         % solve for the C matrix
%         Amat=[A11 A12; A21 A22];
%         Bmat=[B1 ; B2];
%         Cmat=Amat\Bmat;
% 
%         % scattering parameters
%         S11thy(widx)=Cmat(1)-1;
%         S21thy(widx)=Cmat(2);
%     end
    if plots==1
        S11thy=Sthy(1:201);
        S21thy=Sthy(202:402);
        figure;
        subplot(1,2,1)
        plot(wf,real(S11thy),wf,imag(S11thy))
        legend('Re{S11}','Im{S11}')
        xlabel('\omega (GHz)')
        subplot(1,2,2)
        plot(wf,real(S21thy),wf,imag(S21thy))
        legend('Re{S21}','Im{S21}')
        xlabel('\omega (GHz)')
    end
elseif myway==0
    % the Hyde way
    for widx=1:length(wf)
        wval=wf(widx);
        intval11=quadgk(@(lamx) Self_Integral(lamx,{eps0*myet,eps0*myet,mu0*mymu,mu0*mymu,wval,pi/a,a,b,d}),0,inf);
        intval12=quadgk(@(lamx) Coupling_Integral(lamx,{eps0*myet,eps0*myet,mu0*mymu,mu0*mymu,wval,pi/a,a,b,d}),0,inf);
        HA11=(a*b/Zntez(widx))+( ( (1j*kxn^2)/(2*pi^2) )*2*intval11);
        HA12=-( (1j*kxn^2)/(2*pi^2) )*2*intval12;
        HA21=-HA12;
        HA22=-HA11;
        HAmat=[HA11 HA12; HA21 HA22];
        HB1=( (a*b)/Zntez(widx) ) - ( ( (1j*kxn^2)/(2*pi^2) )*2*intval11);
        HB2=-( (1j*kxn^2)/(2*pi^2) )*2*intval12;
        HBmat=[HB1; HB2];
        HCmat=HAmat\HBmat;
        HS11thy(widx)=HCmat(1);
        HS21thy(widx)=HCmat(2);
    end
    % plot the scattering coefficients
    if plots==1
        figure;
        subplot(1,2,1)
        plot(wf,real(HS11thy),wf,imag(HS11thy))
        legend('Re{S11}','Im{S11}')
        xlabel('\omega (GHz)')
        subplot(1,2,2)
        plot(wf,real(HS21thy),wf,imag(HS21thy))
        legend('Re{S21}','Im{S21}')
        xlabel('\omega (GHz)')
    end
end
    


