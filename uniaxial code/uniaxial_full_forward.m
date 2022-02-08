%%% full mode uniaxial constitutive parameter extraction
%% configuration/constants

clear all;
clc;
close all;

global widx v wi a b dmat eps0 mu0 numModes includeModes solveCase porttouse alg ket kmut num_int; 

% % % % % % % % % % % % % % % % % % % %
% % constants & global parameters
% % % % % % % % % % % % % % % % % % % %

v=[1 3 1 1 5 3 3 5 5 7 7 7 9 1 1 3 3 9 9 5]; % the first 20 modes
wi=[0 0 2 2 0 2 2 2 2 0 2 2 0 4 4 4 4 2 2 4];
c=3.0e8;           
eps0=8.854e-12;
mu0=pi*4e-7;
ket=[];
kmut=[];
a=0.9*2.54/100;  % inches to m
b=0.4*2.54/100;  % inches to m
wmin=8.2e9*2*pi;
wmax=12.4e9*2*pi;


% % % % % % % % % % % % % % % % % % % %
% % user flags & config options
% % % % % % % % % % % % % % % % % % % %
% for forward calcs
use_real_data=1;
add_noise=1;
noise_snr=50; % in dB

mat_et_real_mult = 5;
mat_et_imag_mult = 1;
mat_ez_real_mult = 5;
mat_ez_imag_mult = 1;


% for reverse calcs
material='high_permittivity_sim';
vers='v1'; % the version of the measurement
% orientation='_O2'; % needs the preceeding '_'
methd='tfwgt'; % tfwgt or nrw
add2casedesc='_forward';  % a string to add to the case description (put leading underscore)
porttouse=3; % 1=use S11 & S21, 2=use S22 & S12, 3=use all
includeModes=1; % dominant mode only
%  includeModes=[1 3 4 14 15]; % indices of the top 5 modes - 99% of the solution
%  includeModes=[1 3 4];
self_cal=0; % 1=provide your own TRL cal files, 0=use VNA cal'd data
errchk=0;    % error bar routine? 1=yes, 0=no

% % solver options 
numds=25;  % downsample the data to this number of points, 0 equates to no downsampling
init_method='initup'; % 'nrw','nrwTruth','initone', or 'initup'
nrwvers='v6'; % version of nrw measurements, if using them
solveCase=4; % this is a switch to test some different cases, where we force
             %  certain symmetries or conditions on eps & mu
             % 1 - isotropic, dielectric, non-magnetic (et=ez=er) & (mut=muz=mu0)
             % 2 - isotropic, non-dielectric, magnetic (et=ez=eps0) & (mut=muz=mur)
             % 3 - isotropic, dielectric, magnetic (et=ez=er) & (mut=muz=mur)
             % 4 - uniaxial, dielectric, non-magnetic (et,ez) & (mut=muz=mu0)
             % 5 - uniaxial, non-dielectric, magnetic (et=ez=eps0) & (mut,muz)
             % 6 - uniaxial, dielectric, magnetic (et,ez) & (mut,muz)
             % 7 - uniaxial with known et and mut (et and mut set to
             %       etguess and mutguess, then solve for ez and muz)
alg='TRR';  % LM=levenburg-marquardt, TRR=Trust-region-reflective
num_int=0;  % use numerical integration for lamx and lamy 


% % post-processing options
makeFile=1;  % output the results to a csv (*.dat) file? 1=yes, 0=no
print=1;

% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
% % get details for the material specified above
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
dmat=[0.4 0.8].*2.54/100;          % for 0.4" and 0.8" layer of material
% initial guesses for reverse solver
myret=mat_et_real_mult; 
myiet=-mat_et_imag_mult;
myrmut=1;  
myimut=-0.05;
myrez=myret;
myiez=myiet;
myrmuz=myrmut;
myimuz=myimut;

numModes=length(includeModes); % number of modes to include

if use_real_data==1
    % use [VARNAME1,VARNAME2,VARNAME3] = importdatafile(FILENAME);
    [wfds,mat_et_real,mat_et_real_stddev] = importdatfile('ui_honeycomb_0_400_solveCase_4_tfwgt_1mode_et_real_paper_final.dat');
    [~,mat_et_imag,mat_et_imag_stddev] = importdatfile('ui_honeycomb_0_400_solveCase_4_tfwgt_1mode_et_imag_paper_final.dat');
    [~,mat_ez_real,mat_ez_real_stddev] = importdatfile('ui_honeycomb_0_400_solveCase_4_tfwgt_1mode_ez_real_paper_final.dat');
    [~,mat_ez_imag,mat_ez_imag_stddev] = importdatfile('ui_honeycomb_0_400_solveCase_4_tfwgt_1mode_ez_imag_paper_final.dat');
    mat_et_real = mat_et_real.*mat_et_real_mult;
    mat_et_imag = mat_et_imag.*(mat_et_imag_mult);
    mat_ez_real = mat_ez_real.*mat_ez_real_mult;
    mat_ez_imag = mat_ez_imag.*(mat_ez_imag_mult);
    wfds=wfds.*1e9.*2.*pi; % the input files are normalized to 1e9 and in f
else
    % calculate frequencies
    wfds=linspace(wmin,wmax,numds); 
    % since we are running the forward case, we must define the material we are
    %   trying to "simulate" (for entire troutine, params are normalized to
    %   eps0 and mu0)
    mat_et_real = ones(length(wfds),1).*10;
    mat_et_imag = ones(length(wfds),1).*(-0.05);
    mat_ez_real = ones(length(wfds),1).*8;
    mat_ez_imag = ones(length(wfds),1).*(-0.01);
end

mat_mut_real = ones(length(wfds),1).*1;
mat_mut_imag = ones(length(wfds),1).*0;
mat_muz_real = ones(length(wfds),1).*1;
mat_muz_imag = ones(length(wfds),1).*0;
fds=wfds./(1e9.*2.*pi);
% guess method = one initial guess + update as we go

etguess=zeros(length(wfds)+1,1);
ezguess=zeros(length(wfds)+1,1);
mutguess=zeros(length(wfds)+1,1);
muzguess=zeros(length(wfds)+1,1);
etguess(1)=(myret+1j*myiet);
ezguess(1)=(myrez+1j*myiez);
mutguess(1)=(myrmut+1j*myimut);
muzguess(1)=(myrmuz+1j*myimuz);

% start recording
diary on;
display('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')
display(['Material is: ' material]);
dmat
display('(in m)')
display('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')
%%  solution routine

% % % % % % % % % % % % % % % % % % % %
% % Calculation of A coefficients and scattering parameters
% % % % % % % % % % % % % % % % % % % %

% pre-allocate the solution vectors
etsol=zeros(length(wfds),1);
ezsol=zeros(length(wfds),1);
mutsol=zeros(length(wfds),1);
muzsol=zeros(length(wfds),1);
tsolvefwd=zeros(length(wfds),1);
tsolverev=zeros(length(wfds),1);
eigerror=zeros(length(wfds)+1,1);
if porttouse == 3
    Scalc=zeros(length(dmat).*4.*2,length(wfds)); % real, imag parts for S11, S12, S21, S22
else
    Scalc=zeros(length(dmat).*2.*2,length(wfds)); % real, imag parts for S11 and S12 or S21 and S22
end
Y=0;  % initialize the solution array
if errchk==1
    stddeltaetreal=zeros(length(wfds),1);
    stddeltaetimag=zeros(length(wfds),1);
    stddeltaezreal=zeros(length(wfds),1);
    stddeltaezimag=zeros(length(wfds),1);
    stddeltamutreal=zeros(length(wfds),1);
    stddeltamutimag=zeros(length(wfds),1);
    stddeltamuzreal=zeros(length(wfds),1);
    stddeltamuzimag=zeros(length(wfds),1);
end

% calculate the forward s-params
for widx=1:length(wfds) % solve at each frequency
    tic;
    wval=wfds(widx);
    display(sprintf('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%'))
    display(sprintf(['Now processing point ' num2str(widx) '/' num2str(length(wfds))]))
    display(sprintf(['Frequency = ' num2str(wval/2/pi) 'GHz']))
    display(['et is ' num2str(mat_et_real(widx)) ' ' num2str(mat_et_imag(widx)) 'i'])
    display(['ez is ' num2str(mat_ez_real(widx)) ' ' num2str(mat_ez_imag(widx)) 'i'])
    display(['mut is ' num2str(mat_mut_real(widx)) ' - ' num2str(mat_mut_imag(widx)) 'i'])
    display(['muz is ' num2str(mat_muz_real(widx)) ' - ' num2str(mat_muz_imag(widx)) 'i'])
    display(sprintf('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%'))
    display(sprintf(' \n'))
    
% fwd calc
    Scalc(:,widx) = Sparams([mat_et_real(widx) mat_et_imag(widx) ...
                            mat_ez_real(widx) mat_ez_imag(widx) ...
                            mat_mut_real(widx) mat_mut_imag(widx) ...
                            mat_muz_real(widx) mat_muz_imag(widx)],wval);
end

if add_noise==1
 % add noise to S11, S21, S12, S22 (real and imag parts) individually for
 % all freqs
     for nidx = 1:size(Scalc,1)
%          Smeas(nidx,:)=Scalc(nidx,:)+ randn(1,size(Scalc(nidx,:),2))*std(Scalc(nidx,:));
%          Smeas(nidx,:)=add_awgn_noise(Scalc(nidx,:),noise_snr); % gaussian white
%         Smeas(nidx,:) =   % Gaussian 
     end
 % plot Scalc and Smeas (with noise)
 myfigure
 title('S11 Param with noise')
 hold on;
 plot(fds, Smeas(1,:),fds,Scalc(1,:))
else
    Smeas=Scalc;
end
 
for widx=1:length(wfds)
   wval=wfds(widx);
    % calculate eps and mu using calculated s-params as measured s-params
   [etsol(widx), ezsol(widx), mutsol(widx), muzsol(widx)] ...
        = runSolver(Smeas(:,widx),wval,etguess(widx),ezguess(widx),mutguess(widx),muzguess(widx));
    
    if strcmp(init_method,'initup')==1
        etguess(widx+1)=etsol(widx);
        ezguess(widx+1)=ezsol(widx);
        mutguess(widx+1)=mutsol(widx);
        muzguess(widx+1)=muzsol(widx);
    end
    
    display('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')
    display(['Solution for point ' num2str(widx) '/' num2str(length(wfds)) ' is: '])
    display(['Et (known) = ' num2str(mat_et_real(widx)) ' ' num2str(mat_et_imag(widx)) 'i' ])
    display(['Et (calc) = ' num2str(etsol(widx))])
    display(['Ez (known) = ' num2str(mat_ez_real(widx)) ' ' num2str(mat_ez_imag(widx)) 'i'])
    display(['Ez (calc) = ' num2str(ezsol(widx))])
    display('(in m)')
    display('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')
%     pause
    tsolvefwd(widx)=toc;
end

% plot the results
myfigure
hold on;
title('Et Comparison')
plot(fds,mat_et_real,fds,real(etsol),fds,mat_et_imag,fds,imag(etsol));
legend('Real Et (known)','Real Et (calc)','Imag Et (known)','Imag Et (calc)')

myfigure
hold on;
title('Ez Comparison')
plot(fds,mat_ez_real,fds,real(ezsol),fds,mat_ez_imag,fds,imag(ezsol));
legend('Real Ez (known)','Real Ez (calc)','Imag Ez (known)','Imag Ez (calc)')