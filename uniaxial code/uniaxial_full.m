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

material='ui_honeycomb_0_400';
vers='v1'; % the version of the measurement
orientation='_O2'; % needs the preceeding '_'
methd='tfwgt'; % tfwgt or nrw
add2casedesc='_recompute';  % a string to add to the case description (put leading underscore)
porttouse=3; % 1=use S11 & S21, 2=use S22 & S12, 3=use all
% includeModes=1; % dominant mode only
includeModes=[1 3 4 14 15]; % indices of the top 5 modes - 99% of the solution
% includeModes=[1 2 3 4];
self_cal=0; % 1=provide your own TRL cal files, 0=use VNA cal'd data
errchk=1;    % error bar routine? 1=yes, 0=no

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
smooth_zvals=0;
smooth_tvals=0;
smeth='sgolay';
sspan=201;

% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
% % get details for the material specified above
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %

% pick the correct thicknesses and initial guess based on the material
[myret, myiet, myrmut, myimut, myrez, myiez, myrmuz, myimuz, dmat]=material_opts(material);

display('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')
display(['Material is: ' material]);
dmat
display('(in m)')
display('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')


% check solveCase and dmat length and import relevant files
if ismember(solveCase,[3 6])==0
    inputFile=[material orientation '_' methd '_' vers '.txt'];
    dmat=dmat(1);
else
    if length(dmat)==1
        display('Error - dmat not long enough');
        return
    end
    inputFile=[material orientation '_' methd '_' vers '_d1.txt'];
    inputFile2=[material orientation '_' methd '_' vers '_d2.txt'];
end

numModes=length(includeModes); % number of modes to include

% % % % % % % % % % % % % % % % % % % %
% % setup for NRW analysis
% % % % % % % % % % % % % % % % % % % %

if strcmp(init_method,'nrw')==1
    nrwInput=[material '_nrw_' nrwvers '.txt'];
end

% % % % % % % % % % % % % % % % % % % %
% % other misc setup 
% % % % % % % % % % % % % % % % % % % %

% the truth data i have for this init method reqiures a slight adjustment to the
%  frequency values
if (strcmp(init_method,'nrwTruth')==1) && (any(regexp(material,'white_stuff')==1))
    wmax=12e9*2*pi;
end

% names for output files
if numModes==1
    casedesc=[material orientation '_' methd '_' vers '_' upper(init_method) 'Init_' num2str(numds) ' points_' num2str(numModes) ' mode_solveCase' num2str(solveCase) '_port' num2str(porttouse) add2casedesc];
    diary([upper(init_method) 'Init_' num2str(numds) 'points_' num2str(numModes) 'mode_solveCase_' num2str(solveCase) '_port' num2str(porttouse) add2casedesc '_log.txt']);
else
    casedesc=[material orientation '_' methd '_' vers '_' upper(init_method) 'Init_' num2str(numds) ' points_' num2str(numModes) ' modes_solveCase' num2str(solveCase) '_port' num2str(porttouse) add2casedesc];
    diary([upper(init_method) 'Init_' num2str(numds) 'points_' num2str(numModes) ' modes_solveCase_' num2str(solveCase) '_port' num2str(porttouse) add2casedesc '_log.txt']);
end

diary on;
      
%% get input data & setup initial guesses

% % % % % % % % % % % % % % % % % % % %
% %  Grab the input data & format it correctly
% % % % % % % % % % % % % % % % % % % %

if self_cal==0

%   % original import code - needs to have sparams in a sepecific order
%     [f,realS11meas1,imagS11meas1,realS21meas1,imagS21meas1,...
%         realS12meas1,imagS12meas1,realS22meas1,imagS22meas1]=fileformat4(inputFile);
      
    [A,svarnames]=fileformat5(inputFile);
    f=A(:,1);
      
    % use the svarnames to assign the correct name to the sparams
    eval(['real' cell2mat(svarnames(1)) 'meas1=A(:,2);']);
    eval(['imag' cell2mat(svarnames(1)) 'meas1=A(:,3);']);
    eval(['real' cell2mat(svarnames(2)) 'meas1=A(:,4);']);
    eval(['imag' cell2mat(svarnames(2)) 'meas1=A(:,5);']);
    eval(['real' cell2mat(svarnames(3)) 'meas1=A(:,6);']);
    eval(['imag' cell2mat(svarnames(3)) 'meas1=A(:,7);']);
    eval(['real' cell2mat(svarnames(4)) 'meas1=A(:,8);']);
    eval(['imag' cell2mat(svarnames(4)) 'meas1=A(:,9);']);
      
elseif self_cal==1
    % if you want to do your own TRL cal
    [f,S11meas1,S21meas1,S12meas1,S22meas1]=TRL('thru.txt','line.txt','reflect.txt',inputFile);
    realS11meas1=real(S11meas1);
    imagS11meas1=imag(S11meas1);
    realS21meas1=real(S21meas1);
    imagS21meas1=imag(S21meas1);
    realS12meas1=real(S12meas1);
    imagS12meas1=imag(S12meas1);
    realS22meas1=real(S22meas1);
    imagS22meas1=imag(S22meas1);
 end
    
wf=2*pi*f;  % for all freq's

% in case the data is outside of normal freq range, truncate it
minidx=find(wf>=wmin,1,'first');  
maxidx=find(wf>=wmax,1,'first');
wf=wf(minidx:maxidx);
realS11meas1=realS11meas1(minidx:maxidx);
imagS11meas1=imagS11meas1(minidx:maxidx);
realS21meas1=realS21meas1(minidx:maxidx);
imagS21meas1=imagS21meas1(minidx:maxidx);
realS12meas1=realS12meas1(minidx:maxidx);
imagS12meas1=imagS12meas1(minidx:maxidx);
realS22meas1=realS22meas1(minidx:maxidx);
imagS22meas1=imagS22meas1(minidx:maxidx);

% downsample the measured data, combine real & imag, transpose to columnns
if numds==0
    numds=length(wf);
end

wfds=linspace(wf(1),wf(end),numds); 
S11meas1=realS11meas1+1j.*imagS11meas1;
S11m1ds=interp1(wf,S11meas1,wfds);  % downsample
S11m1ds=S11m1ds(:);  % i like column vectors
S21meas1=realS21meas1+1j.*imagS21meas1; 
S21m1ds=interp1(wf,S21meas1,wfds);
S21m1ds=S21m1ds(:);
S12meas1=realS12meas1+1j.*imagS12meas1; 
S12m1ds=interp1(wf,S12meas1,wfds);
S12m1ds=S12m1ds(:);
S22meas1=realS22meas1+1j.*imagS22meas1;
S22m1ds=interp1(wf,S22meas1,wfds);  
S22m1ds=S22m1ds(:); 

% if using the two-thickness method, get second set of measurements
% if length(dmat)==2
if ismember(solveCase,[3 6])==1
    if self_cal==0  % use VNA cal'd data
%       % original import code - needs to have sparams in a specific order
%         [f,realS11meas2,imagS11meas2,realS21meas2,imagS21meas2,...
%            realS12meas2,imagS12meas2,realS22meas2,imagS22meas2]...
%            =fileformat4(inputFile2); % grab data from input file
        [A,svarnames]=fileformat5(inputFile2);      
        % use the svarnames to assign the correct name to the sparams
        eval(['real' cell2mat(svarnames(1)) 'meas2=A(:,2);']);
        eval(['imag' cell2mat(svarnames(1)) 'meas2=A(:,3);']);
        eval(['real' cell2mat(svarnames(2)) 'meas2=A(:,4);']);
        eval(['imag' cell2mat(svarnames(2)) 'meas2=A(:,5);']);
        eval(['real' cell2mat(svarnames(3)) 'meas2=A(:,6);']);
        eval(['imag' cell2mat(svarnames(3)) 'meas2=A(:,7);']);
        eval(['real' cell2mat(svarnames(4)) 'meas2=A(:,8);']);
        eval(['imag' cell2mat(svarnames(4)) 'meas2=A(:,9);']);
    elseif self_cal==1 % if you want to do your own TRL cal
        [f,S11meas2,S21meas2,S12meas2,S22meas2]=TRL('thru.txt','line.txt','reflect.txt',inputFile2);
        realS11meas2=real(S11meas2);
        imagS11meas2=imag(S11meas2);
        realS21meas2=real(S21meas2);
        imagS21meas2=imag(S21meas2);
        realS12meas2=real(S12meas2);
        imagS12meas2=imag(S12meas2);
        realS22meas2=real(S22meas2);
        imagS22meas2=imag(S22meas2);
    end
    
    % in case the data is outside of normal freq range, truncate it
    realS11meas2=realS11meas2(minidx:maxidx); 
    imagS11meas2=imagS11meas2(minidx:maxidx);
    realS21meas2=realS21meas2(minidx:maxidx);
    imagS21meas2=imagS21meas2(minidx:maxidx);
    realS12meas2=realS12meas2(minidx:maxidx);
    imagS12meas2=imagS12meas2(minidx:maxidx);
    realS22meas2=realS22meas2(minidx:maxidx);
    imagS22meas2=imagS22meas2(minidx:maxidx);
    
  % Combine real and imaginary parts, transpose to column vector
    S11meas2=realS11meas2+1j.*imagS11meas2;
    S11m2ds=interp1(wf,S11meas2,wfds);
    S11m2ds=S11m2ds(:); 
    S21meas2=realS21meas2+1j.*imagS21meas2;
    S21m2ds=interp1(wf,S21meas2,wfds);
    S21m2ds=S21m2ds(:); 
    S12meas2=realS12meas2+1j.*imagS12meas2;
    S12m2ds=interp1(wf,S12meas2,wfds);
    S12m2ds=S12m2ds(:); 
    S22meas2=realS22meas2+1j.*imagS22meas2;
    S22m2ds=interp1(wf,S22meas2,wfds);
    S22m2ds=S22m2ds(:); 
end


% % % % % % % % % % % % % % % % % % % %
% %  choose init guess method for LSQCurvefit
% % % % % % % % % % % % % % % % % % % %

if strcmp(init_method,'nrw')==1

  % NRW analysis
    [fnrw,realS11ex,imagS11ex,realS21ex,imagS21ex,...
        realS12ex,imagS12ex,realS22ex,imagS22ex]=fileformat4(nrwInput);
    wnrw=2*pi*fnrw;    
    [epsfwd,epsbw,mufwd,mubw]=waveguide_nrw(fnrw,realS11ex,imagS11ex,realS21ex,imagS21ex,...
        realS12ex,imagS12ex,realS22ex,imagS22ex,ls);
    epsNRW=epsbw;
    muNRW=mubw;
  % display the NRW results
    display('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')
    display('NRW processing finished')
    display(['Mean forward eps = ' num2str(mean(epsNRW))])
    display(['Mean forward mu = ' num2str(mean(muNRW))])
    display(['Mean backward eps = ' num2str(mean(epsbw))])
    display(['Mean backward mu = ' num2str(mean(mubw))])
    display('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')
    display(' ')
    
  % initial guesses based on NRW
    etguess=interp1(wnrw,epsNRW,wfds);
    ezguess=interp1(wnrw,epsNRW,wfds);
    mutguess=interp1(wnrw,muNRW,wfds);
    muzguess=interp1(wnrw,muNRW,wfds);
    etguess=etguess(:);
    ezguess=ezguess(:);
    mutguess=mutguess(:);
    muzguess=muzguess(:);

elseif strcmp(init_method,'nrwTruth')==1 % load initial values from a .mat file
    
    load('NRW_vals.mat');

  % initial guesses based on NRW
    wfds=linspace(wnrw(1),wnrw(end),numds);     
    etguess=interp1(wnrw,epstNRW,wfds);
    ezguess=interp1(wnrw,epszNRW,wfds);
    mutguess=interp1(wnrw,muNRW,wfds);
    muzguess=interp1(wnrw,muNRW,wfds);
    etguess=etguess(:);
    ezguess=ezguess(:);
    mutguess=mutguess(:);
    muzguess=muzguess(:);

elseif strcmp(init_method,'initup')==1 % one initial guess + update as we go

    etguess=zeros(length(wfds)+1,1);
    ezguess=zeros(length(wfds)+1,1);
    mutguess=zeros(length(wfds)+1,1);
    muzguess=zeros(length(wfds)+1,1);
    etguess(1)=(myret+1j*myiet);
    ezguess(1)=(myrez+1j*myiez);
    mutguess(1)=(myrmut+1j*myimut);
    muzguess(1)=(myrmuz+1j*myimuz);

elseif strcmp(init_method,'initone')==1 % same initial guess every time
    
    etguess=(myret+1j*myiet)*ones(length(wfds),1);
    ezguess=(myrez+1j*myiez)*ones(length(wfds),1);
    mutguess=(myrmut+1j*myimut)*ones(length(wfds),1);
    muzguess=(myrmuz+1j*myimuz)*ones(length(wfds),1);

end

%%  solution routine

% % % % % % % % % % % % % % % % % % % %
% % Calculation of A coefficients and scattering parameters
% % % % % % % % % % % % % % % % % % % %

% pre-allocate the solution vectors
etsol=zeros(length(wfds),1);
ezsol=zeros(length(wfds),1);
mutsol=zeros(length(wfds),1);
muzsol=zeros(length(wfds),1);
tsolve=zeros(length(wfds),1);
eigerror=zeros(length(wfds)+1,1);
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

for widx=1:length(wfds) % solve at each frequency
    tic;
    wval=wfds(widx);
    
    % sometimes NRW routine gives positive imag eps or mu...
    if imag(etguess(widx))>=0
        etguess(widx)=etguess(widx)';
        ezguess(widx)=ezguess(widx)';
    end
    if imag(mutguess(widx))>=0
        mutguess(widx)=mutguess(widx)';
        muzguess(widx)=muzguess(widx)';
    end
    
    if solveCase==7 % set the transverse values to known values
        ket=etguess(widx);
        kmut=mutguess(widx); 
    end

    display(sprintf('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%'))
    display(sprintf(['Now processing point ' num2str(widx) '/' num2str(length(wfds))]))
    display(sprintf(['Frequency = ' num2str(wval/2/pi) 'GHz']))
    display(sprintf(['Initial et =' num2str(etguess(widx)) ]))
    display(sprintf(['Initial ez =' num2str(ezguess(widx)) ]))
    display(sprintf(['Initial mut =' num2str(mutguess(widx)) ]))
    display(sprintf(['Initial muz =' num2str(muzguess(widx)) ]))
    display(sprintf('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%'))
    display(sprintf(' \n'))
    
  % construct the measured s-parameter matrices and downsample
    if porttouse==1 % Use only S11 and S21
        Smeas=cat(1,real(S11m1ds(widx)),imag(S11m1ds(widx)),real(S21m1ds(widx)),imag(S21m1ds(widx)));
    elseif porttouse==2 % use only S22 and S12
        Smeas=cat(1,real(S12m1ds(widx)),imag(S12m1ds(widx)),real(S22m1ds(widx)),imag(S22m1ds(widx)));
    elseif porttouse==3  % use both
        Smeas=cat(1,real(S11m1ds(widx)),imag(S11m1ds(widx)),real(S21m1ds(widx)),imag(S21m1ds(widx)),...
             real(S12m1ds(widx)),imag(S12m1ds(widx)),real(S22m1ds(widx)),imag(S22m1ds(widx)));
    end
    
    if (length(dmat)==2) && (porttouse==1)   % two thickness method
            Smeas=cat(1,Smeas,real(S11m2ds(widx)),imag(S11m2ds(widx)),real(S21m2ds(widx)),imag(S21m2ds(widx)));
    elseif (length(dmat)==2) && (porttouse==2) 
            Smeas=cat(1,Smeas,real(S12m2ds(widx)),imag(S12m2ds(widx)),real(S22m2ds(widx)),imag(S22m2ds(widx)));
    elseif (length(dmat)==2) && (porttouse==3)   % use both
            Smeas=cat(1,Smeas,real(S11m2ds(widx)),imag(S11m2ds(widx)),real(S21m2ds(widx)),imag(S21m2ds(widx)),...
              real(S12m2ds(widx)),imag(S12m2ds(widx)),real(S22m2ds(widx)),imag(S22m2ds(widx)));
    end

% % % % % % % % % % % % % % % % % % % % %  
% % Solve using lsqcurvefit - unknowns are [ret iet rez iez rmut imut rmuz imuz] 
% % % % % % % % % % % % % % % % % % % % %  
    
    [etsol(widx), ezsol(widx), mutsol(widx), muzsol(widx)] ...
        = runSolver(Smeas,wval,etguess(widx),ezguess(widx),mutguess(widx),muzguess(widx));
    
% % % % % % % % % % % % % % % % % % % % %  
% % error bars routine
% % % % % % % % % % % % % % % % % % % % %  

    if errchk==1
           [stddeltaetreal(widx), stddeltaetimag(widx), stddeltaezreal(widx), stddeltaezimag(widx),...
               stddeltamutreal(widx), stddeltamutimag(widx), stddeltamuzreal(widx), stddeltamuzimag(widx)]...
               = getErrorTerms(etsol(widx),ezsol(widx),mutsol(widx),muzsol(widx),Smeas,...
                wval,etguess(widx),ezguess(widx),mutguess(widx),muzguess(widx));

    end 
    
  % update initial guesses, if needed       
    if strcmp(init_method,'initup')==1
        etguess(widx+1)=etsol(widx);
        ezguess(widx+1)=ezsol(widx);
        mutguess(widx+1)=mutsol(widx);
        muzguess(widx+1)=muzsol(widx);
    end
    
    tsolve(widx)=toc;
    
  % output the final values
    display(sprintf('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%'))
    display(sprintf(['Final values for Frequency = ' num2str(wval/2/pi) 'GHz']))
    display(sprintf(['et = ' num2str(etsol(widx)) ]))
    display(sprintf(['ez = ' num2str(ezsol(widx)) ]))
    display(sprintf(['mut = ' num2str(mutsol(widx)) ]))
    display(sprintf(['muz = ' num2str(muzsol(widx)) ]))
    display(sprintf(['Time to solution = ' num2str(tsolve(widx)) 's']))
    display(sprintf('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%'))
    display(sprintf(' \n'));
end % end of w loop

%% final diagnostics
% % % % % % % % % % % % % % % % % % % % %  
% all done - print some diagnostics
% % % % % % % % % % % % % % % % % % % % %  

ttsolve=sum(tsolve);
hrs=floor(ttsolve/3600);
mins=floor((ttsolve-hrs*3600)/60);
secs=floor(ttsolve-hrs*3600-mins*60);
tstring=[num2str(hrs) ' hours, ' num2str(mins) ' mins, ' num2str(secs) ' secs' ];
meanet=mean(etsol);
meanez=mean(ezsol);
meanmut=mean(mutsol);
meanmuz=mean(muzsol);

display(sprintf(['Finished - total solution time = ' tstring]))
display(sprintf(['Avg et = ' num2str(meanet) ]))
display(sprintf(['Avg ez = ' num2str(meanez) ]))
display(sprintf(['Avg mut = ' num2str(meanmut) ]))
display(sprintf(['Avg muz = ' num2str(meanmuz) ]))


close all; % close the LSQ diagnostics

% if we want to smooth out the data
if smooth_zvals==1
    ezsolbak=ezsol;
    ezsol=smooth(ezsol,sspan,smeth);
end
if smooth_tvals==1
    etsolbak=etsol;
    etsol=smooth(etsol,sspan,smeth);
end

% save the workspace
save([casedesc '_workspace.mat']);

% % generate pretty pictures

% figure out the best limits, according to solution values
pad=0.2;

etminval=round(10*( (1+pad)*min(imag(etsol) ) ) )/10;
ezminval=round(10*( (1+pad)*min(imag(ezsol) ) ) )/10;
etmaxval=round(10*( (1+pad)*max(real(etsol) ) ) )/10;
ezmaxval=round(10*( (1+pad)*max(real(ezsol) ) ) )/10;
epslims=[min([etminval ezminval]) max([etmaxval ezmaxval])];

mutminval=round(10*( (1+pad)*min(imag(mutsol) ) ) )/10;
muzminval=round(10*( (1+pad)*min(imag(muzsol) ) ) )/10;
mutmaxval=round(10*( (1+pad)*max(real(mutsol) ) ) )/10;
muzmaxval=round(10*( (1+pad)*max(real(muzsol) ) ) )/10;
mulims=[min([mutminval muzminval]) max([mutmaxval muzmaxval])];

if errchk==1
    etminval=round(10*( (1+pad)*(epslims(1)-(2*min(stddeltaetimag ))) ) )/10;
    ezminval=round(10*( (1+pad)*(epslims(1)-(2*min(stddeltaezimag ))) ) )/10;
    etmaxval=round(10*( (1+pad)*(epslims(2)+(2*max(stddeltaetreal ))) ) )/10;
    ezmaxval=round(10*( (1+pad)*(epslims(2)+(2*max(stddeltaezreal ))) ) )/10;
    epslims_err=[min([etminval ezminval]) max([etmaxval ezmaxval])];

    mutminval=round(10*( (1+pad)*(mulims(1)-(2*min(stddeltamutimag ))) ) )/10;
    muzminval=round(10*( (1+pad)*(mulims(1)-(2*min(stddeltamuzimag ))) ) )/10;
    mutmaxval=round(10*( (1+pad)*(mulims(2)+(2*max(stddeltamutreal ))) ) )/10;
    muzmaxval=round(10*( (1+pad)*(mulims(2)+(2*max(stddeltamuzreal ))) ) )/10;
    mulims_err=[min([mutminval muzminval]) max([mutmaxval muzmaxval])];
end

sidebyside_pos_1=[0.1344    0.4   0.3319    0.4290];
sidebyside_pos_2=[0.5267    0.4   0.3319    0.4290];
sidebyside_annotation=[ 0.3184    0.1894    0.1336    0.0833];
sidebyside_leg=[0.4991    0.1239    0.1461    0.1935];

if strcmp(init_method,'nrw')==1 || strcmp(init_method,'nrwTruth')==1
    myfigure;
    orient landscape;
    subplot('Position',[0.0999    0.5514    0.3347    0.3412]);
    plot(wfds/2/pi/1e9,real(etsol),'-o',wfds/2/pi/1e9,imag(etsol),'-o',...
        wfds/2/pi/1e9,real(ezsol),'-o',wfds/2/pi/1e9,imag(ezsol),'-o');
    title('Extracted \epsilon_r Values')
    leg1=legend('Re(\sigma_t)','Imag(\sigma_t)','Re(\sigma_z)','Imag(\sigma_z)');
    xlabel('f (GHz)');
    ylim(epslims);
    set(leg1,'Position',[0.4619    0.6417    0.0994    0.1835])

    subplot('Position',[0.5891    0.5514   0.3347    0.3412]);
    plot(wfds/2/pi/1e9,real(mutsol),'-o',wfds/2/pi/1e9,imag(mutsol),'-o',...
        wfds/2/pi/1e9,real(muzsol),'-o',wfds/2/pi/1e9,imag(muzsol),'-o');
    title('Extracted \mu_r Values')
    xlabel('f (GHz)');
    ylim(mulims);
    
    subplot('Position',[0.0988    0.0776    0.3347    0.3412]);
    plot(wfds/2/pi/1e9,real(etsol),'-o',wfds/2/pi/1e9,imag(etsol),'-o',...
        wnrw/2/pi/1e9,real(epstNRW),'-o',wnrw/2/pi/1e9,imag(epstNRW),'-o');
    title('Comparison of \epsilon_t Values')
    leg2=legend('Re - tFWGT','Im - tFWGT','Re - NRW','Im - NRW');
    xlabel('f (GHz)');
    ylim(epslims);
    set(leg2,'Position',[0.4409    0.2032    0.1250    0.1324])

    subplot('Position',[0.5914    0.0776    0.3347    0.3412]); 
    plot(wfds/2/pi/1e9,real(mutsol),'-o',wfds/2/pi/1e9,imag(mutsol),'-o',...
        wnrw/2/pi/1e9,real(muNRW),'-o',wnrw/2/pi/1e9,imag(muNRW),'-o');
    title('Comparison of \mu_t Values')
    ylim(mulims);
    xlabel('f (GHz)');

    mtit(gcf,regexprep(casedesc,'_',' '),'FontSize',14,'yoff',0.06,'xoff',-0.02);
    
elseif strcmp(init_method,'nrw')==0 && (solveCase==3 || solveCase==6) % dielectric, magnetic

    myfigure;
    orient landscape;
    
    subplot('Position',sidebyside_pos_1);
    plot(wfds/2/pi/1e9,real(etsol),'-o',wfds/2/pi/1e9,real(ezsol),'-o',...
        wfds/2/pi/1e9,imag(etsol),'-o',wfds/2/pi/1e9,imag(ezsol),'-o');
    title('Extracted \epsilon Values')
    leg1=legend('Real - \sigma_t','Real - \sigma_z','Imag - \sigma_t','Imag - \sigma_z');
    xlabel('f (GHz)');
    ylim(epslims);
    
    subplot('Position',sidebyside_pos_2);
    plot(wfds/2/pi/1e9,real(mutsol),'-o',wfds/2/pi/1e9,real(muzsol),'-o',...
        wfds/2/pi/1e9,imag(mutsol),'-o',wfds/2/pi/1e9,imag(muzsol),'-o');
    title('Extracted \mu Values')
    ylim(mulims);
    xlabel('f (GHz)');
    
    set(leg1,'Position',sidebyside_annotation)
    annotation('textbox',sidebyside_leg,'string',{['Avg \epsilon_t = ' ...
        num2str(meanet)], ['Avg \epsilon_z = ' num2str(meanez)],['Avg \epsilon_t = ' ...
        num2str(meanmut)], ['Avg \mu_z = ' num2str(meanmuz)]},'FitBoxToText',...
        'on','VerticalAlignment','middle','BackgroundColor','w');
    mtit(gcf,regexprep(casedesc,'_',' '),'FontSize',14,'yoff',0.1,'xoff',-0.02);
    if errchk==1
        % one figure for permittivity
        myfigure;
        orient landscape;
        
        subplot('Position',sidebyside_pos_1);
        hold on;
        errorbar(wfds/2/pi/1e9,real(etsol),2*stddeltaetreal,'b');
        errorbar(wfds/2/pi/1e9,imag(etsol),2*stddeltaetimag,'Color',[0 0.4 0]);       
        title('Extracted \epsilon_t Values')
        leg1=legend('Real','Imag ');
        xlabel('f (GHz)');
        ylim(epslims_err);
        hold off;
        
        subplot('Position',sidebyside_pos_2);
        hold on;
        errorbar(wfds/2/pi/1e9,real(ezsol),2*stddeltaezreal,'b');
        errorbar(wfds/2/pi/1e9,imag(ezsol),2*stddeltaezimag,'Color',[0 0.4 0]);
        title('Extracted \epsilon_z Values')
        ylim(epslims_err);
        xlabel('f (GHz)');
        hold off;
        
        set(leg1,'Position',sidebyside_leg)
        annotation('textbox',sidebyside_annotation,'string',{['Avg \epsilon_t = ' ...
            num2str(meanet)], ['Avg \epsilon_z = ' num2str(meanez)]},'FitBoxToText',...
            'on','VerticalAlignment','middle','BackgroundColor','w');
        mtit(gcf,[regexprep(casedesc,'_',' ') ' eps uncertainty'],'FontSize',14,'yoff',0.1,'xoff',-0.02);
        
        % another figure for permeability
        myfigure;
        orient landscape;
        
        subplot('Position',sidebyside_pos_1);
        hold on;
        errorbar(wfds/2/pi/1e9,real(mutsol),2*stddeltamutreal,'b');
        errorbar(wfds/2/pi/1e9,imag(mutsol),2*stddeltamutimag,'Color',[0 0.4 0]);       
        title('Extracted \mu_t Values')
        leg1=legend('Real','Imag ');
        xlabel('f (GHz)');
        ylim(mulims_err);
        hold off;
        
        subplot('Position',sidebyside_pos_2);
        hold on;
        errorbar(wfds/2/pi/1e9,real(muzsol),2*stddeltamuzreal,'b');
        errorbar(wfds/2/pi/1e9,imag(muzsol),2*stddeltamuzimag,'Color',[0 0.4 0]);
        title('Extracted \mu_z Values')
        ylim(mulims_err);
        xlabel('f (GHz)');
        hold off;
        
        set(leg1,'Position',sidebyside_leg)
        annotation('textbox',sidebyside_annotation,'string',{['Avg \mu_t = ' ...
            num2str(meanmut)], ['Avg \mu_z = ' num2str(meanmuz)]},'FitBoxToText',...
            'on','VerticalAlignment','middle','BackgroundColor','w');
        mtit(gcf,[regexprep(casedesc,'_',' ') ' mu uncertainty'],'FontSize',14,'yoff',0.1,'xoff',-0.02);
    end
    
elseif strcmp(init_method,'nrw')==0 && (solveCase==2 || solveCase==5) % non-dielectric, magnetic
    myfigure;
    orient landscape;
    
    subplot('Position',sidebyside_pos_1);
    plot(wfds/2/pi/1e9,real(mutsol),'-o',wfds/2/pi/1e9,imag(mutsol),'-o');
    title('Extracted \mu_t Values')
    leg1=legend('Real','Imag');
    xlabel('f (GHz)');
    ylim(mulims);
    

    subplot('Position',sidebyside_pos_2);
    plot(wfds/2/pi/1e9,real(muzsol),'-o',wfds/2/pi/1e9,imag(muzsol),'-o');
    title('Extracted \mu_z Values')
    ylim(mulims);
    xlabel('f (GHz)');
    
    annotation('textbox',sidebyside_annotation,'string',{['Avg \mu_t = ' ...
        num2str(meanmut)], ['Avg \mu_z = ' num2str(meanmuz)]},'FitBoxToText',...
        'on','VerticalAlignment','middle','BackgroundColor','w');
    set(leg1,'Position',sidebyside_leg)
    mtit(gcf,regexprep(casedesc,'_',' '),'FontSize',14,'yoff',0.1,'xoff',-0.02);
    if errchk==1
        myfigure;
        orient landscape;

        subplot(Position',sidebyside_pos_1);
        hold on;
        errorbar(wfds/2/pi/1e9,real(mutsol),2*stddeltamutreal,'b');
        errorbar(wfds/2/pi/1e9,imag(mutsol),2*stddeltamutimag,'Color',[0 0.4 0]);      
        title('Extracted \mu_t Values')
        leg1=legend('Real','Imag');
        xlabel('f (GHz)');
        ylim(mulims_err);
        hold off;
        
        subplot('Position',sidebyside_pos_2);
        hold on;
        errorbar(wfds/2/pi/1e9,real(muzsol),2*stddeltamuzreal,'b');
        errorbar(wfds/2/pi/1e9,imag(muzsol),2*stddeltamuzimag,'Color',[0 0.4 0]);
        title('Extracted \mu_z Values')
        ylim(mulims_err);
        xlabel('f (GHz)');
        hold off;
        
        annotation('textbox',sidebyside_annotation,'string',{['Avg \mu_t = ' ...
            num2str(meanmut)], ['Avg \mu_z = ' num2str(meanmuz)]},'FitBoxToText',...
            'on','VerticalAlignment','middle','BackgroundColor','w');
        set(leg1,'Position',sidebyside_leg)        
        mtit(gcf,[regexprep(casedesc,'_',' ') ' uncertainty'],'FontSize',14,'yoff',0.1,'xoff',-0.02);
    end
    
    
elseif strcmp(init_method,'nrw')==0 && (solveCase==1 || solveCase==4) % dielectric, non-magnetic
    myfigure;
    orient landscape;
    
    subplot('Position',sidebyside_pos_1);
    plot(wfds/2/pi/1e9,real(etsol),'-o',wfds/2/pi/1e9,imag(etsol),'-o');
    title('Extracted \epsilon_t Values')
    leg1=legend('Real','Imag');
    xlabel('f (GHz)');
    ylim(epslims);

    subplot('Position',sidebyside_pos_2);
    plot(wfds/2/pi/1e9,real(ezsol),'-o',wfds/2/pi/1e9,imag(ezsol),'-o');
    title('Extracted \epsilon_z Values')
    ylim(epslims);
    xlabel('f (GHz)');
    
    annotation('textbox',sidebyside_annotation,'string',{['Avg \epsilon_t = ' ...
      num2str(meanet)], ['Avg \epsilon_z = ' num2str(meanez)]},'FitBoxToText','on',...
      'VerticalAlignment','middle','BackgroundColor','w');
    set(leg1,'Position',sidebyside_leg)
    mtit(gcf,regexprep(casedesc,'_',' '),'FontSize',14,'yoff',-0.06,'xoff',-0.02);
    
    if errchk==1
        myfigure;
        orient landscape;

        subplot('Position',sidebyside_pos_1);
        hold on;
        errorbar(wfds/2/pi/1e9,real(etsol),2*stddeltaetreal,'b');
        errorbar(wfds/2/pi/1e9,imag(etsol),2*stddeltaetimag,'Color',[0 0.4 0]);
        title('Extracted \epsilon_t Values')
        leg1=legend('Real','Imag');
        xlabel('f (GHz)');
        ylim(epslims_err);
        hold off;
        
        subplot('Position',sidebyside_pos_2);
        hold on;
        errorbar(wfds/2/pi/1e9,real(ezsol),2*stddeltaezreal,'b');
        errorbar(wfds/2/pi/1e9,imag(ezsol),2*stddeltaezimag,'Color',[0 0.4 0]);
        title('Extracted \epsilon_z Values')
        ylim(epslims_err);
        xlabel('f (GHz)');
        hold off;
        
        annotation('textbox',sidebyside_annotation,'string',{['Avg \epsilon_t = ' ...
          num2str(meanet)], ['Avg \epsilon_z = ' num2str(meanez)]},'FitBoxToText','on',...
          'VerticalAlignment','middle','BackgroundColor','w');
        set(leg1,'Position',sidebyside_leg)
        mtit(gcf,[regexprep(casedesc,'_',' ') ' uncertainty'],'FontSize',14,'yoff',-0.06,'xoff',-0.02);        
    end
end

fds=wfds(:)/2/pi/1e9;

if print==1
    printfigs; 
end


if makeFile==1
    if ismember(solveCase,[1 3 4 6])==1 && errchk==0
       M=[fds, real(etsol)];
       csvwrite([material '_solveCase_' num2str(solveCase) '_tfwgt_' num2str(numModes) 'mode_et_real' add2casedesc '.dat'],M)
       M=[fds, imag(etsol)];
       csvwrite([material '_solveCase_' num2str(solveCase) '_tfwgt_' num2str(numModes) 'mode_et_imag' add2casedesc '.dat'],M)      
       M=[fds, real(ezsol)];
       csvwrite([material '_solveCase_' num2str(solveCase) '_tfwgt_' num2str(numModes) 'mode_ez_real' add2casedesc '.dat'],M)  
       M=[fds, imag(ezsol)];
       csvwrite([material '_solveCase_' num2str(solveCase) '_tfwgt_' num2str(numModes) 'mode_ez_imag' add2casedesc '.dat'],M)  
       filename=[material '_solveCase_' num2str(solveCase) '_tfwgt_' num2str(numModes) 'mode_eps_stats' add2casedesc '.dat'];
       fid=fopen(filename,'w');
       ftitle='et re avg, et im avg, ez re avg, ez im avg';
       fprintf(fid,'%s\n',ftitle);
       results=horzcat(mean(real(etsol)),mean(imag(etsol)),mean(real(ezsol)),mean(imag(ezsol)));
       dlmwrite(filename,results,'-append');
       fclose(fid);
    elseif ismember(solveCase,[1 3 4 6])==1 && errchk==1
       M=[fds, real(etsol), 2*stddeltaetreal];
       csvwrite([material '_solveCase_' num2str(solveCase) '_tfwgt_' num2str(numModes) 'mode_et_real' add2casedesc '.dat'],M)
       M=[fds, imag(etsol), 2*stddeltaetimag];
       csvwrite([material '_solveCase_' num2str(solveCase) '_tfwgt_' num2str(numModes) 'mode_et_imag' add2casedesc '.dat'],M)      
       M=[fds, real(ezsol), 2*stddeltaezreal];
       csvwrite([material '_solveCase_' num2str(solveCase) '_tfwgt_' num2str(numModes) 'mode_ez_real' add2casedesc '.dat'],M)  
       M=[fds, imag(ezsol), 2*stddeltaezimag];
       csvwrite([material '_solveCase_' num2str(solveCase) '_tfwgt_' num2str(numModes) 'mode_ez_imag' add2casedesc '.dat'],M)  
       filename=[material '_solveCase_' num2str(solveCase) '_tfwgt_' num2str(numModes) 'mode_eps_stats' add2casedesc '.dat'];
       fid=fopen(filename,'w');
       ftitle='et re avg, et im avg, ez re avg, ez im avg';
       fprintf(fid,'%s\n',ftitle);
       results=horzcat(mean(real(etsol)),mean(imag(etsol)),mean(real(ezsol)),mean(imag(ezsol)));
       dlmwrite(filename,results,'-append');
    end
    if ismember(solveCase,[2 3 5 6])==1 && errchk==0
       M=[fds, real(mutsol)];
       csvwrite([material '_solveCase_' num2str(solveCase) '_tfwgt_' num2str(numModes) 'mode_mut_real' add2casedesc '.dat'],M)
       M=[fds, imag(mutsol)];
       csvwrite([material '_solveCase_' num2str(solveCase) '_tfwgt_' num2str(numModes) 'mode_mut_imag' add2casedesc '.dat'],M)      
       M=[fds, real(muzsol)];
       csvwrite([material '_solveCase_' num2str(solveCase) '_tfwgt_' num2str(numModes) 'mode_muz_real' add2casedesc '.dat'],M)
       M=[fds, imag(muzsol)];
       csvwrite([material '_solveCase_' num2str(solveCase) '_tfwgt_' num2str(numModes) 'mode_muz_imag' add2casedesc '.dat'],M)  
       filename=[material '_solveCase_' num2str(solveCase) '_tfwgt_' num2str(numModes) 'mode_mu_stats' add2casedesc '.dat'];
       fid=fopen(filename,'w');
       ftitle='mut re avg, mut im avg, muz re avg, muz im avg';
       fprintf(fid,'%s\n',ftitle);
       results=horzcat(mean(real(mutsol)),mean(imag(mutsol)),mean(real(muzsol)),mean(imag(muzsol)));
       dlmwrite(filename,results,'-append');
       fclose(fid);
    elseif ismember(solveCase,[2 3 5 6])==1 && errchk==1
       M=[fds, real(mutsol), 2*stddeltamutreal];
       csvwrite([material '_solveCase_' num2str(solveCase) '_tfwgt_' num2str(numModes) 'mode_mut_real' add2casedesc '.dat'],M)
       M=[fds, imag(mutsol), 2*stddeltamutimag];
       csvwrite([material '_solveCase_' num2str(solveCase) '_tfwgt_' num2str(numModes) 'mode_mut_imag' add2casedesc '.dat'],M)      
       M=[fds, real(muzsol), 2*stddeltamuzreal];
       csvwrite([material '_solveCase_' num2str(solveCase) '_tfwgt_' num2str(numModes) 'mode_muz_real' add2casedesc '.dat'],M)
       M=[fds, imag(muzsol), 2*stddeltamuzimag];
       csvwrite([material '_solveCase_' num2str(solveCase) '_tfwgt_' num2str(numModes) 'mode_muz_imag' add2casedesc '.dat'],M)  
       filename=[material '_solveCase_' num2str(solveCase) '_tfwgt_' num2str(numModes) 'mode_mu_stats' add2casedesc '.dat'];
       fid=fopen(filename,'w');
       ftitle='mut re avg, mut im avg, muz re avg, muz im avg';
       fprintf(fid,'%s\n',ftitle);
       results=horzcat(mean(real(mutsol)),mean(imag(mutsol)),mean(real(muzsol)),mean(imag(muzsol)));
       dlmwrite(filename,results,'-append');
       fclose(fid);
    end 
end

% outdated way of outputting results
% if  makeFile==1
%     % make a file with the relevant mean statistics
%     filename='tFWGT_material_stats.csv';
%     fid=fopen(filename,'w');  % create a results file
%     ftitle='avg et, avg ez, avg mut, avg muz';
%     fprintf(fid,'%s\r',ftitle);
%     dlmwrite(filename,[meanet meanez meanmut meanmuz],'-append');
%     fclose(fid);
% 
%     % save the individual values
%     csvwrite('etVals.csv',etsol);
%     csvwrite('ezVals.csv',ezsol);
%     csvwrite('mutVals.csv',mutsol);
%     csvwrite('muzVals.csv',muzsol);
% end

diary off;

% % % % % % % % % % % % % % % % % % % % %  
% change log
% % % % % % % % % % % % % % % % % % % % %  
% 20140514
%     - changed the filenames for the *.dat output files to be more
%     consistent and usefully sortable
% 20140512
%     - fixed the way error bars plots are made 
%     - updated the filename generation for pgfplots files
% 20140424
%     - added stats file outuput (for pgfplots)
% 20140418
%     - added ability to use quad2d rather than CIT integrals for lamy integrals
% 20140516
%     - fixed a couple of errors in Sparams case 5 integrals (11 and 12)
%     - added logic to output *.dat files 
%     - limits for errorbar plots are now calculated separately from normal plots
%        (helps with large error bars)
% 20140415
%     - changed logic for errchk plots 
%     - made error bars 2*sigma
%
% 20140225:
%     - updated S-parameters import to account for any order in the 
%             cti file
%     - added print switch
%     - reorganized options to be a little more intuitive 
%     - started the change log
%     - added ylimits to make the plots on the same scale
%     - added the average dotted line
%     - added the numds=0 option for no downsampling (requires no knowledge
%     of number of input points)
% 20140107:
%     - finally fixed the errorbars code!!!!
%
% 20140814:
%   - changed from kz notation to gamma notation (trying to fix multimode)
