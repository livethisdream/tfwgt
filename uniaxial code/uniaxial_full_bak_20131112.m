%%% full mode uniaxial constitutive parameter extraction
%% configuration/constants

clear all;
clc;
close all;

global widx v wi a b dmat eps0 mu0 numModes includeModes solveCase; 

% % % % % % % % % % % % % % % % % % % %
% % constants & global parameters
% % % % % % % % % % % % % % % % % % % %

dmat=[3.12 6.24]./1000; % mm to m - for 2 layers of fgm125
% dmat=[3.97 7.94]./1000; % mm to m - for 2 layers of white stuff
% dmat=[1.455 2.91]./1000; % for rcard
% dmat=[3.08 6.16]./1000; % for tile
% dmat=[3.52 7.04]./1000; % for particle board
% dmat=[1.38 2.76]./1000; % for lexan
% dmat=18.5/1000; % for white_particle sandwich
% dmat=16.28/1000; % for tile_particle
% dmat=17.9/1000; % for tile_white_stuff
% dmat=13.32/1000; % for particle_lexan
% dmat=9.3550/1000; % for rcard sandwich v2
v=[1 3 1 1 5 3 3 5 5 7 7 7 9 1 1 3 3 9 9 5]; % the first 20 modes
wi=[0 0 2 2 0 2 2 2 2 0 2 2 0 4 4 4 4 2 2 4];
c=3.0e8;           
eps0=8.854e-12;
mu0=pi*4e-7;
a=0.9*2.54/100;  % inches to m
b=0.4*2.54/100;  % inches to m
wmin=8.2e9*2*pi;
wmax=12.4e9*2*pi;

% % % % % % % % % % % % % % % % % % % %
% % user flags & config options
% % % % % % % % % % % % % % % % % % % %
material='fgm125';
methd='tfwgt'; % tfwgt or nrw
vers='v4'; % needs the preceeding v
init_method='nrw'; % 'nrw','nrwTruth','initone', or 'initup'
nrwvers='v6'; % no preceeding v
add2casedesc = '_testing';  % a string to add to the case description
includeModes=1; % dominant mode only
% includeModes=[1 3 4 14 15]; % the top 5 modes - 99% of the solution
numds=25;  % downsample the data to this number of points
errchk=0;    % error bar routine? 1=yes, 0=no
makefile=0;  % output the results to a csv file? 1=yes, 0=no
solveCase=6; % this is a switch to test some different cases, where we force
             %  certain symmetries or conditions on eps & mu
             % 1 - isotropic, dielectric, non-magnetic (et=ez=er) & (mut=muz=mu0)
             % 2 - isotropic, non-dielectric, magnetic (et=ez=eps0) & (mut=muz=mur)
             % 3 - isotropic, dielectric, magnetic (et=ez=er) & (mut=muz=mur)
             % 4 - uniaxial, dielectric, non-magnetic (et,ez) & (mut=muz=mu0)
             % 5 - uniaxial, non-dielectric, magnetic (et=ez=eps0) & (mut,muz)
             % 6 - uniaxial, dielectric, magnetic (et,ez) & (mut,muz)
if length(dmat)==1
    inputFile=[material '_' methd '_' vers '.txt'];
else
    inputFile=[material '_' methd '_' vers '_d1.txt'];
    inputFile2=[material '_' methd '_' vers '_d2.txt'];
end

numModes=length(includeModes); % number of modes to include

% % % % % % % % % % % % % % % % % % % %
% % setup for NRW analysis
% % % % % % % % % % % % % % % % % % % %
if strcmp(init_method,'nrw')==1
    nrwInput=[material '_nrw_' nrwvers '.txt'];
    ls=dmat(1); % for NRW analysis
end

% % % % % % % % % % % % % % % % % % % %             
% % initialize initial guesses
% % % % % % % % % % % % % % % % % % % %

% % good guesses for fgm125
% myret=7.4988;  
% myiet=-0.015129;
% myrmut=0.50384;  
% myimut=-0.92066;

% % good guess for white stuff
% myret=2.9; 
% myiet=-0.05;
% myrmut=1;  
% myimut=-0.05;

% % good guess for rcard
% myret=30;
% myiet=-30;
% myrmut=1;
% myimut=0;

% % guess for rcard sandwich
% myret=11.4;
% myiet=-10.76;
% myrmut=1;
% myimut=0;
% myrez=3.75;
% myiez=-0.12;

% guess for tile & particle & lexan
myret=4;
myiet=-0.05;
myrmut=1;
myimut=0;

% % give trans and long the same initial guess 
myrez=myret;
myiez=myiet;
myrmuz=myrmut;
myimuz=myimut;

% names for output files
if numModes==1
    casedesc=[material '_' methd '_' vers '_' upper(init_method) 'Init_' num2str(numds) ' points_' num2str(numModes) ' mode_solveCase' num2str(solveCase) add2casedesc];
    diary([upper(init_method) 'Init_' num2str(numds) 'pts_' num2str(numModes) 'mode_solveCase_' num2str(solveCase) add2casedesc '_log.txt']);
else
    casedesc=[material '_' methd '_' vers '_' upper(init_method) 'Init_' num2str(numds) ' points_' num2str(numModes) ' modes_solveCase' num2str(solveCase) add2casedesc];
    diary([upper(init_method) 'Init_' num2str(numds) 'pts_' num2str(numModes) ' modes_solveCase_' num2str(solveCase) add2casedesc '_log.txt']);
end

diary on;

%% get input data & setup initial guesses

% % % % % % % % % % % % % % % % % % % %
% %  Grab the input data & format it correctly
% % % % % % % % % % % % % % % % % % % %

[f,realS11meas1,imagS11meas1,realS21meas1,imagS21meas1,...
    realS12meas1,imagS12meas1,realS22meas1,imagS22meas1]=fileformat4(inputFile);

wf=2*pi*f;  % for all freq's
minidx=find(wf>=wmin,1,'first');  % in case the data is outside of normal freq range
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

% downsample the measured data
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
S22m1ds=interp1(wf,S22meas1,wfds);  % downsample
S22m1ds=S22m1ds(:);  % i like column vectors

% if using the two-thickness method, get second set of measurements
if length(dmat)==2
    [f,realS11meas2,imagS11meas2,realS21meas2,imagS21meas2,...
    realS12meas2,imagS12meas2,realS22meas2,imagS22meas2]...
        =fileformat4(inputFile2); % grab data from input file
    realS11meas2=realS11meas2(minidx:maxidx); % in case the data is outside of normal freq range
    imagS11meas2=imagS11meas2(minidx:maxidx);
    realS21meas2=realS21meas2(minidx:maxidx);
    imagS21meas2=imagS21meas2(minidx:maxidx);
    realS12meas2=realS12meas2(minidx:maxidx);
    imagS12meas2=imagS12meas2(minidx:maxidx);
    realS22meas2=realS22meas2(minidx:maxidx);
    imagS22meas2=imagS22meas2(minidx:maxidx);
    S11meas2=realS11meas2+1j.*imagS11meas2;
    S11m2ds=interp1(wf,S11meas2,wfds);
    S11m2ds=S11m2ds(:); % make it a column vector
    S21meas2=realS21meas2+1j.*imagS21meas2;
    S21m2ds=interp1(wf,S21meas2,wfds);
    S21m2ds=S21m2ds(:); % make it a column vector
    S12meas2=realS12meas2+1j.*imagS12meas2;
    S12m2ds=interp1(wf,S12meas2,wfds);
    S12m2ds=S12m2ds(:); % make it a column vector
    S22meas2=realS22meas2+1j.*imagS22meas2;
    S22m2ds=interp1(wf,S22meas2,wfds);
    S22m2ds=S22m2ds(:); % make it a column vector
end


% % % % % % % % % % % % % % % % % % % %
% %  choose init guess method for LSQCurvefit
% % % % % % % % % % % % % % % % % % % %

if strcmp(init_method,'nrw')==1

  % NRW analysis
    [fnrw,realS11ex,imagS11ex,realS21ex,imagS21ex,...
        realS12ex,imagS12ex,realS22ex,imagS22ex]=fileformat4(nrwInput);
%    ls=3.12/1000; % in m
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
    etguess=interp1(wnrw,epsNRW,wfds);
    ezguess=interp1(wnrw,epsNRW,wfds);
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
Y=0;  % initialize this
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
    display(sprintf('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%'))
    display(sprintf(['Frequency = ' num2str(wval/2/pi) 'GHz']))
    display(sprintf(['Initial et =' num2str(etguess(widx)) ]));
    display(sprintf(['Initial ez =' num2str(ezguess(widx)) ]));
    display(sprintf(['Initial mut =' num2str(mutguess(widx)) ]));
    display(sprintf(['Initial muz =' num2str(muzguess(widx)) ]));
    display(sprintf('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%'));
    display(sprintf(' \n'));
    
  % construct the measured s-parameter matrices and downsample
    if length(dmat)==1 % one thickness method
        Smeas=cat(1,real(S11m1ds(widx)),imag(S11m1ds(widx)),real(S21m1ds(widx)),imag(S21m1ds(widx)),...
              real(S12m1ds(widx)),imag(S12m1ds(widx)),real(S22m1ds(widx)),imag(S22m1ds(widx)));
    else   % two thickness method
        Smeas=cat(1,real(S11m1ds(widx)),imag(S11m1ds(widx)),real(S21m1ds(widx)),imag(S21m1ds(widx)),...
              real(S12m1ds(widx)),imag(S12m1ds(widx)),real(S22m1ds(widx)),imag(S22m1ds(widx)),...
              real(S11m2ds(widx)),imag(S11m2ds(widx)),real(S21m2ds(widx)),imag(S21m2ds(widx)),...
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
            = getErrorTerms(etsol(widx),ezsol(widx),mutsol(widx),muzsol(widx),Smeas,wval,etguess(widx),ezguess(widx),mutguess(widx),muzguess(widx));
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
display(sprintf(['Finished - total solution time = ' tstring]))
display(sprintf(['Avg et = ' num2str(mean(etsol)) ]))
display(sprintf(['Avg ez = ' num2str(mean(ezsol)) ]))
display(sprintf(['Avg mut = ' num2str(mean(mutsol)) ]))
display(sprintf(['Avg muz = ' num2str(mean(muzsol)) ]))

if  makefile==1
    % make a file with the relevant mean statistics
    filename='tFWGT_material_stats.csv';
    fid=fopen(filename,'w');  % create a results file
    ftitle='avg et, avg ez, avg mut, avg muz';
    fprintf(fid,'%s\r',ftitle);
    dlmwrite(filename,[mean(etsol) mean(ezsol) mean(mutsol) mean(muzsol)],'-append');
    fclose(fid);

    % save the individual values
    csvwrite('etVals.csv',etsol);
    csvwrite('ezVals.csv',ezsol);
    csvwrite('mutVals.csv',mutsol);
    csvwrite('muzVals.csv',muzsol);
end

close all; % close the LSQ diagnostics

% generate pretty pictures
if strcmp(init_method,'nrw')==1 || strcmp(init_method,'nrwTruth')==1
    myfigure;
    orient landscape;
    subplot('Position',[ 0.1073    0.5894    0.3347    0.3412]);
    plot(wfds/2/pi/1e9,real(etsol),wfds/2/pi/1e9,imag(etsol),wfds/2/pi/1e9,real(ezsol),wfds/2/pi/1e9,imag(ezsol));
    title('Extracted \epsilon_r Values')
    leg1=legend('Re(\sigma_t)','Imag(\sigma_t)','Re(\sigma_z)','Imag(\sigma_z)');
    xlabel('f (GHz)');
    set(leg1,'Position',[0.4640    0.6774    0.0952    0.1835])

    subplot('Position',[0.5891    0.5838    0.3347    0.3412]);
    plot(wfds/2/pi/1e9,real(mutsol),wfds/2/pi/1e9,imag(mutsol),wfds/2/pi/1e9,real(muzsol),wfds/2/pi/1e9,imag(muzsol));
    title('Extracted \mu_r Values')
    xlabel('f (GHz)');

    subplot('Position',[0.0988    0.1100    0.3347    0.3412]);
    plot(wfds/2/pi/1e9,real(etsol),wfds/2/pi/1e9,imag(etsol),wnrw/2/pi/1e9,real(epsNRW),wnrw/2/pi/1e9,imag(epsNRW));
    title('Comparison of \epsilon_t Values')
    leg2=legend('Re - tFWGT','Im - tFWGT','Re - NRW','Im - NRW');
    xlabel('f (GHz)');
    set(leg2,'Position',[0.4400    0.2389    0.1250    0.1324])

    subplot('Position',[0.5914    0.1100    0.3347    0.3412]); 
    plot(wfds/2/pi/1e9,real(mutsol),wfds/2/pi/1e9,imag(mutsol),wnrw/2/pi/1e9,real(muNRW),wnrw/2/pi/1e9,imag(muNRW));
    title('Comparison of \mu_t Values')
    xlabel('f (GHz)');

    mtit(gcf,regexprep(casedesc,'_',' '),'FontSize',14,'yoff',0.02,'xoff',-0.02);
    
elseif strcmp(init_method,'nrw')==0 && (solveCase==3 || solveCase==6) % dielectric, magnetic
    myfigure;
    orient landscape;
    subplot('Position',[ 0.1073    0.5894    0.3347    0.3412]);
    plot(wfds/2/pi/1e9,real(etsol),wfds/2/pi/1e9,imag(etsol));
    title('Extracted \epsilon_t Values')
    leg1=legend('Real','Imag');
    xlabel('f (GHz)');
    set(leg1,'Position',[0.4640    0.6774    0.0952    0.1835])

    subplot('Position',[0.5891    0.5838    0.3347    0.3412]);
    plot(wfds/2/pi/1e9,real(mutsol),wfds/2/pi/1e9,imag(mutsol));
    title('Extracted \mu_t Values')
    xlabel('f (GHz)');

    subplot('Position',[0.0988    0.1100    0.3347    0.3412]);
    plot(wfds/2/pi/1e9,real(ezsol),wfds/2/pi/1e9,imag(ezsol));
    title('Extracted \epsilon_z Values')
    leg2=legend('Real','Imag');
    xlabel('f (GHz)');
    set(leg2,'Position',[0.4400    0.2389    0.1250    0.1324])

    subplot('Position',[0.5914    0.1100    0.3347    0.3412]); 
    plot(wfds/2/pi/1e9,real(muzsol),wfds/2/pi/1e9,imag(muzsol));
    title('Extracted \mu_z Values')
    xlabel('f (GHz)');

    mtit(gcf,regexprep(casedesc,'_',' '),'FontSize',14,'yoff',0.02,'xoff',-0.02);
    
elseif strcmp(init_method,'nrw')==0 && (solveCase==2 || solveCase==5) % non-dielectric, magnetic
    myfigure;
    orient landscape;
    
    subplot('Position',[0.06    0.3776    0.4    0.5]);
    plot(wfds/2/pi/1e9,real(mutsol),wfds/2/pi/1e9,imag(mutsol));
    title('Extracted \mu_t Values')
    leg1=legend('Real','Imag');
    xlabel('f (GHz)');
    set(leg1,'Position',[0.4341    0.0905    0.1461    0.1935])

    subplot('Position',[0.5351    0.3776    0.4   0.5]);
    plot(wfds/2/pi/1e9,real(muzsol),wfds/2/pi/1e9,imag(muzsol));
    title('Extracted \mu_z Values')
    xlabel('f (GHz)');
    
    mtit(gcf,regexprep(casedesc,'_',' '),'FontSize',14,'yoff',0.1,'xoff',-0.02);
    
elseif strcmp(init_method,'nrw')==0 && (solveCase==1 || solveCase==4) % dielectric, non-magnetic
    myfigure;
    orient landscape;
    
    subplot('Position',[0.06    0.3776    0.4    0.5]);
    plot(wfds/2/pi/1e9,real(etsol),wfds/2/pi/1e9,imag(etsol));
    title('Extracted \epsilon_t Values')
    leg1=legend('Real','Imag');
    xlabel('f (GHz)');
    set(leg1,'Position',[0.4341    0.0905    0.1461    0.1935])

    subplot('Position',[0.5351    0.3776    0.4   0.5]);
    plot(wfds/2/pi/1e9,real(ezsol),wfds/2/pi/1e9,imag(ezsol));
    title('Extracted \epsilon_z Values')
    xlabel('f (GHz)');
    
    
    mtit(gcf,regexprep(casedesc,'_',' '),'FontSize',14,'yoff',0.1,'xoff',-0.02);
end

if errchk==1
    myfigure;
    orient landscape;

    subplot('Position',[ 0.1073    0.5894    0.3347    0.3412]);
    hold on;
    errorbar(wfds/2/pi/1e9,real(etsol),stddeltaetreal,'b');
    errorbar(wfds/2/pi/1e9,imag(etsol),stddeltaetimag,'g');
    title('Extracted \epsilon_t Values')
    xlabel('f (GHz)');
    leg1=legend('Re(\sigma_t)','Imag(\sigma_t)','Location','NorthEast');
    hold off;
    set(leg1,'Position',[0.4640    0.6774    0.0952    0.1835])

    subplot('Position',[0.0988    0.1100    0.3347    0.3412]);
    hold on;
    errorbar(wfds/2/pi/1e9,real(ezsol),stddeltaezreal,'r');
    errorbar(wfds/2/pi/1e9,imag(ezsol),stddeltaezimag,'c');
    title('Extracted \epsilon_z Values')
    xlabel('f (GHz)');
    hold off;

    subplot('Position',[0.5891    0.5838    0.3347    0.3412]);
    hold on;
    errorbar(wfds/2/pi/1e9,real(mutsol),stddeltamutreal,'b');
    errorbar(wfds/2/pi/1e9,imag(mutsol),stddeltamutimag,'g');
    title('Extracted \mu_t Values')
    xlabel('f (GHz)');
    hold off; 

    subplot('Position',[0.5914    0.1100    0.3347    0.3412]); 
    hold on;
    errorbar(wfds/2/pi/1e9,real(muzsol),stddeltamuzreal,'r');
    errorbar(wfds/2/pi/1e9,imag(muzsol),stddeltamuzimag,'c');
    title('Extracted \mu_z Values')
    leg2=legend('Re(\sigma_z)','Imag(\sigma_z)','Location','NorthEast');
    xlabel('f (GHz)');
    set(leg2,'Position',[0.4400    0.2389    0.1250    0.1324])
    hold off;
    
    
    mtit(gcf,regexprep(casedesc,'_',' '),'FontSize',14,'yoff',0.02,'xoff',-0.02);
end

printfigs; 

% save the workspace
if numModes==1
    save([upper(init_method) 'Init_' num2str(numds) 'pts_' num2str(numModes) 'mode_solveCase_' num2str(solveCase) '_workspace.mat']);
else
    save([upper(init_method) 'Init_' num2str(numds) 'pts_' num2str(numModes) 'modes_solveCase_' num2str(solveCase) '_workspace.mat']);
end

diary off;