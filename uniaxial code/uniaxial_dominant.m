%%% dominant mode uniaxial constitutive parameter extraction
%% configuration/constants

clear all;
clc;
close all;

% global myrez myiez myrmuz myimuz
global myway widx vm vn a b dmat m n eps0 mu0 deltamn Zntez Zmtez Zntmz Zmtmz ;
global Mhxmtez Mhxntez Mhymtez Mhyntez Mhxmtmz Mhxntmz Mhymtmz Mhyntmz;

% user flags
mailyes=0;
errchk=0;
myway=1;  % my way or the Hyde way
makefile=0;
numds=25;  % downsample the data to this number of points
casedesc='Init (bad guess), Full';
% % input file - the frequency ranges should be equal, or there's a problem
[~,realS11meas1,imagS11meas1,realS21meas1,imagS21meas1,...
    realS12meas1,imagS12meas1,realS22meas1,imagS22meas1]=fileformat4('fgm125_1_v4.txt');
[f,realS11meas2,imagS11meas2,realS21meas2,imagS21meas2,...
    realS12meas2,imagS12meas2,realS22meas2,imagS22meas2]=fileformat4('fgm125_2_v4.txt');

% define constants
dmat=[3.12 6.24]./1000; % mm to m - for 2 layers
m=1;
n=1;
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

% waveguide parameters
kxm=vm*pi/a;
kym=wm*pi/b;
kxn=vn*pi/a;
kyn=wn*pi/b;
wf=2*pi*f;  % for all freq's
numpts=length(f);  % need to know how many frequency points
wfds=linspace(wf(1),wf(end),numds);  % to downsample the data
k0=sqrt(wfds.^2.*eps0.*mu0);
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

Zmtez=wfds.*mu0./kzm;  % TEZ m index
Zntez=wfds.*mu0./kzn;  %  TEZ n index

Zmtmz=kzm./(wfds.*eps0);  % TMZ m index 
Zntmz=kzn./(wfds.*eps0);  % TMZ n index

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


% downsample the measured data
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

% a set of initial values - from Maj Hyde
myret=7.4988; 
myiet=-0.015129;
% myret=5;  % really bad guesses
% myiet=-0.1;
myrez=myret;
myiez=myiet;
myrmut=0.50384;
myimut=-0.92066;
% myrmut=1;  % really bad guesses
% myimut=-0.4;
myrmuz=myrmut;
myimuz=myimut;


% % just the first guess - update as we go
% etguess=zeros(length(wfds)+1,1);
% ezguess=zeros(length(wfds)+1,1);
% mutguess=zeros(length(wfds)+1,1);
% muzguess=zeros(length(wfds)+1,1);
% etguess(1)=(myret+1j*myiet);
% ezguess(1)=(myrez+1j*myiez);
% mutguess(1)=(myrmut+1j*myimut);
% muzguess(1)=(myrmuz+1j*myimuz);

% all one guess
etguess=(myret+1j*myiet)*ones(length(wfds),1);
ezguess=(myrez+1j*myiez)*ones(length(wfds),1);
mutguess=(myrmut+1j*myimut)*ones(length(wfds),1);
muzguess=(myrmuz+1j*myimuz)*ones(length(wfds),1);

% %  NRW data
[f,realS11ex,imagS11ex,realS21ex,imagS21ex,...
    realS12ex,imagS12ex,realS22ex,imagS22ex]=fileformat4('fgm125_nrw_6.txt');
% ls=3.12/1000; % in m
ls=2.99/1000;
[epsNRW,epsbw,muNRW,mubw]=waveguide_nrw(f,realS11ex,imagS11ex,realS21ex,imagS21ex,...
    realS12ex,imagS12ex,realS22ex,imagS22ex,ls);
% 
% % display the NRW results
% display('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')
% display('NRW processing finished')
% display(['Mean forward eps = ' num2str(mean(epsNRW))])
% display(['Mean forward mu = ' num2str(mean(muNRW))])
% display(' ')
% display(['Mean backward eps = ' num2str(mean(epsbw))])
% display(['Mean backward mu = ' num2str(mean(mubw))])
% display('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')
% 
% % initial guesses based on NRW
% etguess=interp1(wf,epsNRW,wfds);
% ezguess=interp1(wf,epsNRW,wfds);
% mutguess=interp1(wf,muNRW,wfds);
% muzguess=interp1(wf,muNRW,wfds);
% etguess=etguess(:);
% ezguess=ezguess(:);
% mutguess=mutguess(:);
% muzguess=muzguess(:);

% if using the two-thickness method
if length(dmat)==2
    S11meas2=realS11meas2+1j.*imagS11meas2;
    S11m2ds=interp1(wf,S11meas2,wfds);
    S11m2ds=S11m2ds(:); % make it a column vector
    S21meas2=realS21meas2+1j.*imagS21meas2;
    S21m2ds=interp1(wf,S21meas2,wfds);
    S21m2ds=S21m2ds(:);
    S12meas2=realS12meas2+1j.*imagS12meas2;
    S12m2ds=interp1(wf,S12meas2,wfds);
    S12m2ds=S12m2ds(:);
    S22meas2=realS22meas2+1j.*imagS22meas2;
    S22m2ds=interp1(wf,S22meas2,wfds);
    S22m2ds=S22m2ds(:); % make it a column vector
end


%% Calculation of A coefficients and scattering parameters

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

for widx=1:length(wfds)
    tic;
    wval=wfds(widx);
    display('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')
    display(['Frequency = ' num2str(wval/2/pi) 'GHz'])
    display(['Initial et =' num2str(etguess(widx))]);
    display(['Initial ez =' num2str(ezguess(widx))]);
    display(['Initial mut =' num2str(mutguess(widx))]);
    display(['Initial muz =' num2str(muzguess(widx))]);
    display('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%');
    display(' ');
    
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

% % % solve using lsqcurvefit
% unknowns are in the form [ret iet rez iez rmut imut rmuz imuz] - use
%  full NRW data for initial guesses or use one inital value and update as
%  we go

% % use the TRR, but a known ez and muz
%     options=optimset('Display','off','PlotFcns',{@optimplotfunccount @optimplotx @optimplotfval @optimplotstepsize});
%     lb=-10;
%     ub=10;
%     Y=lsqcurvefit(@Sparams,[real(etguess(widx)) imag(etguess(widx))...
%         real(mutguess(widx)) imag(mutguess(widx))],...
%         wval,Smeas,[0 lb 0 lb],[ub 0 ub 0],options);

% % use the TRR (trust-reflective-region algorithm - upper and lower bounds ensure a real material
        options=optimset('Display','off','PlotFcns',{@optimplotfunccount @optimplotx @optimplotfval @optimplotstepsize});
        lb=-10;
        ub=10; % restrict the lsq search
        Y=lsqcurvefit(@Sparams_dominant,[real(etguess(widx)) imag(etguess(widx)) real(ezguess(widx)) imag(ezguess(widx)) ...
            real(mutguess(widx)) imag(mutguess(widx)) real(muzguess(widx)) imag(muzguess(widx))],...
            wval,Smeas,[0 lb 0 lb 0 lb 0 lb],[ub 0 ub 0 ub 0 ub 0],options);

% % use LM algorithm (doesn't accept bounds)  
%     options=optimset('Display','off','PlotFcns',{@optimplotfunccount @optimplotx @optimplotfval @optimplotstepsize},...
%         'Algorithm',{'levenberg-marquardt',1e-6});
%     Y=lsqcurvefit(@Sparams,[real(etguess(widx)) imag(etguess(widx)) real(ezguess(widx)) imag(ezguess(widx)) ...
%             real(mutguess(widx)) imag(mutguess(widx)) real(muzguess(widx)) imag(muzguess(widx))],...
%             wval,Smeas,[],[],options);

% solve using nlinfit?

% save the final values
    etsol(widx)=Y(1)+1j*Y(2);
    ezsol(widx)=Y(3)+1j*Y(4);
    mutsol(widx)=Y(5)+1j*Y(6);
    muzsol(widx)=Y(7)+1j*Y(8);
    tsolve(widx)=toc;
    
% % do we want error bars?
        if errchk==1
           dmatbak=dmat; % backup old dmat values
           vard=(5.0800e-05)^2;
           deltaS=eps^(1/3);
           dmat(1)=dmatbak(1)+deltaS;  % increment the first thickness by just a tad
           
         % solve for the sparams with a new d1 and original d2
           Ydeltad1=lsqcurvefit(@Sparams,[real(etguess(widx)) imag(etguess(widx)) real(ezguess(widx)) imag(ezguess(widx)) ...
            real(mutguess(widx)) imag(mutguess(widx)) real(muzguess(widx)) imag(muzguess(widx))],...
            wval,Smeas,[0 lb 0 lb 0 lb 0 lb],[ub 0 ub 0 ub 0 ub 0],options);
           Ydeltad2=0;
           if length(dmatbak)==2  % solve for the sparams with a new d2 and original d1
              dmat=[dmatbak(1) dmatbak(2)+deltaS];  % increment second thickness by just a tad
              Ydeltad2=lsqcurvefit(@Sparams,[real(etguess(widx)) imag(etguess(widx)) real(ezguess(widx)) imag(ezguess(widx)) ...
                real(mutguess(widx)) imag(mutguess(widx)) real(muzguess(widx)) imag(muzguess(widx))],...
                wval,Smeas,[0 lb 0 lb 0 lb 0 lb],[ub 0 ub 0 ub 0 ub 0],options); 
           end
           dmat=dmatbak;  % restore the old dmat value
         
         % now, do the error for variations in the S-parameters  
           deltaSetreal=zeros(length(dmat)*8,1);
           deltaSetimag=zeros(length(dmat)*8,1);
           deltaSezreal=zeros(length(dmat)*8,1);
           deltaSezimag=zeros(length(dmat)*8,1);
           deltaSmutreal=zeros(length(dmat)*8,1);
           deltaSmutimag=zeros(length(dmat)*8,1);
           deltaSmuzreal=zeros(length(dmat)*8,1);
           deltaSmuzimag=zeros(length(dmat)*8,1);
           SmeasDelta=Smeas*ones(1,8*length(dmat))+deltaS*eye(8*length(dmat)); % increment the measured S by a small number
           for sidx=1:length((dmat))*8
              YdeltaS=lsqcurvefit(@Sparams,[real(etguess(widx)) imag(etguess(widx)) real(ezguess(widx)) imag(ezguess(widx)) ...
                real(mutguess(widx)) imag(mutguess(widx)) real(muzguess(widx)) imag(muzguess(widx))],...
                wval,SmeasDelta(:,sidx),[0 lb 0 lb 0 lb 0 lb],[ub 0 ub 0 ub 0 ub 0],options);
            % calculate the derivatives numerically 
              deltaSetreal(sidx)=(YdeltaS(1)-Y(1))/deltaS; 
              deltaSetimag(sidx)=(YdeltaS(2)-Y(2))/deltaS; 
              deltaSezreal(sidx)=(YdeltaS(3)-Y(3))/deltaS; 
              deltaSezimag(sidx)=(YdeltaS(4)-Y(4))/deltaS; 
              deltaSmutreal(sidx)=(YdeltaS(5)-Y(5))/deltaS; 
              deltaSmutimag(sidx)=(YdeltaS(6)-Y(6))/deltaS; 
              deltaSmuzreal(sidx)=(YdeltaS(7)-Y(7))/deltaS; 
              deltaSmuzimag(sidx)=(YdeltaS(8)-Y(8))/deltaS; 
           end
           svarmat=[var(real(S11meas1)); var(imag(S11meas1)); var(real(S21meas1)); var(imag(S21meas1));...
               var(real(S12meas1)); var(imag(S12meas1)); var(real(S22meas1)); var(imag(S22meas1))];
           if length(dmat)==2
           svarmat=cat(1,svarmat,[var(real(S11meas2)); var(imag(S11meas2)); var(real(S21meas2)); var(imag(S21meas2));...
               var(real(S12meas2)); var(imag(S12meas2)); var(real(S22meas2)); var(imag(S22meas2))]);
           end
           % calculate the final standard deviation
            stddeltaetreal(widx)=sqrt(sum((deltaSetreal.^2).*(svarmat))+( ( (Ydeltad1(1)-Y(1) ) / deltaS )^2)*(vard)+(((Ydeltad2(1)-Y(1))/deltaS)^2)*(vard));
            stddeltaetimag(widx)=sqrt(sum((deltaSetimag.^2).*(svarmat))+(((Ydeltad1(2)-Y(2))/deltaS)^2)*(vard)+(((Ydeltad2(2)-Y(2))/deltaS)^2)*(vard));
            stddeltaezreal(widx)=sqrt(sum((deltaSezreal.^2).*(svarmat))+(((Ydeltad1(3)-Y(3))/deltaS)^2)*(vard)+(((Ydeltad2(3)-Y(3))/deltaS)^2)*(vard));
            stddeltaezimag(widx)=sqrt(sum((deltaSezimag.^2).*(svarmat))+(((Ydeltad1(4)-Y(4))/deltaS)^2)*(vard)+(((Ydeltad2(4)-Y(4))/deltaS)^2)*(vard));
            stddeltamutreal(widx)=sqrt(sum((deltaSmutreal.^2).*(svarmat))+(((Ydeltad1(5)-Y(5))/deltaS)^2)*(vard)+(((Ydeltad2(5)-Y(5))/deltaS)^2)*(vard));
            stddeltamutimag(widx)=sqrt(sum((deltaSmutimag.^2).*(svarmat))+(((Ydeltad1(6)-Y(6))/deltaS)^2)*(vard)+(((Ydeltad2(6)-Y(6))/deltaS)^2)*(vard));
            stddeltamuzreal(widx)=sqrt(sum((deltaSmuzreal.^2).*(svarmat))+(((Ydeltad1(7)-Y(7))/deltaS)^2)*(vard)+(((Ydeltad2(7)-Y(7))/deltaS)^2)*(vard));
            stddeltamuzimag(widx)=sqrt(sum((deltaSmuzimag.^2).*(svarmat))+(((Ydeltad1(8)-Y(8))/deltaS)^2)*(vard)+(((Ydeltad2(8)-Y(8))/deltaS)^2)*(vard));
        end % end error bars routine
    
% % update initial guesses, if not using full NRW data
%     etguess(widx+1)=etsol(widx);
%     ezguess(widx+1)=ezsol(widx);
%     mutguess(widx+1)=mutsol(widx);
%     muzguess(widx+1)=muzsol(widx);

% output the final values
    display('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')
    display(['Final values for Frequency = ' num2str(wval/2/pi) 'GHz'])
    display(['et = ' num2str(etsol(widx))])
    display(['ez = ' num2str(ezsol(widx))])
    display(['mut = ' num2str(mutsol(widx))])
    display(['muz = ' num2str(muzsol(widx))])
    display(['Time to solution = ' num2str(tsolve(widx)) 's'])
    display('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')
    display(' ');
end % end of w loop

% all done - print some diagnostics
display(['Finished - total solution time = ' num2str(sum(tsolve)) 's'])
display(['Avg et = ' num2str(mean(etsol))])
display(['Avg ez = ' num2str(mean(ezsol))])
display(['Avg mut = ' num2str(mean(mutsol))])
display(['Avg muz = ' num2str(mean(muzsol))])

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
if errchk==0
    myfigure;
    orient landscape;
    subplot('Position',[ 0.1073    0.5894    0.3347    0.3412]);
    plot(wfds/2/pi,real(etsol),wfds/2/pi,imag(etsol),wfds/2/pi,real(ezsol),wfds/2/pi,imag(ezsol));
    title('Extracted \epsilon_r Values')
    leg1=legend('Re(\sigma_t)','Imag(\sigma_t)','Re(\sigma_z)','Imag(\sigma_z)');
    xlabel('f (GHz)');
    set(leg1,'Position',[0.4640    0.6774    0.0952    0.1835])
    
    subplot('Position',[0.5891    0.5838    0.3347    0.3412]);
    plot(wfds/2/pi,real(mutsol),wfds/2/pi,imag(mutsol),wfds/2/pi,real(muzsol),wfds/2/pi,imag(muzsol));
    title('Extracted \mu_r Values')
    xlabel('f (GHz)');
    
    subplot('Position',[0.0988    0.1100    0.3347    0.3412]);
    plot(wfds/2/pi,real(etsol),wfds/2/pi,imag(etsol),wf/2/pi,real(epsNRW),wf/2/pi,imag(epsNRW));
    title('Comparison of \epsilon_t Values')
    leg2=legend('Re - tFWGT','Im - tFWGT','Re - NRW','Im - NRW');
    xlabel('f (GHz)');
    set(leg2,'Position',[0.4400    0.2389    0.1250    0.1324])

    subplot('Position',[0.5914    0.1100    0.3347    0.3412]); 
    plot(wfds/2/pi,real(mutsol),wfds/2/pi,imag(mutsol),wf/2/pi,real(muNRW),wf/2/pi,imag(muNRW));
    title('Comparison of \mu_t Values')
    xlabel('f (GHz)');
    
    mtit(gcf,casedesc,'FontSize',14,'yoff',0.02,'xoff',-0.02);
else
    myfigure;
    orient landscape;
    subplot(2,2,1);
    hold on;
    errorbar(wfds/2/pi,real(etsol),stddeltaetreal,'b');
    errorbar(wfds/2/pi,imag(etsol),stddeltaetimag,'g');
    title('Extracted \epsilon_t Values')
    legend('Re(\epsilon_t)','Imag(\epsilon_t)','Location','NorthEastOutside');
    hold off;
    
    subplot(2,2,2);    
    hold on;
    errorbar(wfds/2/pi,real(ezsol),stddeltaezreal,'r');
    errorbar(wfds/2/pi,imag(ezsol),stddeltaezimag,'c');
    title('Extracted \epsilon_z Values')
    legend('Re(\epsilon_z)','Imag(\epsilon_z)','Location','NorthEastOutside');
    xlabel('f (GHz)');
    hold off;

    subplot(2,2,3);
    hold on;
    errorbar(wfds/2/pi,real(mutsol),stddeltamutreal,'b');
    errorbar(wfds/2/pi,imag(mutsol),stddeltamutimag,'g');
    title('Extracted \mu_t Values')
    legend('Re(\mu_t)','Imag(\mu_t)','Location','NorthEastOutside');
    hold off; 
    
    subplot(2,2,4);
    hold on;
    errorbar(wfds/2/pi,real(muzsol),stddeltamuzreal,'r');
    errorrbar(wfds/2/pi,imag(muzsol),stddeltamuzimag,'c');
    title('Extracted \mu_z Values')
    legend('Re(\mu_z)','Imag(\mu_z)','Location','NorthEastOutside');
    xlabel('f (GHz)');
    hold off;
end

printfigs; 

% send an email to tell me you're done - outlook
if mailyes==1
    sendolmail('livethisdream@gmail.com','Simulation Complete',['Extraction of constitutive parameters is complete' 10 ...
    'Parameters are as follows' 10 'Avg et = ' num2str(mean(etsol)) 10 'Avg ez = ' num2str(mean(ezsol)) 10 ...
    'Avg mut = ' num2str(mean(mutsol)) 10 'Avg muz = ' num2str(mean(muzsol))])
elseif mailyes==2  % gmail version
    maildone('Simulation Complete',['Extraction of constitutive parameters is complete' 10 ...
    'Parameters are as follows' 10 'Avg et = ' num2str(mean(etsol)) 10 'Avg ez = ' num2str(mean(ezsol)) 10 ...
    'Avg mut = ' num2str(mean(mutsol)) 10 'Avg muz = ' num2str(mean(muzsol))])
end
