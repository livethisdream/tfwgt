function [stddeltaetreal,stddeltaetimag,stddeltaezreal,stddeltaezimag,...
    stddeltamutreal,stddeltamutimag,stddeltamuzreal,stddeltamuzimag]...
    = getErrorTerms(etsol,ezsol,mutsol,muzsol,...
    Smeas,wval,etguess,ezguess,mutguess,muzguess)

% % % % % % % % % % % % % % % % % % % % %  
% % error bars function
% % % % % % % % % % % % % % % % % % % % %  

global dmat solveCase;
   
% % % % % % % % % % % % % % % % % % % % %  
% % error due to variations in thicknesses d1 & d2
% % % % % % % % % % % % % % % % % % % % %

dmatbak=dmat; % backup old dmat values
vard=(5.0800e-05)^2;
deltaS=eps^(1/3);
dmat(1)=dmatbak(1)+deltaS;  % increment both thicknesses by just a tad

% solve for the sparams with a new d1 and original d2
[etDeltad1, ezDeltad1, mutDeltad1, muzDeltad1] = runSolver(Smeas,wval,etguess,ezguess,mutguess,muzguess);
% only if single thickness method:
etDeltad2 = 0;
ezDeltad2 = 0; 
mutDeltad2 = 0; 
muzDeltad2 = 0; 

if length(dmatbak)==2  % 2 thickness method solve with a new d2 & original d1
  dmat=[dmatbak(1) dmatbak(2)+deltaS];
  [etDeltad2,ezDeltad2,mutDeltad2,muzDeltad2] = runSolver(Smeas,wval,etguess,ezguess,mutguess,muzguess);
end
dmat=dmatbak;  % restore the old dmat value(s)

% % % % % % % % % % % % % % % % % % % % %  
% % error for variations in the S-parameters  
% % % % % % % % % % % % % % % % % % % % %  

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
  [etDeltaS,ezDeltaS,mutDeltaS,muzDeltaS] = runSolver(SmeasDelta(:,sidx),wval,etguess,ezguess,mutguess,muzguess);
% calculate the derivatives numerically 
  deltaSetreal(sidx)=(real(etDeltaS)-real(etsol))/deltaS; 
  deltaSetimag(sidx)=(imag(etDeltaS)-imag(etsol))/deltaS; 
  deltaSezreal(sidx)=(real(ezDeltaS)-real(ezsol))/deltaS; 
  deltaSezimag(sidx)=(imag(ezDeltaS)-imag(ezsol))/deltaS; 
  deltaSmutreal(sidx)=(real(mutDeltaS)-real(mutsol))/deltaS; 
  deltaSmutimag(sidx)=(imag(mutDeltaS)-imag(mutsol))/deltaS; 
  deltaSmuzreal(sidx)=(real(muzDeltaS)-real(muzsol))/deltaS; 
  deltaSmuzimag(sidx)=(imag(muzDeltaS)-imag(muzsol))/deltaS; 
end

% find the std deviations from the measured trans/refl coefficients
% Smeas is [real(S11m1) imag(S11m1) real(S21m1) imag(S21m1) 
%            real(S12m1) imag(S12m1) real(S22m1) imag(S22m1)
%            real(S11m2) imag(S11m2) real(S21m2) imag(S21m2) 
%            real(S12m2) imag(S12m2) real(S22m2) imag(S22m2)]
S11m1ds=Smeas(1)+1j*Smeas(2);
S21m1ds=Smeas(3)+1j*Smeas(4);
S12m1ds=Smeas(5)+1j*Smeas(6);
S22m1ds=Smeas(7)+1j*Smeas(8);
[u_21_m1,u_11_m1]=s_uncertainty(S21m1ds,S11m1ds); 
[u_12_m1,u_22_m1]=s_uncertainty(S12m1ds,S22m1ds);
sstdmat=[real(u_11_m1); imag(u_11_m1); real(u_21_m1); imag(u_21_m1); ...
    real(u_12_m1); imag(u_12_m1); real(u_22_m1); imag(u_22_m1)];
if length(dmat)==2
    S11m2ds=Smeas(9)+1j*Smeas(10);
    S21m2ds=Smeas(11)+1j*Smeas(12);
    S12m2ds=Smeas(13)+1j*Smeas(14);
    S22m2ds=Smeas(15)+1j*Smeas(16);
    [u_21_m2,u_11_m2]=s_uncertainty(S21m2ds,S11m2ds); 
    [u_12_m2,u_22_m2]=s_uncertainty(S12m2ds,S22m2ds);
    sstdmat=cat(1,sstdmat,[real(u_11_m2); imag(u_11_m2); real(u_21_m2); imag(u_21_m2); ...
    real(u_12_m2); imag(u_12_m2); real(u_22_m2); imag(u_22_m2)]);
end

% calculate the final standard deviation
stddeltaetreal=sqrt(sum((deltaSetreal.^2).*(sstdmat).^2)+(((real(etDeltad1)-real(etsol))/deltaS)^2)*(vard)+(((real(etDeltad2)-real(etsol))/deltaS)^2)*(vard));
stddeltaetimag=sqrt(sum((deltaSetimag.^2).*(sstdmat).^2)+(((imag(etDeltad1)-imag(etsol))/deltaS)^2)*(vard)+(((imag(etDeltad2)-imag(etsol))/deltaS)^2)*(vard));
stddeltaezreal=sqrt(sum((deltaSezreal.^2).*(sstdmat).^2)+(((real(ezDeltad1)-real(ezsol))/deltaS)^2)*(vard)+(((real(ezDeltad2)-real(ezsol))/deltaS)^2)*(vard));
stddeltaezimag=sqrt(sum((deltaSezimag.^2).*(sstdmat).^2)+(((imag(ezDeltad1)-imag(ezsol))/deltaS)^2)*(vard)+(((imag(ezDeltad2)-imag(ezsol))/deltaS)^2)*(vard));
stddeltamutreal=sqrt(sum((deltaSmutreal.^2).*(sstdmat).^2)+(((real(mutDeltad1)-real(mutsol))/deltaS)^2)*(vard)+(((real(mutDeltad2)-real(mutsol))/deltaS)^2)*(vard));
stddeltamutimag=sqrt(sum((deltaSmutimag.^2).*(sstdmat).^2)+(((imag(mutDeltad1)-imag(mutsol))/deltaS)^2)*(vard)+(((imag(mutDeltad2)-imag(mutsol))/deltaS)^2)*(vard));
stddeltamuzreal=sqrt(sum((deltaSmuzreal.^2).*(sstdmat).^2)+(((real(muzDeltad1)-real(muzsol))/deltaS)^2)*(vard)+(((real(muzDeltad2)-real(muzsol))/deltaS)^2)*(vard));
stddeltamuzimag=sqrt(sum((deltaSmuzimag.^2).*(sstdmat).^2)+(((imag(muzDeltad1)-imag(muzsol))/deltaS)^2)*(vard)+(((imag(muzDeltad2)-imag(muzsol))/deltaS)^2)*(vard));