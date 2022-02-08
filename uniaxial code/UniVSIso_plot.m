load('isoSol.mat')
myfigure; 
orient landscape;
plot(wfds/2/pi/1e9,real(etsol),wfds/2/pi/1e9,imag(etsol),wfds/2/pi/1e9,real(etisosol),wfds/2/pi/1e9,imag(etisosol))
legend('Re(\epsilon_t) - uniaxial','Im(\epsilon_t) - uniaxial','Re(\epsilon_t) - isotropic','Im(\epsilon_t) - isotropic')
title('Comparison of Uniaxial (1 Mode, SolveCase3) vs Isotropic (1 Mode) Permittivity - FGM125'); xlabel('Frequency (GHz)');
myfigure; 
orient landscape;
plot(wfds/2/pi/1e9,real(mutsol),wfds/2/pi/1e9,imag(mutsol),wfds/2/pi/1e9,real(mutisosol),wfds/2/pi/1e9,imag(mutisosol))
legend('Re(\mu_t) - uniaxial','Im(\mu_t) - uniaxial','Re(\mu_t) - isotropic','Im(\mu_t) - isotropic')
title('Comparison of Uniaxial (1 Mode, SolveCase3) vs Isotropic (1 Mode) Permeability - FGM125'); xlabel('Frequency (GHz)');
printfigs