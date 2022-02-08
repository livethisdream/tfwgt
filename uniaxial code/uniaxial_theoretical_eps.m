% clear all;
close all;

casedesc='Theoretical Parameters for Lexan and Particle Board Layered Material';
inputFile='tile_nrw_v1.txt';
lsa=3.08/1000;  % HAVE TO CHANGE THIS FOR EVERY MATERIAL!
[f,realS11ex,imagS11ex,realS21ex,imagS21ex,...
    realS12ex,imagS12ex,realS22ex,imagS22ex]=fileformat4(inputFile);
[epsfwda,epsbwa,mufwda,mubwa]=waveguide_nrw(f,realS11ex,imagS11ex,realS21ex,imagS21ex,realS12ex,imagS12ex,realS22ex,imagS22ex,lsa);
epsnrwa=epsfwda;

inputFile='particle_nrw_v1.txt';
lsb=3.52/1000;  % HAVE TO CHANGE THIS FOR EVERY MATERIAL!
[f,realS11ex,imagS11ex,realS21ex,imagS21ex,...
    realS12ex,imagS12ex,realS22ex,imagS22ex]=fileformat4(inputFile);
[epsfwdb,epsbwb,mufwdb,mubwb]=waveguide_nrw(f,realS11ex,imagS11ex,realS21ex,imagS21ex,realS12ex,imagS12ex,realS22ex,imagS22ex,lsb);
epsnrwb=epsbwb;

% ezguess=( (1./epsfwdb)-((epsfwda-epsfwdb)./(epsfwda.*epsfwdb)).*(lsa/(lsa+lsb)) ).^(-1);
% etguess=epsfwdb + (epsfwda-epsfwdb).*(lsa/(lsa+lsb));
% ezguess=( (1./epsbwb)-((epsbwa-epsbwb)./(epsbwa.*epsbwb)).*(lsa/(lsa+lsb)) ).^(-1);
% etguess=epsbwb + (epsbwa-epsbwb).*(lsa/(lsa+lsb));
ezguess=( (1./epsnrwb)-((epsnrwa-epsnrwb)./(epsnrwa.*epsnrwb)).*(lsa/(lsa+lsb)) ).^(-1);
etguess=epsnrwb + (epsnrwa-epsnrwb).*(lsa/(lsa+lsb));

myfigure;
subplot(2,1,1);
plot(f/1e9,real(etguess),f/1e9,imag(etguess));
title('Transverse');
ylabel('\epsilon_t');
xlabel('Frequency (GHz)');
legend('real','imag');

subplot(2,1,2);
plot(f/1e9,real(ezguess),f/1e9,imag(ezguess));
title('Longitudinal');
ylabel('\epsilon_z');
xlabel('Frequency (GHz)');
legend('real','imag');

mtit(gcf,casedesc,'FontSize',14,'yoff',0.03,'xoff',-0.02);

% printfigs