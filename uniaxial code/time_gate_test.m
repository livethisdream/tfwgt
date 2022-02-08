% time gate s-parameters

clear all;
clc;
close all;

[~,realS11meas1,imagS11meas1,realS21meas1,imagS21meas1,...
    realS12meas1,imagS12meas1,realS22meas1,imagS22meas1]=fileformat4('material_d1.txt');
[f,realS11meas2,imagS11meas2,realS21meas2,imagS21meas2,...
    realS12meas2,imagS12meas2,realS22meas2,imagS22meas2]=fileformat4('material_d2.txt');

S11m1=realS11meas1+1j*imagS11meas1;
S21m1=realS21meas1+1j*imagS21meas1;
S12m1=realS12meas1+1j*imagS12meas1;
S22m1=realS22meas1+1j*imagS22meas1;

S11m2=realS11meas2+1j*imagS11meas2;
S21m2=realS21meas2+1j*imagS21meas2;
S12m2=realS12meas2+1j*imagS12meas2;
S22m2=realS22meas2+1j*imagS22meas2;

% S11m1zp=padarray(S11m1,2^15);
x=1:length(S11m1);
S11m1zp=S11m1;
w=kaiser(length(S11m1zp));
S11m1td=ifftshift(ifft(S11m1zp.*w,2^(nextpow2(length(S11m1)))));
plot(abs(S11m1td));
% xlim([0 10^4]);
df=4.2e9/201;
dt=1/4.2e9;
gate=zeros(length(S11m1td),1);
gate(124:133)=1;
hold on; 
plot(gate);
gateft=fft(gate);
gateftwin=gate.*kaiser(length(gate));
gatewin=ifft(gateftwin);
S11m1tdgated=S11m1td.*gatewin;
S11m1fdgated=fft(S11m1tdgated,length(S11m1));
S11m1gated=S11m1fdgated./w;
figure
plot(abs(S11m1gated));