myfigure;
orient landscape;

colors={'blue',[0 0.4 0],'red','black','magenta'};
lspec={'-','--'};

subplot(2,1,1);
hold on;
plot(wfds/2/pi/1e9,real(etsol_1mode),'Color',colors{1},'LineStyle',lspec{1})
plot(wfds/2/pi/1e9,real(etsol_2mode),'Color',colors{2},'LineStyle',lspec{1})
plot(wfds/2/pi/1e9,real(etsol_3mode),'Color',colors{3},'LineStyle',lspec{1})
plot(wnrw/2/pi/1e9,real(epstNRW),'Color',colors{4},'LineStyle',lspec{1})

plot(wfds/2/pi/1e9,imag(etsol_1mode),'Color',colors{1},'LineStyle',lspec{2})
plot(wfds/2/pi/1e9,imag(etsol_2mode),'Color',colors{2},'LineStyle',lspec{2})
plot(wfds/2/pi/1e9,imag(etsol_3mode),'Color',colors{3},'LineStyle',lspec{2})
plot(wnrw/2/pi/1e9,imag(epstNRW),'Color',colors{4},'LineStyle',lspec{2})

title('Transverse');
legend('Real - TE10',...
    'Real - TE10, TE30',...
    'Real - TE10, TE30, TE12, TM12',...
    'Real - Free Space',...
    'Imag - TE10',...
    'Imag - TE10, TE30',...
    'Imag - TE10, TE30, TE12, TM12',...
    'Imag - Free Space',...
    'Location','NorthEastOutside')

subplot(2,1,2);
hold on;
plot(wfds/2/pi/1e9,real(ezsol_1mode),'Color',colors{1},'LineStyle',lspec{1})
plot(wfds/2/pi/1e9,real(ezsol_2mode),'Color',colors{2},'LineStyle',lspec{1})
plot(wfds/2/pi/1e9,real(ezsol_3mode),'Color',colors{3},'LineStyle',lspec{1})
plot(wnrw/2/pi/1e9,real(epszNRW),'Color',colors{4},'LineStyle',lspec{1})

plot(wfds/2/pi/1e9,imag(ezsol_1mode),'Color',colors{1},'LineStyle',lspec{2})
plot(wfds/2/pi/1e9,imag(ezsol_2mode),'Color',colors{2},'LineStyle',lspec{2})
plot(wfds/2/pi/1e9,imag(ezsol_3mode),'Color',colors{3},'LineStyle',lspec{2})
plot(wnrw/2/pi/1e9,imag(epszNRW),'Color',colors{4},'LineStyle',lspec{2})


title('Longitudinal');
legend('Real - TE10',...
    'Real - TE10, TE30',...
    'Real - TE10, TE30, TE12, TM12',...
    'Real - Free Space',...
    'Imag - TE10',...
    'Imag - TE10, TE30',...
    'Imag - TE10, TE30, TE12, TM12',...
    'Imag - Free Space',...
    'Location','NorthEastOutside')
mtit(gcf,['UI Honeycomb (0.400) - Hyde method'],'FontSize',14,'yoff',0.03,'xoff',-0.02);

