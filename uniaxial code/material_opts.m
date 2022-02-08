function [myret, myiet, myrmut, myimut, myrez, myiez, myrmuz, myimuz, dmat]=material_opts(material)

material_choices={'fgm125','white_stuff','rcard','tile','particle','lexan',...
    'white_particle','tile_particle','tile_particle_thick','tile_white_stuff',...
    'particle_lexan','rcard_sandwhich','honeycomb','fgm125_white_stuff','fgm40',...
    'white_honeycomb','ui_honeycomb','uniaxial_square','uniaxial_square_2_75',...
    'white_honeycomb_125','white_honeycomb_500','ui_honeycomb_0_400'};
midx=find(ismember(material_choices,material));
if isempty(midx)==1
    midx=0;
end
switch midx
    case 0
        dmat=input('Enter dmat (as a matrix with []):   ');
        myret=2.9; 
        myiet=-0.05;
        myrmut=1;  
        myimut=-0.05;
        myrez=myret;
        myiez=myiet;
        myrmuz=myrmut;
        myimuz=myimut;
    case 1
        dmat=[3.12 6.24]./1000;  % mm to m - for 2 layers of fgm125
        % good guesses for fgm125
        myret=7.4988;  
        myiet=-0.015129;
        myrmut=0.50384;  
        myimut=-0.92066;
        myrez=myret;
        myiez=myiet;
        myrmuz=myrmut;
        myimuz=myimut;
        ls=dmat(1); % for NRW analysis
    case 2
        dmat=[3.97 7.94]./1000;  % mm to m - for 2 layers of white stuff
        % good guess for white stuff
        myret=2.9; 
        myiet=-0.05;
        myrmut=1;  
        myimut=-0.05;
        myrez=myret;
        myiez=myiet;
        myrmuz=myrmut;
        myimuz=myimut;
    case 3
        dmat=[1.455 2.91]./1000; % for rcard
        % good guess for rcard
        myret=30;
        myiet=-30;
        myrmut=1;
        myimut=0;
        myrez=myret;
        myiez=myiet;
        myrmuz=myrmut;
        myimuz=myimut;
    case 4
        dmat=[3.08 6.16]./1000;  % for tile
        myret=4;
        myiet=-0.05;
        myrmut=1;
        myimut=0;
        myrez=myret;
        myiez=myiet;
        myrmuz=myrmut;
        myimuz=myimut;
    case 5
        dmat=[3.52 7.04]./1000;  % for particle board
        myret=4;
        myiet=-0.05;
        myrmut=1;
        myimut=0;
        myrez=myret;
        myiez=myiet;
        myrmuz=myrmut;
        myimuz=myimut;
    case 6
        dmat=[1.38 2.76]./1000;  % for lexan
        myret=2;
        myiet=-0.05;
        myrmut=1;
        myimut=0;
        myrez=myret;
        myiez=myiet;
        myrmuz=myrmut;
        myimuz=myimut;
    case 7
        dmat=18.5/1000;          % for white_particle sandwich
        myret=3;
        myiet=-0.05;
        myrmut=1;
        myimut=0;
        myrez=myret;
        myiez=myiet;
        myrmuz=myrmut;
        myimuz=myimut;
    case 8
        dmat=16.28/1000;         % for tile_particle
        myret=4;
        myiet=-0.05;
        myrmut=1;
        myimut=0;
        myrez=myret;
        myiez=myiet;
        myrmuz=myrmut;
        myimuz=myimut;
    case 9
        dmat=49.72/1000;         % for tile_particle_thick
        myret=4;
        myiet=-0.05;
        myrmut=1;
        myimut=0;
        myrez=myret;
        myiez=myiet;
        myrmuz=myrmut;
        myimuz=myimut;
    case 10
        dmat=17.9/1000;          % for tile_white_stuff
        myret=3;
        myiet=-0.05;
        myrmut=1;
        myimut=0;
        myrez=myret;
        myiez=myiet;
        myrmuz=myrmut;
        myimuz=myimut;
    case 11
        dmat=13.32/1000;         % for particle_lexan
        % guess for tile & particle & lexan
        myret=4;
        myiet=-0.05;
        myrmut=1;
        myimut=0;
        myrez=myret;
        myiez=myiet;
        myrmuz=myrmut;
        myimuz=myimut;
    case 12
        dmat=9.3550/1000;        % for rcard sandwich v2
        % guess for rcard sandwich
        myret=11.4;
        myiet=-10.76;
        myrmut=1;
        myimut=0;
        myrez=3.75;
        myiez=-0.12;
        myrez=myret;
        myiez=myiet;
        myrmuz=myrmut;
        myimuz=myimut;
    case 13
        dmat=0.5*2.54/100;       % for honeycomb from emmerson-cuming
        myret=4;
        myiet=-0.5;
        myrmut=1;
        myimut=0;
        myrez=myret;
        myiez=myiet;
        myrmuz=myrmut;
        myimuz=myimut;
    case 14
        dmat=[10.21 17.3]./1000; % for fgm125_white_stuff
        myret=7.4988;  
        myiet=-0.015129;
        myrmut=0.50384;  
        myimut=-0.92066;
        myrez=myret;
        myiez=myiet;
        myrmuz=myrmut;
        myimuz=myimut;
    case 15
        dmat=[1 2]./1000;        % for 2 layers of fgm40
        % guess for fgm40
        myret=40;
        myiet=-0.1;
        myrmut=0.5;
        myimut=-0.05;
        myrez=myret;
        myiez=myiet;
        myrmuz=myrmut;
        myimuz=myimut;
    case 16
        dmat=[6.47]/1000;        % for 1 layer of white honeycomb
        myret=2.9; 
        myiet=-0.05;
        myrmut=1;  
        myimut=-0.05;
        myrez=myret;
        myiez=myiet;
        myrmuz=myrmut;
        myimuz=myimut;
    case 17
        dmat=6.34/1000;          % for 1/4" layer of ui honeycomb
        myret=2.9; 
        myiet=-0.05;
        myrmut=1;  
        myimut=-0.05;
        myrez=myret;
        myiez=myiet;
        myrmuz=myrmut;
        myimuz=myimut;
    case 18
        dmat=6.34/1000;          % for 1 layer of uniaxial white squares (1mm cells)
        myret=2.9; 
        myiet=-0.05;
        myrmut=1;  
        myimut=-0.05;
        myrez=myret;
        myiez=myiet;
        myrmuz=myrmut;
        myimuz=myimut;
    case 19
        dmat=6.35/1000;          % for 1 layer of uniaxial white squares (2.75mm cells)
        myret=2.9; 
        myiet=-0.05;
        myrmut=1;  
        myimut=-0.05;
        myrez=myret;
        myiez=myiet;
        myrmuz=myrmut;
        myimuz=myimut;
    case 20
        dmat=[3.175]/1000;       % for 1/8" layer of white honeycomb
        myret=2.9; 
        myiet=-0.05;
        myrmut=1;  
        myimut=-0.05;
        myrez=myret;
        myiez=myiet;
        myrmuz=myrmut;
        myimuz=myimut;
    case 21
        dmat=[12.7]/1000;        % for 1/2" layer of white honeycomb
        myret=2.9; 
        myiet=-0.05;
        myrmut=1;  
        myimut=-0.05;
        myrez=myret;
        myiez=myiet;
        myrmuz=myrmut;
        myimuz=myimut;
    case 22
        dmat=0.4*2.54/100;          % for 0.4" layer of ui honeycomb
        myret=2.9; 
        myiet=-0.75;
        myrmut=1;  
        myimut=-0.05;
        myrez=myret;
        myiez=myiet;
        myrmuz=myrmut;
        myimuz=myimut;
end