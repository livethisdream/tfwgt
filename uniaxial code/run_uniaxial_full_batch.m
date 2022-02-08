% controlling script to run batch jobs for uniaxial_full

global widx m_mat n_mat a b dmat eps0 mu0 numModes includeModes q_mat solveCase porttouse alg ket kmut num_int; 

%%%%%% UI Honeycomb 0.400 %%%%%%%%
cd('ui_honeycomb_0_400_tfwgt_final')

%%%% v1 port 1
material='ui_honeycomb_0_400';
vers='v1'; % the version of the measurement
orientation='_O2'; % needs the preceeding '_'
porttouse=1; % 1=use S11 & S21, 2=use S22 & S12, 3=use all
includeModes=[1]; % dominant mode only
uniaxial_full
display(['Finished processing ' casedesc])

%%%% v1 port 1 - 2 modes (10,12)
includeModes=[1 3 4]; 
uniaxial_full
display(['Finished processing ' casedesc])

%%%% v1 port 1 - 3 modes (10,12,14)
includeModes=[1 3 4 14 15]; 
uniaxial_full
display(['Finished processing ' casedesc])

%%%%%% White Honeycomb %%%%%%%%%

cd('../white_honeycomb_tfwgt_final')

%%%% v3 port 3
material='white_honeycomb';
vers='v3'; % the version of the measurement
orientation='_O1'; % needs the preceeding '_'
porttouse=3; % 1=use S11 & S21, 2=use S22 & S12, 3=use all
includeModes=[1]; % dominant mode only
uniaxial_full
display(['Finished processing ' casedesc])

%%%% v3 port 3 - 2 modes (10,12)
includeModes=[1 3 4]; 
uniaxial_full
display(['Finished processing ' casedesc])

%%%% v3 port 3 - 3 modes (10,12,14)
includeModes=[1 3 4 14 15]; 
uniaxial_full
display(['Finished processing ' casedesc])