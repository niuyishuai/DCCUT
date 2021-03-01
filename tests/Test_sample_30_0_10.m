%%
parpool('local',30);
%%
%%%%%%%%%%%%%%%%%%%%%%%%
% SECOND EXAMPLE
% 10-30 pure
%%%%%%%%%%%%%%%%%%%%%%%%
clear MIP;
load('matlab_30.mat')
%%
% test using DCCUT
ops=DCCUT_options;
ops.t=1000;
ops.maxiterDCCUT=inf;
ops.maxiterDCA=300;
ops.maxtime=120;
ops.verbose=1; % 1 display details 0 silence
ops.plot=0; % 1 plot iterations, 0 no plot
ops.paralLAP=true; % parallel LAP

ops.method = 2;
ops.numOfDCA = 1; 
ops.numOfLAP = 30; %feature('numCores');
fprintf('Method = %d, LAP = %d, ParDCA = 1\n',ops.method,ops.numOfLAP);

[xopt,yopt,fopt,info] = DCCUT(MIP,ops);
