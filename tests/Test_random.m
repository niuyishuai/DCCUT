%%
parpool('local',30);
%%
%%%%%%%%%%%%%%%%%%%%%%%%
% 3th EXAMPLE
% Random MBLP problems
%%%%%%%%%%%%%%%%%%%%%%%%
clear MIP;
MIP = CreateMIP(50,50,300);

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

ops.method = 1;
ops.numOfDCA = 1; 
ops.numOfLAP = 30;
fprintf('Method = %d, LAP = %d, ParDCA = 1\n',ops.method,ops.numOfLAP);

[xopt,yopt,fopt,info] = DCCUT(MIP,ops);

%%
% test using gurobi
[xopt_g,fopt_g,exitflag_g]=gurobimilp(MIP);
