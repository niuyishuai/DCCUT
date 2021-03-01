%%
parpool('local',30);
%%
%%%%%%%%%%%%%%%%%%%%%%%%%
% 1st Example 10-10 pure
%%%%%%%%%%%%%%%%%%%%%%%%%
clear MIP;
MIP.C1 = [-2.2696;-0.3942;3.6018;-0.9882;-3.6844;2.9946;-0.6465;-0.4548;3.3601;-2.1507];
MIP.C2 = [];
MIP.A=[-0.9760   -0.2597    0.3479   -0.2245   -1.3819   -0.1699    1.8862    0.7490   -1.0811   -1.6861
   -0.4179   -0.5955   -0.3371   -0.3062    0.8015   -1.4063    1.5717   -1.2475   -0.1122    1.4859
   -0.7621    1.3052    0.5319    0.9497    1.6959    1.1078    0.6512    1.4769   -0.1982   -0.8773
   -0.6498    0.0076   -0.3066    1.8852    1.1138   -0.7372    0.6152   -0.5799   -0.0836   -0.0105
    0.0810   -0.6236   -0.8002    0.8480    1.6002   -1.5630   -0.5826    0.9058   -1.2923   -1.5472
   -1.0805    1.5700    1.0696   -0.0523   -0.5775    1.9393    0.7844    1.3766   -0.4787    0.6557
    1.9829   -0.0374    1.3243    1.4588   -1.8527   -0.7158    0.0847    1.4510   -0.5313    0.2882
   -1.7562    1.8777    1.2762   -1.1542    1.4269    1.9436   -0.9917    1.8222   -0.2950   -1.4279
   -0.5849    0.8489   -0.9042    0.8495   -0.6398   -1.2533   -0.5929    1.5026   -1.0576   -1.6837
   -0.8693    1.5267    1.2275   -0.6247    0.9537   -0.8005   -1.2797   -0.3857    0.5152    0.0165];
MIP.B = [];
MIP.b = [0.5288
    0.7537
    0.7413
    0.1773
    0.3005
    0.2894
    0.7848
    0.2411
    0.5476
    0.2180];
MIP.fopt = 0; % best solution for closed gap
%%
% test using DCCUT
ops=DCCUT_options;
ops.t=500;
ops.maxiterDCCUT=inf;
ops.maxiterDCA=300;
ops.maxtime=120;
ops.verbose=1; % 1 display details 0 silence
ops.plot=0; % 1 plot iterations, 0 no plot
ops.paralLAP=true; % parallel LAP

ops.method = 1;
ops.numOfDCA = 1; 
ops.numOfLAP = 30; %feature('numCores');
fprintf('Method = %d, LAP = %d, ParDCA = 1\n',ops.method,ops.numOfLAP);

[xopt,yopt,fopt,info] = DCCUT(MIP,ops);

%%
% test using gurobi
[xopt_g,fopt_g,exitflag_g]=gurobimilp(MIP);

%%
% turning parameter numOfLAP
ops.paralLAP=true;
ops.verbose=0; % 1 display details 0 silence
ops.plot=0; 
ops.method=1;
ops.numOfDCA = 1;
range=1:20;
timelst=zeros(1,length(range));
for i=1:length(range)
    ops.numOfLAP = range(i);
    [xopt,yopt,fopt,info] = DCCUT(MIP,ops);
    fprintf('method = %d, LAP = %d, ParDCA = 1, time = %.3f\n',ops.method,ops.numOfLAP,info.time);
    timelst(i)=info.time;
end
figure;
plot(range,timelst,'r-o','LineWidth',1.5);
xlabel('numOfLAP');
ylabel('time (sec)');

%%
% turning parameter numOfDCA
ops.isparal=false;
ops.verbose=0; % 1 display details 0 silence
ops.plot=0; 
ops.method=1;
ops.numOfLAP = 10;
range=1:20;
timelst=zeros(1,length(range));
for i=1:length(range)
    ops.numOfDCA = range(i);
    [xopt,yopt,fopt,info] = DCCUT(MIP,ops);
    fprintf('method = %d, LAP = %d, ParDCA = %d, time = %.3f\n',ops.method,ops.numOfLAP,ops.numOfDCA,info.time);
    timelst(i)=info.time;
end
figure;
plot(range,timelst,'r-o','LineWidth',1.5);
xlabel('numOfDCA');
ylabel('time (sec)');
