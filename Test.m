%TEST NEW
fp=fopen('youyu1.txt','at+');
%parpool('local',2);
load('matlab_30.mat')

C1 = MIP.C1;
C2 = MIP.C2;
A = MIP.A;
B = MIP.B;
b = MIP.b;

intlinprog([C1;C2],1:30,A

ops.box_lb = zeros(30,1);
ops.box_ub = ones(30,1);

par.verbose=0; % 1 display details 0 silence
par.plot=1; % 1 plot iterations, 0 no plot
par.max_iter = 200;

par.method = 1;
par.numOfDCA = 1; %feature('numCores');
par.numOfLAP = 2;
[x1,y1,fval1,info1] = DC_CUT_gurobi(C1,C2,A,B,b,ops,par);

par.method = 2;
par.numOfDCA = 1; %feature('numCores');
par.numOfLAP = 2;
[x1,y1,fval1,info1] = DC_CUT_gurobi(C1,C2,A,B,b,ops,par);

par.method = 3;
par.t = 500;
par.numOfDCA = 1; %feature('numCores');
par.numOfLAP = 2;
[x1,y1,fval1,info1] = DC_CUT_gurobi(C1,C2,A,B,b,ops,par);

par.method = 4;
par.numOfDCA = 2; %feature('numCores');
par.numOfLAP = 2;
[x1,y1,fval1,info1] = DC_CUT_gurobi(C1,C2,A,B,b,ops,par);

par.method = 5;
par.numOfDCA = 2; %feature('numCores');
par.numOfLAP = 2;
[x1,y1,fval1,info1] = DC_CUT_gurobi(C1,C2,A,B,b,ops,par);

par.method = 6;
par.t = 500;
par.numOfDCA = 2; %feature('numCores');
par.numOfLAP = 2;
[x1,y1,fval1,info1] = DC_CUT_gurobi(C1,C2,A,B,b,ops,par);

%par.method = 2;
%par.numOfDCA = 2; %feature('numCores');
%par.numOfLAP = 2;
%[x2,y2,fval2,info2] = DC_CUT_gurobi(C1,C2,A,B,b,ops,par);
%fprintf(fp,'\n');
%if x1'*(1-x1)== 0
%    x = [x1;y1];
%    len = length(x);
%    for i = 1:len
%        fprintf(fp,'%2d \t',x(i));
%    end
%end
%fprintf(fp,'\n');
 
 %%

load('matlab_30.mat')

C1 = MIP.C1;
C2 = MIP.C2;
A = MIP.A;
B = MIP.B;
b = MIP.b;

model.modelsense = 'min';
model.obj = [C1;C2];
model.A = sparse([A,B]);
model.rhs = b;
model.lb = zeros(30,1);
model.ub = ones(30,1);
model.sense = repmat('<',size(A,1),1);
model.vtype = repmat('B',length(C1)+length(C2),1);

problem = gurobi(model);
%%
par.params.OutputFlag = 0;

A0 = [-2 -1;-2/3 1;1 0;0 1; -1 0;0 -1];
b0 = [-5/2;-1/6;0;0;-1;-1];
x0 = [3/4;1];
j=1;
[w,b] = LandP(A0,b0,x0,j);
%%

A0 = [-4 12;-12 -4;1 0;0 1; -1 0;0 -1];
b0 = [-1;-13;0;0;-1;-1];
x0 = [1;1/4];
j=2;
[w,b] = LandP(A0,b0,x0,j);