% Test_LAP
%% problem
c=[-1;-1];
n=2;
x = sdpvar(n,1);
K=[-4*x(1)+12*x(2)>=-1;-12*x(1)-4*x(2)>=-13;x>=0;x<=1];
f=c'*x;
plot(K);
V = [0 0; 0 1; 0.25 0; 0.75 1; 1 0.25]; % each row is a vertex
ops=sdpsettings('solver','gurobi','verbose',0);
params.Method=-1;
params.OutputFlag = 0;

%%
% run this part to add LAP cut step by step
hold on;
obj=export(K); % get constraints
optimize(K,f,ops); % optimize linear relaxation
x0=double(x); % get vertex solution
scatter(x0(1),x0(2),'filled','b'); % draw vertex
if norm(x0-round(x0))==0 % check optimality
    fprintf('* Optimized! the optimal solution is:\n');
    disp(x0'); 
    return;
end
[A,b]=computeBaseIneq(x0,obj.A,obj.rhs);

[lhs, rhs] = GMIR(A,[],b,1) % create GMIR cut
K = K + [lhs*x<=rhs]; % add cut to K
plot(K); % plot K
