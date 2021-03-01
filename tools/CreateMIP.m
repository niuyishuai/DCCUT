function MIP = CreateMIP(n,m,q)
% MIP = CreateMIP(n,m,q)
% create MBLP problem
% min {C1'*x + C2'*y : A*x + B*y <= b; x binary, 0 <= y <= uby}
% where n is the number of binary variables, m is the number of continuous variables,
% q is the number of linear constraints. 
% For pure 0-1 problem, m=0.

% dimensions
MIP.n = n; % number of integer variables
MIP.m = m; % number of continuous variables
MIP.q = q; % number of linear constraints (A is a q x n matrix)

% coef for objective function
MIP.C1 = randi([-10,10],n,1);
MIP.C2 = 10*rand(m,1);

% linear constraints
MIP.A = sparse(randi([-10,10],q,n));
MIP.B = sparse(10*randn(q,m));
MIP.b = randi([0,10*(m+n)],q,1);

% bounds for binary and continuous variables
MIP.lbx = zeros(n,1);
MIP.ubx = ones(n,1);
MIP.lby = zeros(m,1);
MIP.uby = randi([1,5],m,1);
MIP.type='MILP';
MIPInfo(MIP);
end

function MIPInfo(MIP)
    fprintf('**********************************************************\n');
    fprintf('Model informations:\n');
    fprintf('Variables:(Binary variables: %d; Continuous variables: %d)\n',MIP.n, MIP.m);
    fprintf('Constraints:(Linear inequalities: %d)\n', length(MIP.b));
    fprintf('Problem type: MBLP\n');
    fprintf('**********************************************************\n');
end