% ************   DC cutting plane algorithm ************
% Problem (P):
%      min   C1*x+C2*y
%      s.t.  Ax+By <= b,
%            x \in {0,1}^n,
%      we assume that 0<=y<=ub is inserted in Ax+By <= b. 
%
% Enquivalent problem (DCt):
%     min    C1^T*x + C2^T*y + t*p(x,y)
%     s.t.   A*x + B*y <= b
%            x \in [0,1]^n
%
% where p(x) is a penalty function defined by
% p = sum(min(x,1-x))
%
% Syntax:
%   [x_opt,y_opt,fval,info] = DC_CUT_gurobi(C1,C2,A,B,beq,ops,[par])
% Usage: 
%        1. [x_opt,y_opt,fval] = DC_CUT_gurobi(C1,C2,A,B,beq,ops,[par])
%        2. [x_opt,y_opt,fval,info] = DC_CUT_gurobi(C1,C2,A,B,beq,ops,[par])
% INPUTS:
%    Mixed-Integer program (P)
%    C1: n-dimensional column vector
%    C2: q-dimensional column vector
%    A, B, b: matrices for defining linear constraints
%    
%    [ops.box_lb,ops.box_ub]: The domain for randomly choosing points
%    
%    par: a struct containing different values with respect to DCCUT algorithms
%       verbose: verbose =1 means showing the iteration process, =0 otherwise, [default = 1]
%       method:
%           method = 1 means adding LAP cuts when infeasible, and adding
%           dccut2 when feasible
%           method = 2 means when infeasible, adding dccut1 when well constructed, otherwise
%           adding LAP cuts, and when feasible, adding dccut2
%           method = 3 means when infeasible, adding LAP cuts, and adding
%           dccut1 if well constructed, when feasible, adding dccut2
%
%           method = 4 means the parallel version of method = 1
%           method = 5 means the parallel version of method = 2
%           method = 6 means the parallel version of method = 3
%
%       plot:    plot = 1 means plotting, =0 otherwise, [default = 0]
%       t:       t is the penalty parameter for the primal mixed binary problem, [default = 100]
%       max_iter: [default = 200]
%       numOfDCA: numOfDCA is the workers for parallelization [default = 1]
%       numOfLAP: numOfLAP is the maximal number for adding Lift-and-project cut at once, [default = 1]
%
% OUTPUTS:
%  *x_opt: optimal solution for integer variables.
%  *y_opt: optimal solution for continuous variables.
%  *fval: optimal value
%  *info: output informations
%    info.cuts: cuts informations
%    info.iter: number of iterations
%    info.nbdca: number of dca
%    info.time: cpu time
%    info.exitflag: -1: infeasible or unbounded, 1: optimized, 2: maxiter exceeded
% ========================================================
% The DC decompsition of f:
% f(u) = g(u) - h(u)
% g(u) = 0
% h(u) = -[C1^T*x + C2^T*y + t*p(x,y)]
% ========================================================
% written by Yi-Shuai Niu 2018/09/03
function [x_opt,y_opt,fval,info] = DC_CUT_gurobi(C1,C2,A,B,b,ops,par)
fp=fopen('youyu1.txt','at+');
if nargin < 6 || isempty(C1) || isempty(b) || isempty(ops)
    error ('Usage: [x_opt,y_opt,fval,info] = DC_CUT_gurobi(C1,C2,A,B,b,ops,par)');
end
if nargin < 7
    par = struct();
end

par = initialize_par_DCCUT(par);

max_iter = par.max_iter;
t = par.t;
numOfDCA = par.numOfDCA;
numOfLAP = par.numOfLAP;
eps = par.eps;
verbose = par.verbose;
method = par.method;
plot = par.plot;
params = par.params;

pval = @(x) sum(min(x,1-x));

% initialize outputs
x_opt = nan;
y_opt = nan;
fval = nan;

if ~isempty(nargout) || verbose
    info.cuts.dc1 = 0; % number of dc cuts (global cuts)
    info.cuts.dc2 = 0; % number of dc cuts (local cuts)
    info.cuts.lap1 = 0; % number of LAP cuts
    info.cuts.lap2 = 0; % number of LAP cuts
                             % number of iterations
    info.time = 0; % cpu time
    info.nbdca = 0; % number of dca
    info.exitflag = -1; % -1: infeasible or unbounded, 1: optimized, 2: maxiter exceeded
end

box_lb = ops.box_lb;
box_ub = ops.box_ub;

nrows = max(size(A,1),size(B,1));
% check and set right size for A and B
if isempty(A)
    A=zeros(nrows,0);
end
if isempty(B)
    B=zeros(nrows,0);
end
s1 = size(A,2);
s2 = size(B,2);
c=[C1;C2]; % linear objective function
A = [A;eye(s1);-eye(s1)];
B = [B;zeros(s1,s2);zeros(s1,s2)];
b = [b;ones(s1,1);zeros(s1,1)];

%if ~isempty(C2)
%    F = @ (x,y) c'*[x;y] + t*sum(min(x,1-x));
%else
%    F = @ (x) c'*x + t*sum(min(x,1-x));
%end

% call gurobi to solve this problem
model = setmodel(c,A,B,b);
problem = gurobi(model,params);
status = problem.status;

%[u,LB,flag] = linprog(c,[A,B],b,[],[],[],[],options);
%[u,LB,flag] = call_Gurobi(c,[A,B],b,[],[],[],[],[],options);

if ~strcmp(status,'OPTIMAL') % if the relaxation is infeasible or unbounded, then no solution found
    info.time=toc;
    return;
else
    u = problem.x;
    if  norm(round(u(1:s1))-u(1:s1),1) < 1e-8
        u = problem.x;
        x_opt = u(1:s1);         % judge whether u is feasible
        y_opt = u(s1+1:end);     % if yes, return u as the (global) minimizer
        fval = problem.objval;
        info.exitflag = 1;
        info.time = toc;
        return;
    end
end

% initialize upper bound (UB) as infinity
UB = inf;
LB = problem.objval;
if verbose == 1
    showcopyrights(fp);
    fprintf(fp,'%5s\t %8s\t %8s\t %8s\t %8s\t %8s\t %8s\t %8s\t %8s %8s\t\n','iter','ub','lb','cuts_dc1','cuts_dc2','cuts_lap1','cuts_lap2','nbdca','gap(%)','time');
end
if plot==1
    figure;
end

tic
info.iter = 0;

while (LB <= UB - eps && info.iter < max_iter)
    LB_old = LB;
    UB_old = UB;
    
    info.iter = info.iter + 1;
     
    if method == 1
        % ==================================================================
        % If (x,y) is infeasible, then add LAP cuts,
        % Otherwise if (x,y) is feasible, then
        %    add DC-CUT
        % ==================================================================
        
        % Use DCA to get a local minimizer
        % We set the relaxation solution u as the initial point of DCA
        [x,y] = my_DCA(C1,C2,t,A,B,b,u,params);
        info.nbdca=info.nbdca+1;
        
        % If the minimizer (x,y) is infeasible
        if norm(round(x)-x) > 1e-8
            % add a lift and project cut
            index = find(abs(x-round(x)) >=1e-8,numOfLAP);
            num = length(index);
            CUTS_A = cell(1,num);
            CUTS_B = cell(1,num);
            CUTS_b = cell(1,num);
            nb_lapcut=0;
                
            parfor i = 1:num
                j = index(i);
                [nu,zeta] = LandP(-[A,B],-b,[x;y],j,params);
                CUTS_A{i} = -nu(1:s1)';
                CUTS_B{i} = -nu(s1+1:end)';
                CUTS_b{i} = -zeta;   
            end
            info.cuts.lap1 = nb_lapcut+num;
            %Let us check all lift_and_project cuts returned by parallel procedure
            for i=1:num
                A=[A;CUTS_A{i}];
                B=[B;CUTS_B{i}];
                b=[b;CUTS_b{i}];
            end
        else    % If (x,y) is feasible
            % update the upper bound (UB) if the new feasible solution is better
            UB1 = c'*[x;y];
            if UB1 < UB
                UB = UB1;
                x_opt = x;
                y_opt = y;
                fval = UB;
            end
            
            % add DC cut
            J0 = find(x <= 1/2);
            J1 = find(x > 1/2);
            A1(J0) = -1;
            A1(J1) = 1;
            b1 = length(J1) -1;
            
            % update A,B,b
            A = [A;A1];
            B = [B;zeros(1,s2)];
            b = [b;b1];
            info.cuts.dc2=info.cuts.dc2+1;
            info.time = toc;
        end
             
    elseif method == 2
        % ==================================================================
        % If (x,y) is infeasible, then 
        %     add DC-CUT if p(x,y) is not integer, and all entries of x are not 1/2
        %     otherwise, add Lift-and-project cut
        % If (x,y) is feasible, then
        %    add DC-CUT
        % ==================================================================
        
        % Use DCA to get a local minimizer
        % We set the relaxation solution u as the initial point of DCA
        [x,y] = my_DCA(C1,C2,t,A,B,b,u,params);
        info.nbdca=info.nbdca+1;
        
        % If the minimizer (x,y) is infeasible
        if norm(round(x)-x) > 1e-8
            % Add DC-CUT if p(x) is not integer, and all entries of x are not 1/2
            pval = sum(min(x,1-x));
            if abs(pval-round(pval)) >= 1e-8 && norm(x -0.5,inf) > 1e-8
                out = check_construction(A,B,b,u,params);
                if out == 1
                    A1=zeros(1,length(x));
                    J0 = find(x <= 1/2);
                    J1 = find(x > 1/2);
                    A1(J0) = -1;
                    A1(J1) = 1;
                    b1 = length(J1) - ceil(sum(x(J0))+sum(1-x(J1)));
                
                    % Update A,B,b
                    A = [A;A1];
                    B = [B;zeros(1,s2)];
                    b = [b;b1];
                    info.cuts.dc1 = info.cuts.dc1+1;
                end
            else
                % otherwise, i.e., if p(u) is not integer or any entry of x is 1/2, 
                % add a lift and project cut
                index = find(abs(x-round(x)) >=1e-8,numOfLAP);
                %A = MIP.A;
                %B = MIP.B;
                %b = MIP.b;
                num = length(index);
                CUTS_A = cell(1,num);
                CUTS_B = cell(1,num);
                CUTS_b = cell(1,num);
                nb_lapcut=0;
                
                parfor i = 1:num
                    j = index(i);
                    [nu,zeta] = LandP(-[A,B],-b,[x;y],j,params);
                    CUTS_A{i} = -nu(1:s1)';
                    CUTS_B{i} = -nu(s1+1:end)';
                    CUTS_b{i} = -zeta;   
                end
                info.cuts.lap1 = nb_lapcut+num;
                %Let us check all lift_and_project cuts returned by parallel procedure
                for i=1:num
                    A=[A;CUTS_A{i}];
                    B=[B;CUTS_B{i}];
                    b=[b;CUTS_b{i}];
                end
            end
            
        else   % If (x,y) is feasible
            
            % update the upper bound (UB) if the new feasible solution is better
            UB1 = c'*[x;y];
            if UB1 < UB
                UB = UB1;
                x_opt = x;
                y_opt = y;
                fval = UB;
            end
            
            % add DC cut
            J0 = find(x <= 1/2);
            J1 = find(x > 1/2);
            A1(J0) = -1;
            A1(J1) = 1;
            b1 = length(J1) -1;
            
            % update A,B,b
            A = [A;A1];
            B = [B;zeros(1,s2)];
            b = [b;b1];
            info.cuts.dc2=info.cuts.dc2+1;
        end
        info.time = toc;
    
    elseif method == 3
        % ==================================================================
        % If (x,y) is infeasible, add LAP cuts
        %     then add DC CUT if p(x,y) is not integer, all entries of x
        %     are not 1/2, and is well constructed
        %     otherwise, 
        % If (x,y) is feasible, then
        %    add DC-CUT
        % ==================================================================
        
        % Use DCA to get a local solution
        % We set the relaxation solution u as the initial point of DCA
        [x,y] = my_DCA(C1,C2,t,A,B,b,u,params);
        info.nbdca=info.nbdca+1;
        
        % If the minimizer (x,y) is infeasible
        if norm(round(x)-x) > 1e-8
            % Add LAP cuts
            index = find(abs(x-round(x)) >=1e-8,numOfLAP);
            num = length(index);
            CUTS_A = cell(1,num);
            CUTS_B = cell(1,num);
            CUTS_b = cell(1,num);
            nb_lapcut=0;
                
            parfor i = 1:num
                j = index(i);
                [nu,zeta] = LandP(-[A,B],-b,[x;y],j,params);
                CUTS_A{i} = -nu(1:s1)';
                CUTS_B{i} = -nu(s1+1:end)';
                CUTS_b{i} = -zeta;   
            end
            info.cuts.lap1 = nb_lapcut+num;
            %Let us check all lift_and_project cuts returned by parallel procedure
            for i=1:num
                A=[A;CUTS_A{i}];
                B=[B;CUTS_B{i}];
                b=[b;CUTS_b{i}];
            end
            % Add DC-CUT if p(x) is not integer, all entries of x are not
            % 1/2, and is well constructed
            pval = sum(min(x,1-x));
            if abs(pval-round(pval)) >= 1e-8 && norm(x -0.5,inf) > 1e-8
                [x1,y1] = my_DCA(C1,C2,t,A,B,b,[x;y],params);
                
                
                [out,u1] = check_construction(A,B,b,[x;y],t,params)

                if c'*[x1;y1] + t*sum(min(x1,1-x1)) < c'*[x;y] + t*sum(min(x,1-x))    
                    A1=zeros(1,length(x));
                    J0 = find(x <= 1/2);
                    J1 = find(x > 1/2);
                    A1(J0) = -1;
                    A1(J1) = 1;
                    b1 = length(J1) - ceil(sum(x(J0))+sum(1-x(J1)));
                
                    % Update A,B,b
                    A = [A;A1];
                    B = [B;zeros(1,s2)];
                    b = [b;b1];
                    info.cuts.dc1 = info.cuts.dc1+1;
                end
            end
            
        else   % If (x,y) is feasible
            
            % update the upper bound (UB) if the new feasible solution is better
            UB1 = c'*[x;y];
            if UB1 < UB
                UB = UB1;
                x_opt = x;
                y_opt = y;
                fval = UB;
            end
            % add dccut2
            J0 = find(x <= 1/2);
            J1 = find(x > 1/2);
            A1(J0) = -1;
            A1(J1) = 1;
            b1 = length(J1) -1;
            
            % update A,B,b
            A = [A;A1];
            B = [B;zeros(1,s2)];
            b = [b;b1];
            info.cuts.dc2=info.cuts.dc2+1;
        end
        info.time = toc;
    
    elseif method==4
        % ==================================================================
        % This is the parallel version of method = 1
        % ==================================================================
            
        % initial number of cuts to be created parallelly
        CUTS_A=cell(1,numOfDCA);
        CUTS_B=cell(1,numOfDCA);
        CUTS_b=cell(1,numOfDCA);
       
        NEWFEASIBLES=cell(1,numOfDCA);
        NEWUBS=cell(1,numOfDCA);
        
        DCA_points = [u,MultiRandChoose(box_lb,box_ub,numOfDCA-1)];
        
        nb_dca = zeros(1,numOfDCA);
        %nb_dccut1 = zeros(1,numOfDCA);
        nb_dccut2 = zeros(1,numOfDCA);
        nb_lapcut1 = zeros(1,numOfDCA);
        
        % Use DCA to get a local solution;
        % We choose the relaxation solution u and some points in the box;
        parfor i = 1:numOfDCA
            [x,y] = my_DCA(C1,C2,t,A,B,b,DCA_points(:,i),params);
            nb_dca(i) = nb_dca(i) + 1;
            % If the minimizer (x,y) is infeasible
            if norm(round(x)-x,1) >= 1e-8 
                % add a LAP cut
                index = find(abs(x-round(x)) >=1e-4,numOfLAP);
                num = size(index,1);
                A_1=cell(1,num);
                B_1=cell(1,num);
                b_1=cell(1,num);
                    
                for k = 1:num
                    [nu,zeta] = LandP(-[A,B],-b,[x;y],index(k),params);
                    A_1{k} = -nu(1:s1)';
                    B_1{k} = -nu(s1+1:end)';
                    b_1{k} = -zeta; 
                end
                    
                CUTS_A{i} = [];
                CUTS_B{i} = [];
                CUTS_b{i} = [];
                    
                for k=1:num
                    CUTS_A{i} = [CUTS_A{i};A_1{k}];
                    CUTS_B{i} = [CUTS_B{i};B_1{k}];
                    CUTS_b{i} = [CUTS_b{i};b_1{k}];
                end
                nb_lapcut1(i) = nb_lapcut1(i)+num;
            
            else  %if (x,y) is feasible
                
                % update the upper bound (UB) if the new feasible solution is better
                UB1 = c'*[x;y];
                if UB1 < UB
                    NEWUBS{i}=UB1; % save better upper bounds
                    NEWFEASIBLES{i}=[x;y];
                end
                
                % add DC cut
                A1=zeros(1,length(x));
                J0 = x <= 1/2;
                J1 = find(x >  1/2);
                A1(J0) = -1;
                A1(J1) = 1;
                b1 = length(J1) -1;
                
                % Add dc cut
                CUTS_A{i}=A1;
                CUTS_B{i}=zeros(1,s2);
                CUTS_b{i}=b1;
                nb_dccut2(i)=nb_dccut2(i)+1;
            end
        end
        
        
        % Let us check all possible cuts returned by parallel procedure
        for i=1:numOfDCA
            % add cuts into A,B and b
            if ~isempty(CUTS_A{i})
                A=[A;CUTS_A{i}];
                B=[B;CUTS_B{i}];
                b=[b;CUTS_b{i}];
            end
            % update upper bounds
            if ~isempty(NEWUBS{i})
                if NEWUBS{i}<UB
                UB = NEWUBS{i};
                v = NEWFEASIBLES{i};
                x_opt = v(1:s1);
                y_opt = v(s1+1:end);
                fval = UB;
                end
            end
        end
        info.cuts.lap1 = info.cuts.lap1 + sum(nb_lapcut1);
        %info.cuts.dc1 = info.cuts.dc1 + sum(nb_dccut1);
        info.cuts.dc2 = info.cuts.dc2 + sum(nb_dccut2);
        info.nbdca = info.nbdca + sum(nb_dca);
        info.time = toc;
        
    elseif method==5
        % ==================================================================
        % This is the parallel version of method = 2
        % ==================================================================
            
        % initial number of cuts to be created parallelly
        CUTS_A=cell(1,numOfDCA);
        CUTS_B=cell(1,numOfDCA);
        CUTS_b=cell(1,numOfDCA);
        %nb_dca = 0;
        %nb_dccut1=0;
        %nb_dccut2=0;
        %nb_lapcut1=0;
        %nb_lapcut2=0;
        NEWFEASIBLES=cell(1,numOfDCA);
        NEWUBS=cell(1,numOfDCA);
        
        % Set the changed constrains in each iteration 
        %A=MIP.A;
        %B=MIP.B;
        %b=MIP.b;
        
        DCA_points = [u,MultiRandChoose(box_lb,box_ub,numOfDCA-1)];
        
        nb_dca = zeros(1,numOfDCA);
        nb_dccut1 = zeros(1,numOfDCA);
        nb_dccut2 = zeros(1,numOfDCA);
        nb_lapcut1 = zeros(1,numOfDCA);
        
        % Use DCA to get a local minimizer;
        % We choose the relaxation solution u and some points in the box;
        parfor i = 1:numOfDCA
            [x,y] = my_DCA(C1,C2,t,A,B,b,DCA_points(:,i),params);
            nb_dca(i) = nb_dca(i) + 1;
            % If the minimizer (x,y) is infeasible
            if norm(round(x)-x,1) >= 1e-8 
                % Add a DC-CUT if p(u) is not integer and all
                % entries of x are not 1/2
                pval=sum(min(x,1-x));
                if norm(x -0.5,inf) > 1e-8 && abs(pval-round(pval)) >= 1e-8
                    A1=zeros(1,length(x));
                    J0 = find(x <= 1/2);
                    J1 = find(x > 1/2);
                    A1(J0) = -1;
                    A1(J1) = 1;
                    b1 = length(J1) - ceil(sum(x(J0))+sum(1-x(J1)));
                    
                    % Add DC-CUT
                    CUTS_A{i}=A1;
                    CUTS_B{i}=zeros(1,s2);
                    CUTS_b{i}=b1;
                    nb_dccut1(i)=nb_dccut1(i)+1;
                
                else % add a Lift-and-project cut
                    index = find(abs(x-round(x)) >=1e-4,numOfLAP);
                    num = size(index,1);
                    A_1=cell(1,num);
                    B_1=cell(1,num);
                    b_1=cell(1,num);
                    
                    for k = 1:num
                        [nu,zeta] = LandP(-[A,B],-b,[x;y],index(k),params);
                        A_1{k} = -nu(1:s1)';
                        B_1{k} = -nu(s1+1:end)';
                        b_1{k} = -zeta; 
                    end
                    
                    CUTS_A{i} = [];
                    CUTS_B{i} = [];
                    CUTS_b{i} = [];
                    
                    for k=1:num
                        CUTS_A{i} = [CUTS_A{i};A_1{k}];
                        CUTS_B{i} = [CUTS_B{i};B_1{k}];
                        CUTS_b{i} = [CUTS_b{i};b_1{k}];
                    end
                    nb_lapcut1(i) = nb_lapcut1(i)+num;
                end
            else  %if (x,y) is feasible
                
                % update the upper bound (UB) if the new feasible solution is better
                UB1 = c'*[x;y];
                if UB1 < UB
                    NEWUBS{i}=UB1; % save better upper bounds
                    NEWFEASIBLES{i}=[x;y];
                end
                
                % add DC cut
                A1=zeros(1,length(x));
                J0 = x <= 1/2;
                J1 = find(x >  1/2);
                A1(J0) = -1;
                A1(J1) = 1;
                b1 = length(J1) -1;
                
                % Add dc cut
                CUTS_A{i}=A1;
                CUTS_B{i}=zeros(1,s2);
                CUTS_b{i}=b1;
                nb_dccut2(i)=nb_dccut2(i)+1;
            end
        end
        
        
        % Let us check all possible cuts returned by parallel procedure
        for i=1:numOfDCA
            % add cuts into A,B and b
            if ~isempty(CUTS_A{i})
                A=[A;CUTS_A{i}];
                B=[B;CUTS_B{i}];
                b=[b;CUTS_b{i}];
            end
            % update upper bounds
            if ~isempty(NEWUBS{i})
                if NEWUBS{i}<UB
                UB = NEWUBS{i};
                v = NEWFEASIBLES{i};
                x_opt = v(1:s1);
                y_opt = v(s1+1:end);
                fval = UB;
                end
            end
        end
        info.cuts.lap1 = info.cuts.lap1 + sum(nb_lapcut1);
        info.cuts.dc1 = info.cuts.dc1 + sum(nb_dccut1);
        info.cuts.dc2 = info.cuts.dc2 + sum(nb_dccut2);
        info.nbdca = info.nbdca + sum(nb_dca);
        info.time = toc;
        
    elseif method==6
        % ==================================================================
        % This is the parallel version of method = 3
        % ==================================================================
            
        % initial number of cuts to be created parallelly
        CUTS_A=cell(1,numOfDCA);
        CUTS_B=cell(1,numOfDCA);
        CUTS_b=cell(1,numOfDCA);
       
        NEWFEASIBLES=cell(1,numOfDCA);
        NEWUBS=cell(1,numOfDCA);
           
        DCA_points = [u,MultiRandChoose(box_lb,box_ub,numOfDCA-1)];
        
        nb_dca = zeros(1,numOfDCA);
        nb_dccut1 = zeros(1,numOfDCA);
        nb_dccut2 = zeros(1,numOfDCA);
        nb_lapcut1 = zeros(1,numOfDCA);
        
        % Use DCA to get a local minimizer;
        % We choose the relaxation solution u and some points in the box;
        parfor i = 1:numOfDCA
            [x,y] = my_DCA(C1,C2,t,A,B,b,DCA_points(:,i),params);
            nb_dca(i) = nb_dca(i) + 1;
            
            % If the minimizer (x,y) is infeasible
            if norm(round(x)-x,1) >= 1e-8 
                % Firstly add LAP cuts, and then add DC-CUT if p(u) is not integer, all
                % entries of x are not 1/2, and well constructed
                
                % add a Lift-and-project cut
                index = find(abs(x-round(x)) >=1e-4,numOfLAP);
                num = size(index,1);
                A_1=cell(1,num);
                B_1=cell(1,num);
                b_1=cell(1,num);
                    
                for k = 1:num
                    [nu,zeta] = LandP(-[A,B],-b,[x;y],index(k),params);
                    A_1{k} = -nu(1:s1)';
                    B_1{k} = -nu(s1+1:end)';
                    b_1{k} = -zeta; 
                end
                    
                CUTS_A{i} = [];
                CUTS_B{i} = [];
                CUTS_b{i} = [];
                    
                for k=1:num
                    CUTS_A{i} = [CUTS_A{i};A_1{k}];
                    CUTS_B{i} = [CUTS_B{i};B_1{k}];
                    CUTS_b{i} = [CUTS_b{i};b_1{k}];
                end
                nb_lapcut1(i) = nb_lapcut1(i)+num;
                
                % Add dccut2 if well constructed
                pval=sum(min(x,1-x));
                if norm(x -0.5,inf) > 1e-8 && abs(pval-round(pval)) >= 1e-8
                    A1=zeros(1,length(x));
                    J0 = find(x <= 1/2);
                    J1 = find(x > 1/2);
                    A1(J0) = -1;
                    A1(J1) = 1;
                    b1 = length(J1) - ceil(sum(x(J0))+sum(1-x(J1)));
                    
                    % Add DC-CUT
                    CUTS_A{i}=[CUTS_A{i};A1];
                    CUTS_B{i}=[CUTS_B{i};zeros(1,s2)];
                    CUTS_b{i}=[CUTS_b{i};b1];
                    nb_dccut1(i)=nb_dccut1(i)+1;
                end
                
            else  %if (x,y) is feasible
                
                % update the upper bound (UB) if the new feasible solution is better
                UB1 = c'*[x;y];
                if UB1 < UB
                    NEWUBS{i}=UB1; % save better upper bounds
                    NEWFEASIBLES{i}=[x;y];
                end
                
                % add DC cut
                A1=zeros(1,length(x));
                J0 = x <= 1/2;
                J1 = find(x >  1/2);
                A1(J0) = -1;
                A1(J1) = 1;
                b1 = length(J1) -1;
                
                % Add dc cut
                CUTS_A{i}=A1;
                CUTS_B{i}=zeros(1,s2);
                CUTS_b{i}=b1;
                nb_dccut2(i)=nb_dccut2(i)+1;
            end
        end
        
        
        % Let us check all possible cuts returned by parallel procedure
        for i=1:numOfDCA
            % add cuts into A,B and b
            if ~isempty(CUTS_A{i})
                A=[A;CUTS_A{i}];
                B=[B;CUTS_B{i}];
                b=[b;CUTS_b{i}];
            end
            % update upper bounds
            if ~isempty(NEWUBS{i})
                if NEWUBS{i}<UB
                UB = NEWUBS{i};
                v = NEWFEASIBLES{i};
                x_opt = v(1:s1);
                y_opt = v(s1+1:end);
                fval = UB;
                end
            end
        end
        info.cuts.lap1 = info.cuts.lap1 + sum(nb_lapcut1);
        info.cuts.dc1 = info.cuts.dc1 + sum(nb_dccut1);
        info.cuts.dc2 = info.cuts.dc2 + sum(nb_dccut2);
        info.nbdca = info.nbdca + sum(nb_dca);
        info.time = toc;
    end

    

    % This problem could be infeasible since the cutting planes are added
    % change -> [u,LB,flag] = call_Gurobi(c,[A,B],b,[],[],[],[],[],options);
    model = setmodel(c,A,B,b);
    %params.OutputFlag = 0;
    problem = gurobi(model,params);
    status = problem.status;
    
    %[u,LB,flag] = linprog(c,[A,B],b,[],[],[],[],options);
    % If the lower bound problem is infeasible, this occurs
    % when the problem is infeasible after added cuts or
    % when the problem has no integer solution
    if ~strcmp(status,'OPTIMAL') % return current best solution
        info.exitflag=1;
        info.time=toc;
        return;
    else  
        u = problem.x;
        LB = problem.objval;
        x = u(1:s1);
        y = u(s1+1:end);
        % when the problem has an infeasible solution,then one can introduce
        % the lift-and-project cut to cut off the point
        if norm(round(x)-x,1) >= 1e-8 
            index = find(abs(x-round(x)) >=1e-8,numOfLAP);
            
            num = size(index,1);
            CUTS_A = cell(1,num);
            CUTS_B = cell(1,num);
            CUTS_b = cell(1,num);
            nb_lapcut2=0;
                
           parfor i = 1:num
                j = index(i);
                [nu,zeta] = LandP(-[A,B],-b,[x;y],j,params);
                CUTS_A{i} = -nu(1:s1)';
                CUTS_B{i} = -nu(s1+1:end)';
                CUTS_b{i} = -zeta;
                
           end
           nb_lapcut2 = nb_lapcut2+num;
                
            % Let us check all LAP cuts returned by parallel procedure
            for i=1:num
                A=[A;CUTS_A{i}];
                B=[B;CUTS_B{i}];
                b=[b;CUTS_b{i}];
            end
            info.cuts.lap2 = info.cuts.lap2 + nb_lapcut2;
            
            % [u,LB,flag] = call_Gurobi(c,[A,B],b,[],[],[],[],[],options);
            model = setmodel(c,A,B,b);
            params.OutputFlag = 0;
            problem = gurobi(model,params);
            status = problem.status;
            % [u,LB,flag] = linprog(c,[A,B],b,[],[],[],[],options);
        end
        
        if ~strcmp(status,'OPTIMAL') % return current best solution
            info.exitflag=1;
            info.time=toc;
            return;
        elseif  norm(round(x)-x,1) < 1e-8 
            u = problem.x;
            LB = problem.objval;
            if LB < UB
                % if u is a better feasible solution
                % then Stop, update and return the best solution
                UB = LB;
                x_opt = u(1:s1);
                y_opt = u(s1+1:end);
                fval = UB;
                info.exitflag=1;
                info.time=toc;
                return;
            end
        end
    end
    
    
    if verbose==1
        fprintf(fp,'%5d\t %8.4f\t %8.4f\t %8d\t %8d\t %8d\t %8d\t %8d\t %8.4f\t %8.4f\n',info.iter,UB,LB,info.cuts.dc1,info.cuts.dc2,info.cuts.lap1,info.cuts.lap2,info.nbdca,100*(UB-LB)/(max(abs(UB),abs(LB))+1),info.time);
    end
    
    if plot==1
        myplot(LB_old,LB,UB_old,UB,info.iter);
    end
    
end
% check solution status
% if maxiter exceed, return the current best solution
if info.iter>= max_iter
    info.exitflag=2;
else % else optimal solution found
    info.exitflag=1;
end
info.time = toc;
fclose(fp);
end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% randomly choosing points function
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% a -> lower bound
% b -> upper bound
% n -> the number of chosen points
function y = MultiRandChoose(a,b,n)
s = length(a);
y = zeros(s,n);
for j = 1:n
    for i = 1:s
        y(i,j) = a(i) + (b(i)-a(i))*rand(1);    
    end    
end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Show version and copyrights
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function showcopyrights(fp)
fprintf(fp,'%s\n','*******************************************************************');
fprintf(fp,'%s\n', 'DC_CUT algorithm for solving Mixed Binary Linear Program (ver 1.0)');
fprintf(fp,'%s\n\n','*******************************************************************');
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Plot function
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function myplot(lb_old,lb,ub_old,ub,iter)
hold on;
if iter == 1
    %title('DC CUT ALGORITHM');
    xlabel('Iterations');
    ylabel('UB v.s. LB');
    if ~isinf(lb) && ~isnan(lb)
        plot(iter, lb,'k--s');
    end
    if ~isinf(ub) && ~isnan(ub)
        plot(iter, ub,'k-o');
    end
    drawnow
end
if iter>1
    if ~isinf(lb) && ~isnan(lb) && ~isinf(lb_old) && ~isnan(lb_old)
        plot([iter-1,iter], [lb_old,lb],'k--s');
    end
    if ~isinf(ub) && ~isnan(ub) && ~isinf(ub_old) && ~isnan(ub_old)
        plot([iter-1,iter], [ub_old,ub],'k-o');
    end
    if ~isinf(lb) && ~isnan(lb)
        plot(iter, lb,'k--s');
    end
    if ~isinf(ub) && ~isnan(ub)
        plot(iter, ub,'k-o');
    end
    drawnow
end
end