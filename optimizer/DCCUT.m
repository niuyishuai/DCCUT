%DCCUT
%% DESCRIPTION:
% DC cutting plane algorithms for solving mixed-binary linear program:
%      min   f(x,y) = C1*x+C2*y
%      s.t.  A*x+B*y <= b,
%            x \in {0,1}^n,
%      we assume that 0<=y<=uby is included in Ax+By <= b.
%
% Enquivalent DC formulation:
%     min    f(x,y) + t*p(x)
%     s.t.   A*x + B*y <= b
%            x \in [0,1]^n
%
% where p(x) is a penalty function defined by p = sum(min(x,1-x))
% The DC decompsition of the objective function is given by:
% g(u) = 0
% h(u) = -[C1*x + C2*y + t*p(x)]
%
%% SYNTAX:
%   [xopt,yopt,fopt,info] = DCCUT(MIP,ops)
%
%% INPUTS:
%    -MIP: Mixed-binary linear program (MBLP)
%       MIP.C1: n-dimensional column vector
%       MIP.C2: q-dimensional column vector
%       MIP.A, MIP.B, MIP.b: matrices for linear constraints
%       MIP.uby: if presented, then incorporate them to A,B and b,
%                otherwise, uby is assumed to be included in A,B,b already.
%
%    -ops: options for DCCUT algorithms
%       verbose: 1 show iterations, 0 silence, [default = 1]
%       method: (will be updated later)
%           method = 1 adding dccut1 when when feasible, adding either dccut2 or LAP cut when infeasible.
%           method = 2 adding dccut1 when when feasible, adding both dccut2 and LAP when infeasible.
%           method = 3 adding dccut1 when when feasible, adding dccut2, LAP and GMIR when infeasible.
%           method = 5 adding LAP only
%           method = 11 parallel version of method 1
%           method = 12 parallel version of method 2
%       plot:    1 plotting, 0 otherwise, [default = 0]
%       t:       t is the initial penalty parameter [default = 100]
%       dca.maxiterDCA: max number of iterations for DCA [default = 200]
%       dca.tolx:    tolerence of variables for DCA [default = 1e-3]
%       dca.tolf:    tolerence of objective function for DCA [default = 1e-6]
%       tolgap:      tolerence of absolute gap for stopping DCCUT algorithm [default = 0.01]
%       zero:        option for floating zero [default = 1e-9]
%       numOfDCA: numOfDCA is the number of workers for parallel DCA in each iteration [default = 1]
%       numOfLAP: numOfLAP is the maximal number of LAP cuts at each iteration [default = 1]
%       paralLAP: true for creating LAP cuts in parallel [default = false]
%       params_gurobi: gurobi params, see gurobi
%       maxtime: max time for DCCUT [default = 3600]
%
%% OUTPUTS:
%  xopt: optimal solution for integer variables.
%  yopt: optimal solution for continuous variables.
%  fopt: optimal value
%  info: output information
%    info.cuts: cuts information
%    info.iter: number of iterations
%    info.nbdca: number of dca
%    info.time: computing time
%    info.exitflag: -1: infeasible or unbounded, 1: optimized, 2: maxiter exceeded
%%
%  See also DCCUT_options, gurobi
%
%% COPYRIGHT:
% Copyright since 2018, Yi-Shuai NIU. All Rights Reserved.
% 2018/09/03    initial code
% 2021/01/28    modified

function [xopt,yopt,fopt,info] = DCCUT(MIP,ops)
    if nargin < 1
        disp('see syntax by help DCCUT');
        return;
    end
    % initialize inputs
    C1=MIP.C1;
    C2=MIP.C2;
    A=MIP.A;
    B=MIP.B;
    b=MIP.b;
    bestknown = nan;
    if isfield(MIP,'fopt')
        bestknown = MIP.fopt; % best known optimal value if given
    end
    
    % initialize parameters
    if nargin <2
        ops = initialize_par_DCCUT();
    end
    maxiterDCCUT = ops.maxiterDCCUT;
    numOfDCA = ops.numOfDCA;
    t = ops.t; % penalty parameter
    tolgap = ops.tolgap; % gap tolerence
    method = ops.method; % algorithm
    plotflag = ops.plot;
    verboseflag = ops.verbose;
    params = ops.params_gurobi; % gurobi options
    opsdca = ops.dca; % dca options
    tolzero = ops.zero; % tolerence for zero
    maxtime = ops.maxtime; % max time for DCCUT
    
    % initialize outputs
    xopt = nan;
    yopt = nan;
    fopt = inf;
    info.cuts.dc1 = 0; % number of dc cuts type I (global cuts)
    info.cuts.dc2 = 0; % number of dc cuts type II (local cuts)
    info.cuts.lap = 0; % number of LAP cuts
    info.cuts.gom = 0;
    info.cuts.mir = 0; % number of MIR cuts
    info.iter = 0;
    info.time = 0; % cpu time
    info.nbdca = 0; % number of dca
    info.gap = nan;
    info.clgap = nan;
    info.exitflag = -1; % -1: infeasible or unbounded, 1: optimized, 2: maxiter exceeded
    
    % check and set A, B and b
    nrows = max(size(A,1),size(B,1));
    if isempty(A)
        A=zeros(nrows,0);
    end
    if isempty(B)
        B=zeros(nrows,0);
    end
    s1 = size(A,2); % number of integer variables
    s2 = size(B,2); % number of continuous variables
    c = [C1;C2]; % coefs of linear objective function
    % set lower and upper bounds of x by 0 and 1
    A = [A;eye(s1);-eye(s1)];
    B = [B;zeros(s1,s2);zeros(s1,s2)];
    b = [b;ones(s1,1);zeros(s1,1)];
    % set lower and upper bounds of y
    if s2 > 0
        if isfield(MIP,'uby')
            % set upper bound of y to uby and lower bound of y to 0
            A = [A;zeros(s2,s1);zeros(s2,s1)];
            B = [B;-eye(s2);eye(s2)];
            b = [b;zeros(s2,1);MIP.uby];
        else
            % set lower bound of y to 0
            A = [A;zeros(s2,s1)];
            B = [B;-eye(s2)];
            b = [b;zeros(s2,1)];
        end
    end
    
    % initialize objective upper bound as +inf and lower bound as -inf
    UB = inf;
    LB = -inf;
    if verboseflag == 1
        showcopyrights();
        fprintf('%5s | %10s | %10s | %8s | %8s | %8s | %8s | %8s | %8s | %9s | %8s\n','iter','ub','lb','cuts_dc1','cuts_dc2','cuts_lap','cuts_mir','nbdca','gap(%)','cl.gap(%)','time(s)');
        fprintf([repmat('-',1,122),'\n']);
    end
    if plotflag==1
        figure;
    end
    
    %%%%%% Start DCCUT algorithm %%%%%
    tic;
    while (UB - LB >= tolgap && info.iter < maxiterDCCUT && info.time < maxtime)
        LB_old = LB;
        UB_old = UB;
        
        % call gurobi to solve this problem
        model = setmodel(c,A,B,b);
        problem = gurobi(model,params);
        status = problem.status;
        if strcmp(status,'OPTIMAL')
            u0 = problem.x;
            x0 = u0(1:s1);
            y0 = u0(s1+1:end);
            if info.iter == 0
                initfval = problem.objval;
            end
        else
            % If the lower bound problem is infeasible, this occurs
            % when the original problem is infeasible
            % or the problem with cuts is infeasible.
            if info.iter == 0
                % if the relaxation is infeasible or unbounded, then no solution found
                info.time=toc;
                [info.gap,info.clgap]=updategap(UB,LB,initfval,bestknown);
                printfinalinfo(verboseflag,UB,LB,fopt,info);
                return;
            end
            % return the UB solution.
            info.time=toc;
            info.exitflag=1;
            LB = UB;
            [info.gap,info.clgap]=updategap(UB,LB,initfval,bestknown);
            printfverbose(verboseflag,UB,LB,info);
            plotiters(plotflag,LB_old,LB,UB_old,UB,info.iter);
            printfinalinfo(verboseflag,UB,LB,fopt,info);
            return;
        end
        
        % try to update lower bound
        if problem.objval > LB
            LB = problem.objval;
        end
        
        % check if the solution of the linear relaxation is feasible
        % u0 in S
        if  x0'*(1-x0) <= ops.zero
            % return better solution
            info.time=toc;
            fval = problem.objval;
            if fval < fopt
                % if the solution is better then update best solution
                xopt = x0;
                yopt = y0;
                fopt = fval;
                info.exitflag = 1;
                UB = fval;
            end
            [info.gap,info.clgap]=updategap(UB,LB,initfval,bestknown);
            printfverbose(verboseflag,UB,LB,info);
            plotiters(plotflag,LB_old,LB,UB_old,UB,info.iter);
            printfinalinfo(verboseflag,UB,LB,fopt,info);
            return;
        else
            %if rrr
            % add random pertubation on inital point u
            %u = u+randn(length(u),1);
            switch method
                case 1
                    % ==================================================================
                    % If (x,y) is feasible, then add a DC cut type I
                    % Otherwise if p(x) is not integer, and all entries of x are not 1/2
                    % then add a DC cut type II. Else, add other cuts
                    % ==================================================================
                    
                    % Use DCA to get a local minimizer
                    % (Linear relaxation for initial point of DCA)
                    [x,y] = DCA(model,t,s1,u0,params,opsdca);
                    info.nbdca=info.nbdca+1;
                    
                    % If (x,y) is feasible
                    if x'*(1-x) <= ops.zero
                        % add DC cut
                        [A1,B1,b1] = addDCcutI(x,s1,s2);
                        % update A,B,b
                        A = [A;A1];
                        B = [B;B1];
                        b = [b;b1];
                        info.cuts.dc1=info.cuts.dc1+1;
                        
                        % update the upper bound (UB) if the new feasible solution is better
                        fval = c'*[x;y];
                        if fval < UB
                            UB = fval;
                            xopt = x;
                            yopt = y;
                            fopt = UB;
                        end
                    else % If the (x,y) is infeasible
                        % use proc-V to get a critical point u_new and
                        % update t (not used yet)
                        pval = sum(min(x,1-x));
                        if abs(pval-round(pval)) >= ops.zero && norm(x-0.5,inf) > ops.zero
                            % add DC cut type II
                            [A1,B1,b1] = addDCcutII(x,s1,s2);
                            % update A,B,b
                            A = [A;A1];
                            B = [B;B1];
                            b = [b;b1];
                            info.cuts.dc2=info.cuts.dc2+1;
                        else
                            % add a lift and project cut
                            [A1,B1,b1] = addLAPcuts(x,y,A,B,b,params,ops);
                            % update A,B,b
                            A = [A;A1];
                            B = [B;B1];
                            b = [b;b1];
                            info.cuts.lap=info.cuts.lap+length(b1);
                        end
                    end
                    info.time = toc;
                case 2
                    % ==================================================================
                    % If (x,y) is feasible, then add DC cut type I
                    % Otherwise add LAP cut and if p(x) is not integer, and all entries of x are not 1/2
                    % then add DC cut type II
                    % ==================================================================
                    % Use DCA to get a local minimizer
                    % (Linear relaxation for initial point of DCA)
                    [x,y] = DCA(model,t,s1,u0,params,opsdca);
                    info.nbdca=info.nbdca+1;
                    
                    % If (x,y) is feasible
                    if x'*(1-x) <= ops.zero
                        % add Valide ineq
                        %[A1,B1,b1] =addValidIneq(x,s1,s2);
                        % add DC cut
                        [A1,B1,b1] = addDCcutI(x,s1,s2);
                        % update A,B,b
                        A = [A;A1];
                        B = [B;B1];
                        b = [b;b1];
                        info.cuts.dc1=info.cuts.dc1+1;
                        
                        % update the upper bound (UB) if the new feasible solution is better
                        fval = c'*[x;y];
                        if fval < UB
                            UB = fval;
                            xopt = x;
                            yopt = y;
                            fopt = UB;
                        end
                    else % If the (x,y) is infeasible
                        % use proc-V to get a critical point u_new and
                        % update t (not used yet)
                        % add a lift and project cut
                        [A1,B1,b1] = addLAPcuts(x,y,A,B,b,params,ops);
                        % update A,B,b
                        A = [A;A1];
                        B = [B;B1];
                        b = [b;b1];
                        info.cuts.lap=info.cuts.lap+length(b1);
                        pval = sum(min(x,1-x));
                        if abs(pval-round(pval)) >= ops.zero && norm(x-0.5,inf) > ops.zero
                            % add DC cut type II
                            [A1,B1,b1] = addDCcutII(x,s1,s2);
                            % update A,B,b
                            A = [A;A1];
                            B = [B;B1];
                            b = [b;b1];
                            info.cuts.dc2=info.cuts.dc2+1;
                        end
                    end
                    info.time = toc;
                case 3
                    % ==================================================================
                    % If (x,y) is feasible, then add DC cut type I
                    % Otherwise, add LAP and GMIC, and if p(x) is not integer, and all entries
                    % of x are not 1/2, then add DC cut type II
                    % ==================================================================
                    
                    % Use DCA to get a local minimizer
                    % (Linear relaxation for initial point of DCA)
                    [x,y,~,idxB_dca] = DCA(model,t,s1,u0,params,opsdca);
                    info.nbdca=info.nbdca+1;
                    
                    % If (x,y) is feasible
                    if x'*(1-x) <= ops.zero
                        % add DC cut
                        [A1,B1,b1] = addDCcutI(x,s1,s2);
                        % update A,B,b
                        A = [A;A1];
                        B = [B;B1];
                        b = [b;b1];
                        info.cuts.dc1=info.cuts.dc1+1;
                        
                        % update the upper bound (UB) if the new feasible solution is better
                        fval = c'*[x;y];
                        if fval < UB
                            UB = fval;
                            xopt = x;
                            yopt = y;
                            fopt = UB;
                        end
                    else % If the (x,y) is infeasible
                        % use proc-V to get a critical point u_new and
                        % update t (not used yet)
                        % add a lift and project cut
                        [A1,B1,b1] = addLAPcuts(x,y,A,B,b,params,ops);
                        info.cuts.lap=info.cuts.lap+length(b1);
                        % add GMIR cut
                        [A2,B2,b2] = addGMIRcuts(x,y,A,B,b);
                        info.cuts.mir=info.cuts.mir+1;
                        % add gomory cut
                        %[A3,B3,b3] = addGOMcuts(x,y,A,b,idxB_dca,ops);
                        %info.cuts.gom=info.cuts.gom+length(b2);
                        % update A,B,b
                        A = [A;A1;A2];
                        B = [B;B1;B2];
                        b = [b;b1;b2];
                        % add DC cut type II
                        pval = sum(min(x,1-x));
                        if abs(pval-round(pval)) >= ops.zero && norm(x-0.5,inf) > ops.zero
                            % add DC cut type II
                            [A1,B1,b1] = addDCcutII(x,s1,s2);
                            % update A,B,b
                            A = [A;A1];
                            B = [B;B1];
                            b = [b;b1];
                            info.cuts.dc2=info.cuts.dc2+1;
                        end
                    end
                    info.time = toc;
                case 5
                    % ==================================================================
                    % LAP only
                    % ==================================================================
                    info.time=toc;
                case 11
                    % ==================================================================
                    % Parallel version of method 1
                    % ==================================================================
                    
                    % initialize number of cuts to be created parallelly
                    CUTS_A=cell(numOfDCA,1);
                    CUTS_B=cell(numOfDCA,1);
                    CUTS_b=cell(numOfDCA,1);
                    
                    XOPT=cell(1,numOfDCA);
                    YOPT=cell(1,numOfDCA);
                    FOPT=inf(1,numOfDCA);
                    r=0.5;
                    DCA_points = [u0,MultiRandChoose_r(u0,r,s1,s2,numOfDCA-1)]; % choose points around u0 in a ball of radius r=0.5
                    
                    nb_dca = zeros(1,numOfDCA);
                    nb_dccut1 = zeros(1,numOfDCA);
                    nb_dccut2 = zeros(1,numOfDCA);
                    nb_lapcut = zeros(1,numOfDCA);
                    
                    % Run DCA in parallel
                    parfor i = 1:numOfDCA
                        [x,y] = DCA(model,t,s1,DCA_points(:,i),params,opsdca);
                        nb_dca(i) = nb_dca(i) + 1;
                        % If (x,y) is feasible
                        if x'*(1-x) <= tolzero
                            % add DC cut
                            [CUTS_A{i},CUTS_B{i},CUTS_b{i}] = addDCcutI(x,s1,s2);
                            nb_dccut1(i) = nb_dccut1(i) + 1;
                            
                            % update the upper bound (UB) if the new feasible solution is better
                            fval = c'*[x;y];
                            if fval < UB
                                XOPT{i} = x;
                                YOPT{i} = y;
                                FOPT(i) = fval;
                            end
                        else % If the (x,y) is infeasible
                            % use proc-V to get a critical point u_new and
                            % update t (not used yet)
                            pval = sum(min(x,1-x));
                            if abs(pval-round(pval)) >= tolzero && norm(x-0.5,inf) > tolzero
                                % add DC cut type II
                                [CUTS_A{i},CUTS_B{i},CUTS_b{i}] = addDCcutII(x,s1,s2);
                                nb_dccut2(i)=nb_dccut2(i)+1;
                            else
                                % add a lift and project cut
                                [CUTS_A{i},CUTS_B{i},CUTS_b{i}] = addLAPcuts(x,y,A,B,b,params,ops);
                                nb_lapcut(i)=nb_lapcut(i)+length(CUTS_b{i});
                            end
                        end
                    end
                    
                    % check and update best UB, xopt, yopt, fopt, A,B,b
                    % there may exists many redundant cuts, remove them
                    ALL = unique([cell2mat(CUTS_A),cell2mat(CUTS_B),cell2mat(CUTS_b)],'row');
                    A = [A;ALL(:,1:s1)];
                    B = [B;ALL(:,s1+1:end-1)];
                    b = [b;ALL(:,end)];
                    [minfopt,idx]=min(FOPT);
                    if minfopt<inf % if there is new UB
                        fopt=minfopt;
                        xopt = XOPT{idx};
                        yopt = YOPT{idx};
                        UB = fopt;
                    end
                    info.cuts.dc1 = info.cuts.dc1 + sum(nb_dccut1);
                    info.cuts.dc2 = info.cuts.dc2 + sum(nb_dccut2);
                    info.cuts.lap = info.cuts.lap + sum(nb_lapcut);
                    info.nbdca = info.nbdca + sum(nb_dca);
                    info.time = toc;
                case 12
                    % ==================================================================
                    % Parallel version of method 2
                    % ==================================================================
                    
                    % initialize number of cuts to be created parallelly
                    CUTS_A=cell(numOfDCA,1);
                    CUTS_B=cell(numOfDCA,1);
                    CUTS_b=cell(numOfDCA,1);
                    
                    XOPT=cell(1,numOfDCA);
                    YOPT=cell(1,numOfDCA);
                    FOPT=inf(1,numOfDCA);
                    DCA_points = [u0,MultiRandChoose_r(u0,0.1,s1,s2,numOfDCA-1)];
                    
                    nb_dca = zeros(1,numOfDCA);
                    nb_dccut1 = zeros(1,numOfDCA);
                    nb_dccut2 = zeros(1,numOfDCA);
                    nb_lapcut = zeros(1,numOfDCA);
                    
                    % Run DCA in parallel
                    parfor i = 1:numOfDCA
                        [x,y] = DCA(model,t,s1,DCA_points(:,i),params,opsdca);
                        nb_dca(i) = nb_dca(i) + 1;
                        % If (x,y) is feasible
                        if x'*(1-x) <= ops.zero
                            % add DC cut
                            [CUTS_A{i},CUTS_B{i},CUTS_b{i}] = addDCcutI(x,s1,s2);
                            nb_dccut1(i) = nb_dccut1(i) + 1;
                            
                            % update the upper bound (UB) if the new feasible solution is better
                            fval = c'*[x;y];
                            if fval < UB
                                XOPT{i} = x;
                                YOPT{i} = y;
                                FOPT(i) = fval;
                            end
                        else % If the (x,y) is infeasible
                            % use proc-V to get a critical point u_new and
                            % update t (not used yet)
                            % add a lift and project cut
                            [CUTS_A{i},CUTS_B{i},CUTS_b{i}] = addLAPcuts(x,y,A,B,b,params,ops);
                            nb_lapcut(i)=nb_lapcut(i)+length(CUTS_b{i});
                            pval = sum(min(x,1-x));
                            if abs(pval-round(pval)) >= ops.zero && norm(x-0.5,inf) > ops.zero
                                % add DC cut type II
                                [A1,B1,b1] = addDCcutII(x,s1,s2);
                                CUTS_A{i}=[CUTS_A{i};A1];
                                CUTS_B{i}=[CUTS_B{i};B1];
                                CUTS_b{i}=[CUTS_b{i};b1];
                                nb_dccut2(i)=nb_dccut2(i)+1;
                            end
                        end
                    end
                    
                    % check and update best UB, xopt, yopt, fopt, A,B,b
                    A = [A;cell2mat(CUTS_A)];
                    B = [B;cell2mat(CUTS_B)];
                    b = [b;cell2mat(CUTS_b)];
                    [minfopt,idx]=min(FOPT);
                    if minfopt<inf
                        fopt=minfopt;
                        xopt = XOPT{idx};
                        yopt = YOPT{idx};
                        UB = fopt;
                    end
                    info.cuts.dc1 = info.cuts.dc1 + sum(nb_dccut1);
                    info.cuts.dc2 = info.cuts.dc2 + sum(nb_dccut2);
                    info.cuts.lap = info.cuts.lap + sum(nb_lapcut);
                    info.nbdca = info.nbdca + sum(nb_dca);
                    info.time = toc;
            end
            
            % add LAP cuts at initial point
            [A1,B1,b1] = addLAPcuts(x0,y0,A,B,b,params,ops);
            info.cuts.lap=info.cuts.lap+length(b1);
            if method==3
                % add gomory cut
                %[A2,B2,b2] = addGOMcuts(x0,y0,A,b,idxB,ops);
                %info.cuts.gom=info.cuts.gom+length(b2);
                % update A,B,b
                %A = [A;A1;A2];
                %B = [B;B1;B2];
                %b = [b;b1;b2];
                [A2,B2,b2] = addGMIRcuts(x0,y0,A,B,b);
                info.cuts.mir=info.cuts.mir+1;
                % update A,B,b
                A = [A;A1;A2];
                B = [B;B1;B2];
                b = [b;b1;b2];
            else
                A = [A;A1];
                B = [B;B1];
                b = [b;b1];
            end
            
            % show and print
            [info.gap,info.clgap]=updategap(UB,LB,initfval,bestknown);
            printfverbose(verboseflag,UB,LB,info);
            plotiters(plotflag,LB_old,LB,UB_old,UB,info.iter);
            info.iter = info.iter + 1;
        end
    end
    % check solution status
    % if maxiterDCCUT exceed, return the current best solution
    if info.iter>= maxiterDCCUT
        info.exitflag=2;
    elseif info.time >= maxtime
        info.exitflag=3;
    else % else optimal solution found
        info.exitflag=1;
    end
    printfinalinfo(verboseflag,UB,LB,fopt,info);
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Show version and copyrights
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function showcopyrights()
    fprintf('%s\n','*******************************************************************');
    fprintf('%s\n', 'DC_CUT algorithm for solving Mixed Binary Linear Program (ver 1.0)');
    fprintf('%s\n\n','*******************************************************************');
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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% randomly choosing n points around u0 with radius r
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% u0 -> the center
% r -> the radius for y
% n -> the number of chosen points
function y = MultiRandChoose_r(u0,r,s1,s2,n)
    y = [rand(s1,n);repmat(u0(s1+1:end),1,n)];
    lambda = rand(s2,n);
    y(s1+1:end,:)=y(s1+1:end,:)+lambda.*(-r) + (1-lambda).*r;
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% randomly choosing n points in [a,b]
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% a -> lower bound
% b -> upper bound
% n -> the number of chosen points
function y = MultiRandChoose(a,b,n)
    s=length(a);
    y = zeros(s,n);
    for j = 1:n
        lambda = rand(s,1);
        y(:,j)=lambda.*a + (1-lambda).*b;
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Add LAP cuts
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [A1,B1,b1] = addLAPcuts(x,y,A,B,b,params,ops)
    s1=length(x);
    L=1:s1;
    D = abs(x-round(x));
    I = L( D >=0.001);
    [~,idx] = maxk(D(I),ops.numOfLAP); % find index of the biggest numOfLAP elements
    index = I(idx);
    num = length(index);
    CUTS_A = cell(1,num);
    CUTS_B = cell(1,num);
    CUTS_b = cell(1,num);
    
    % we can generate LAP cuts by parallel procedure or not
    if ops.paralLAP % parallel proc
        parfor i = 1:num
            j = index(i);
            [nu,zeta] = LAP(-[A,B],-b,[x;y],j,params);
            CUTS_A{i} = -nu(1:s1);
            CUTS_B{i} = -nu(s1+1:end);
            CUTS_b{i} = -zeta;
        end
    else % non-parallel proc
        for i = 1:num
            j = index(i);
            [nu,zeta] = LAP(-[A,B],-b,[x;y],j,params);
            CUTS_A{i} = -nu(1:s1);
            CUTS_B{i} = -nu(s1+1:end);
            CUTS_b{i} = -zeta;
        end
    end
    A1 = cell2mat(CUTS_A)';
    B1 = cell2mat(CUTS_B)';
    b1 = cell2mat(CUTS_b)';
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Add GMIR cuts
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [A1,B1,b1] = addGMIRcuts(x0,y0,A,B,b)
    s1 = length(x0);
    [w,s]=computeBaseIneq([x0;y0],[A,B],b);
    
    [lhs, rhs] = GMIR(w(1:s1),w(s1+1:end),s,1); % create GMIR cut
    A1 = lhs(1:s1);
    B1 = lhs(s1+1:end);
    b1 = rhs;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Linear part of valide inequality
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function  [c] = getLin(x,s1)
    J0 = x <= 1/2;
    c=-ones(1,s1);
    c(J0) = 1;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Valide inequality
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function  [A1,B1,b1] = addValidIneq(x,s1,s2)
    J0 = x <= 1/2;
    A1=ones(1,s1);
    A1(J0) = -1;
    B1 = zeros(1,s2);
    b1 = s1-sum(J0); % |J1|
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Add DC cut type I
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function  [A1,B1,b1] = addDCcutI(x,s1,s2)
    J0 = x <= 1/2;
    A1=ones(1,s1);
    A1(J0) = -1;
    B1 = zeros(1,s2);
    b1 = s1-sum(J0)-1; % |J1|-1
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Add DC cut type II
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function  [A1,B1,b1] = addDCcutII(x,s1,s2)
    J0 = x <= 1/2;
    J1 = not(J0);
    A1=ones(1,s1);
    A1(J0) = -1;
    B1 = zeros(1,s2);
    b1 = sum(J1)-ceil(sum(x(J0))+sum(1-x(J1)));
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Update gap and closed gap
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [gap,clgap]=updategap(UB,LB,initfval,bestknown)
    % gap (%) is computed by (UB-LB)/((max(abs(UB),abs(LB))+1)
    % closed gap (%) is computed by 100*(LB-initfval)/(bestknown-initfval)
    gap = 100*(UB-LB)/(max(abs(UB),abs(LB))+1);
    clgap = nan;
    if ~isnan(bestknown)
        clgap = 100*(LB-initfval)/(bestknown-initfval);
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Print verbose
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function printfverbose(verboseflag,UB,LB,info)
    if verboseflag==1
        fprintf('%5d | %10.4f | %10.4f | %8d | %8d | %8d | %8d | %8d | %8.2f | %9.2f | %8.4f\n',info.iter,UB,LB,info.cuts.dc1,info.cuts.dc2,info.cuts.lap,info.cuts.mir,info.nbdca,info.gap,info.clgap,info.time);
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Plot iterations
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function plotiters(plotflag,LB_old,LB,UB_old,UB,iter)
    if plotflag==1
        myplot(LB_old,LB,UB_old,UB,iter);
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Print final info
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function printfinalinfo(verboseflag, UB, LB, fopt, info)
    if verboseflag ~=0
        fprintf('*******************************************************************\n');
        fprintf('* Summary:\n');
        fprintf('* Optimal value: %.5e; Solution time: %.3f seconds; Status: %d\n',fopt,info.time, info.exitflag);
        fprintf('* LB: %.5e; UB: %.5e; Gap: %.2f%%; Closed Gap (cl.gap): %.2f%%\n',LB,UB,info.gap,info.clgap);
        fprintf('* Number of iterations: %d; number of DCA: %d\n', info.iter,info.nbdca);
        fprintf('* Number of cuts: DC cut type I (%d); DC cut type II (%d); LAP (%d); MIR (%d).\n',info.cuts.dc1,info.cuts.dc2,info.cuts.lap,info.cuts.mir);
    end
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Add Gomory cuts
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [A1,B1,b1] = addGOMcuts(x,y,A,b,idxB,ops)
    s1=length(x);
    s2=length(y);
    L=1:s1;
    D = abs(x-round(x));
    I = L( D >=0.001); % minimum fractionality
    [~,idx] = maxk(D(I),ops.numOfGOM); % find index of the biggest numOfGOM elements
    index = I(idx);
    
    num = length(index);
    CUTS_A = cell(1,num);
    CUTS_B = cell(1,num);
    CUTS_b = cell(1,num);
    
    % we can generate GOM cuts by parallel procedure or not
    if ops.paralGOM % parallel proc
        parfor i = 1:num
            j = index(i);
            [w,z] = GMIC(x,A,b,idxB,j); % for pure integer only at present
            CUTS_A{i} = w';
            CUTS_B{i} = zeros(s2,length(z));
            CUTS_b{i} = z;
        end
    else % non-parallel proc
        for i = 1:num
            j = index(i);
            [w,z] = GMIC(x,A,b,idxB,j);
            CUTS_A{i} = w';
            CUTS_B{i} = zeros(s2,length(z));
            CUTS_b{i} = z;
        end
    end
    A1 = cell2mat(CUTS_A)';
    B1 = cell2mat(CUTS_B)';
    b1 = cell2mat(CUTS_b)';
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Probability for restart DCA
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function r = ProbaRestart()
    % Todo
end