%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% read miplib files
% min c'x+d'y
% s.t. Ax+By<=b
%      Aeqx+Beqy=beq
%      lb<=(x,y)<=ub;
%      x in Z^{m}_{+}
%
% not that for the use of DCCUT solver, we process
% 1. convert all equalities to inequalities by formulating Aeqx+Beqy=beq
%    as Aeq x+Beq y<=beq and -Aeq x - Beq y<= -beq
% 2. if uby is not finit, then let it bounded by 10.
% 3. ignore the lower bound of y (always be 0).


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function MIP = readMipModel(modelname,reader)
    if nargin < 2
        reader = 'cplex';
    end
    switch reader
        case 'cplex'
            cpx=Cplex;
            cpx.readModel(modelname);
            Model=cpx.Model;
            
            nvars=length(Model.obj);
            idxC=find(Model.ctype=='C'); %index for continuous variables
            idxI=setdiff([1:nvars],idxC); % index for integer variables
            MIP.n=length(idxI);
            MIP.m=length(idxC);
            
            % objective
            MIP.C1=Model.obj(idxI);
            MIP.C2=Model.obj(idxC);
            
            % bounds
            %MIP.lbx=Model.lb(idxI);
            %MIP.ubx=Model.ub(idxI);
            %MIP.lby=Model.lb(idxC);
            uby=Model.ub(idxC);
            % for unbounded case
            uby(uby==inf)=10; % set an upper bound for y if unbounded
            MIP.uby=uby;
            
            % linear constraints
            idxeq=find(Model.rhs==Model.lhs); % index for equality constraints
            idxineq=setdiff(1:length(Model.rhs),idxeq); % index for inequality constraints
            AC=Model.A(:,idxC); % continuous part
            AI=Model.A(:,idxI); % integer part
            MIP.A=AI(idxineq,:);
            MIP.B=AC(idxineq,:);
            MIP.b=Model.rhs(idxineq,1);
            % convert equality to inequality
            if ~isempty(idxeq)
                Aeq = AI(idxeq,:);
                Beq = AC(idxeq,:);
                beq = Model.rhs(idxeq,1);
                
                MIP.A=[MIP.A;Aeq;-Aeq];
                MIP.B=[MIP.B;Beq;-Beq];
                MIP.b=[MIP.b;beq;-beq];
            end
            
            if ~isempty(idxI)
                MIP.type='MILP';
            end
        case 'gurobi'
    end
end
