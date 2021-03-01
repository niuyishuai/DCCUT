function [x, fval, exitflag, output] = gurobimilp(MIP)
%gurobimilp A mixed-integer linear programming interface using Gurobi
%
%   This example is based on the linprog interface defined in the
%   MATLAB Optimization Toolbox. The Optimization Toolbox
%   is a registered trademark of The MathWorks, Inc.
%
%   [x, fval, exitflag, output] = gurobmilp(MIP)
%
%   minimize     f'*x
%   subject to     A*x <= b,
%              lb <=  x <= ub.
%              x \in vtypes
%
%   exitflag containing the status of the optimization. The values for
%   exitflag and corresponding status codes are:
%      1 - OPTIMAL,
%      0 - ITERATION_LIMIT,
%     -2 - INFEASIBLE,
%     -3 - UNBOUNDED.
%

model=MIP2GRU(MIP);

params.outputflag = 1;
result = gurobi(model, params);

if strcmp(result.status, 'OPTIMAL')
    exitflag = 1;
elseif strcmp(result.status, 'ITERATION_LIMIT')
    exitflag = 0;
elseif strcmp(result.status, 'INF_OR_UNBD')
    params.presolve = 0;
    result = gurobi(model, params);
    if strcmp(result.status, 'INFEASIBLE')
        exitflag = -2;
    elseif strcmp(result.status, 'UNBOUNDED')
        exitflag = -3;
    else
        exitflag = nan;
    end
else
    exitflag = nan;
end

if isfield(result, 'x')
    x = result.x;
else
    x = nan(n,1);
end

if isfield(result, 'objval')
    fval = result.objval;
else
    fval = nan;
end
output=result;
end


