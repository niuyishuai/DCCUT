function ops = DCCUT_options()
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
    ops.verbose = 1;
    ops.method = 1;
    ops.plot = 0;
    ops.t = 100;
    ops.maxiterDCCUT = 1000;
    ops.maxtime = 60;
    ops.dca.maxiterDCA = 200;
    ops.dca.tolf = 1e-6;
    ops.dca.tolx = sqrt(ops.dca.tolf);
    ops.tolgap = 0.01;
    ops.zero = 1e-9;
    ops.params_gurobi.Method=-1;
    ops.params_gurobi.OutputFlag = 0;
    ops.numOfDCA = 1;
    ops.numOfLAP = 1;
    ops.paralLAP = false;
    ops.numOfGOM = 1;
    ops.paralGOM = false;
end
