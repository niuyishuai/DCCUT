function premodel = presolvemilp(MIP)
    %presolvemilp using Gurobi
    %
    %   premodel = presolvemilp(MIP)
    %
    %   minimize     f'*x
    %   subject to     A*x <= b,
    %              lb <=  x <= ub.
    %              x \in vtypes
    %
    %
    
    model=MIP2GRU(MIP);
    
    params.outputflag = 0;
    premodel = gurobi_presolve(model, params);
end


