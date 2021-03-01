function model = MIP2GRU(MIP)
    % convert MIP model to Gurobi model
    n=length(MIP.C1); % nb of integer variables
    m=length(MIP.C2); % nb of continuous variables
    
    model.modelsense = 'min';
    model.obj = [MIP.C1;MIP.C2];    
    model.A = sparse([MIP.A,MIP.B]);
    model.rhs = MIP.b;
    model.sense = [repmat('<', size(model.A,1), 1)];    
    model.lb = zeros(n+m,1);
    % set upper bound of y if exists
    if isfield(MIP,'uby')
        model.ub=[ones(n,1);MIP.uby];
    else
        model.ub=[ones(n,1);inf(m,1)];
    end
    model.vtype=[repmat('B', 1, n),repmat('C', 1, m)];
    model.varnames=cell(n+m,1);
    for i=1:n
        model.varnames{i}=['x',num2str(i)];
    end
    for j=1:m
        model.varnames{j+n}=['y',num2str(j)];
    end
end