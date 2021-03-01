function [w,b] = LAP(A0,b0,x0,j,params)
    % A0*x0 >= b0 -A0*x0 <= -b0
    [m,dim] = size(A0);
    
    Aj = A0;
    Aj(:,j) = [];
    aj = A0(:,j);
    Bjt = zeros(m,dim);
    Bjt(:,j) = aj-b0;
    Ajt = A0-Bjt;
    Aeq = [Aj',-Aj'];
    
    Aeq = [Aeq;ones(1,2*m)];
    beq = [zeros(dim-1,1);1];
    f = [Bjt*x0;Ajt*x0-b0];
    
    %[X,~] = call_Gurobi(f,-eye(2*m,2*m),zeros(2*m,1),...
    %Aeq,beq,[],[],[],options);
    model = setmodel(f,Aeq,[],beq);
    problem = gurobi(model,params);
    
    x = problem.x;
    
    u = x(1:m);
    v = x(m+1:end);
    w = Bjt'*u + Ajt'*v;
    b = b0'*v; 
end


