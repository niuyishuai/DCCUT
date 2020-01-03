function [x,y] = my_DCA(C1,C2,t,A,B,b,u,params)

s2 = size(A,2);

done = false;
while  ~done
    u1 = u;
    % var1 and var2: the values of two iterates
    x1 = u1(1:s2); 
    var1 = [C1',C2']*u1 + t*sum(min(x1,1-x1));
    
    % compute the subgradient of h(u) at u1;
    %[v,w] = com_subgrad_h(C1,C2,t,u1);
    v = -ones(s2,1);
    v(find(x1 >= 1/2)) = 1;
    v = t*v;
    % get the coefficients of the objective function
    f = [C1-v;C2];
    
    % call gurobi to solve the program
    % u = call_Gurobi(f,[A,B],b,[],[],[],[],[],options);
    model = setmodel(f,A,B,b);
    
    %params.OutputFlag = 0;
    problem = gurobi(model,params);
    u = problem.x;
    %u = linprog(f,[A,B],b,[],[],[],[],options);
    x = u(1:s2);
    
    var2 =  [C1',C2']*u + t*sum(min(x,1-x));
    
    if abs(var1-var2) < 1e-8 || abs(var2-var1) < 1e-8
        done = true;       % stop check!        
    end
end
y = u(s2+1:end);    
end
