function [xopt,yopt,iter,idxB] = DCA(model,t,sx,u,params_gurobi,ops_dca)
    % DCA algoritm for solving MIP model
    % sx is the size for x    
    done = false;
    uk = u;
    C = model.obj;
    Ct = C';
    iter = 1;
    while  ~done
        xk = uk(1:sx);
        fk = Ct*uk + t*sum(min(xk,1-xk));
        
        % compute the subgradient of h(u) at u;
        vk = -ones(sx,1);
        vk(xk >= 1/2) = 1;
        % get the coefficients of the objective function, note that the
        % coefficients in y is not changed
        model.obj(1:sx) = C(1:sx) -t*vk;
        
        % call gurobi to solve the linear program        
        problem = gurobi(model,params_gurobi);
            
        uopt = problem.x;
        xopt = uopt(1:sx);
        
        fopt =  Ct*uopt + t*sum(min(xopt,1-xopt));
        
        if abs(fk-fopt)/(1+abs(fopt)) < ops_dca.tolf || norm(uopt-u,2)/(1+norm(uopt,2)) < ops_dca.tolx || iter >= ops_dca.maxiterDCA
            done = true;       % stop check!
        else
            uk=uopt;
        end
        iter=iter+1;
    end
    yopt = uopt(sx+1:end);
    if nargout == 4
        idxB = find(problem.vbasis == 0);
    end
end
