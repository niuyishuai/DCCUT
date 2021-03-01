% Example 1 in the paper
x = sdpvar(3,1);
K=[x>=0;x<=1/2];
plot(K);

%%
% g = -2x1-x2-x3. h = -13*sum min{x_i,1-xi}

c=[0;0;0];
t=1;
f=c'*x+t*sum(min(x,1-x));
xk=[0.5;0.5;0.5];
assign(x,xk);
fk=value(f);
iter = 0;
ops=sdpsettings('solver','gurobi','verbose',0);

while(1)
    iter=iter+1;
    u=ones(3,1);
    u(xk<0.5)=-1;
    idx=[xk==0.5];
    %u(xk==0.5)=2*rand(sum(idx),1)-1;
    u(idx)=-1;
    optimize(K,-(-c+t*u)'*x,ops);
    
    xopt = value(x);
    fopt = value(f);
    dx = norm(xopt-xk)/(1+norm(xopt));
    df = abs(fk - fopt)/(abs(fopt)+1);
    fprintf('iter: %d; f: %e; delta_x: %e; delta_f: %e\n',iter, fopt,dx,df);
    %xopt
    if df < 1e-8 && dx < 1e-8
        break
    end
    xk=xopt;
    fk=fopt;
end
xopt