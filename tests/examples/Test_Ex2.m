% Example 2 in the paper
c=[-1;-1];
n=2;
x = sdpvar(n,1);
C1=[-4*x(1)+12*x(2)>=-1;x>=0;x<=1];
C2=[-12*x(1)-4*x(2)>=-13];
K = C1+C2;
plot(K);
V = [0 0; 0 1; 0.25 0; 0.75 1; 1 0.25]; % each row is a vertex
%%
% compute exact penalty parameter t0
alpha0 = min(V*c);
p=sum(min(V,(1-V)),2);
m=min(p(p>0));
minf=min(V(p==0,:)*c);
t0=(minf-alpha0)/m;

% compute big enough parameter t1 for valid inequality
maxf=max(V*c);
minfb=min(V(p>0,:)*c);
M=maxf-minfb;
sigma=1;
t1=M/sigma;

%%
t=4;
f=c'*x+t*sum(min(x,1-x));
xk=[0.25;0];
assign(x,xk);
fk=value(f);
iter = 0;
ops=sdpsettings('solver','gurobi','verbose',0);

fprintf('-------------Start DCA-------------\n');
while(1)
    iter=iter+1;
    u=ones(n,1);
    u(xk<=0.5)=-1;
    optimize(K,-(-c+t*u)'*x,ops);
    
    xopt = value(x);
    fopt = value(f);
    dx = norm(xopt-xk)/(1+norm(xopt));
    df = abs(fk - fopt)/(abs(fopt)+1);
    fprintf('iter: %d; fval: %e; delta_x: %e; delta_f: %e\n',iter, fopt,dx,df);
    %xopt
    if df < 1e-8 && dx < 1e-8
        break
    end
    xk=xopt;
    fk=fopt;
end
fprintf('-------------End DCA-------------\n');
xopt

%%
% test GMIR cut
obj=export(K);
%A=[4 -12;12 4;1 0;0 1];
%B=[];
%b=[1;13;1;1];
b = obj.rhs;
%lambda = [1,1,0,0]';%
ones(length(b),1); %[1;0;0;2];
%lambda = lambda*0.01;
lambda = lambda/(2*lambda'*b);
[lhs, rhs] = GMIR(A,B,b,lambda)
if norm(lhs)>0
    %K = K +[lhs*x<=rhs];
    plot(K+[lhs*x<=rhs]);
end