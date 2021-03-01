function model = transform_to_my_model(model)

idx_geq = model.sense == '>';
idx_eq = model.sense == '=';
model.A(idx_geq,:) = -model.A(idx_geq,:);
model.rhs(idx_geq) = -model.rhs(idx_geq);

A = model.A;
rhs = model.rhs;
A = [A;-A(idx_eq,:)];
model.b = [rhs;-rhs(idx_eq,:)];

model.sense = repmat('<',size(A,1),1);

Iindex = model.vtype=='I';
Bindex = model.vtype=='B';
idx = Iindex | Bindex; % either I or B

model.A = A(:,idx); % integer part
model.B = A(:,~idx); % continuous part
model.C1 = model.obj(idx);
model.C2 = model.obj(~idx);

model.uby = model.ub(~idx);
model.uby(model.uby == Inf) = 1e+10;

end
