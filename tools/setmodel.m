% Set model
function model = setmodel(c,A,B,b)
model = struct();
model.modelsense = 'min';
model.obj = c;
model.A = sparse([A,B]);
model.rhs = b;
model.sense = repmat('<',size(A,1),1);
end
