function par = initialize_par_DCCUT(par)
if ~isfield(par,'verbose')
    par.verbose = 1;
end

if ~isfield(par,'print')
    par.print = 1;
end

if ~isfield(par,'method')
    par.method = 1;
end

if ~isfield(par,'numOfDCA')
    par.numOfDCA = 1;
end

%if isfield(par,'numOfDCA') && par.numOfDCA == 1
%    if par.print
%        fprintf('Changing to parallel version\n');
%    end
%    par.method = 2;
%end

if ~isfield(par,'plot')
    par.plot = 0;
end


if ~isfield(par,'t')
    par.t = 100;
end

if ~isfield(par,'max_iter')
    par.max_iter = 200;
end


    
if ~isfield(par,'params')    
    par.params.OutputFlag = 0;
end
if ~isfield(par,'numOfLAP')
    par.numOfLAP = 1;
end

if ~isfield(par,'eps')
    par.eps = 1e-2;
end

end