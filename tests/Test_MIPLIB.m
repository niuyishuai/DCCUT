%%
parpool('local',30);
%%
pathdir = '/cluster/home/1134203/MIPLIB/';
modelfilelst = {'mas74.mps','neos5','mad','pk1'};
modelfoptlst = [11801.18572, 15, 0.0268, 11];
ops=DCCUT_options;
ops.t=500;
ops.maxiterDCCUT=inf;
ops.maxiterDCA=1000;
ops.verbose=1; % 1 display details 0 silence
ops.plot=0; % 1 plot iterations, 0 no plot

methodlst=[5,1,2,3,5,1,2,3];
nlaplst= 1:2:30;
    
for idx = 1:length(modelfilelst)
    model = gurobi_read([pathdir,modelfilelst{idx}]);
    MIP = transform_to_my_model(model);
    MIP.fopt = modelfoptlst(idx);
    
    for maxtime = [30,60,120]
        ops.maxtime = maxtime;
        fprintf('*************************\n');
        fprintf('maxtime = %4d\n',ops.maxtime);
        fprintf('*************************\n');
        R=cell(length(methodlst),length(nlaplst));
        for i = 1:length(nlaplst)
            for j = 1:length(methodlst)
                if j<=length(methodlst)/2 % parallel algos
                    ops.paralLAP=true;
                else % nonparallel algos
                    ops.paralLAP=false;
                end
                ops.method=methodlst(j);
                ops.numOfLAP=nlaplst(i);
                
                ops.numOfDCA = 1;
                [xopt,yopt,fopt,info] = DCCUT(MIP,ops);
                fprintf('example = %6s, method = %d, LAP = %d, Paral = %d, gap = %.2f, clgap = %.2f, ub = %.5f, time = %.3f\n',modelfilelst{idx},ops.method,ops.numOfLAP,ops.paralLAP,info.gap,info.clgap,fopt,info.time);
                R{j,i} = info;
                R{j,i}.fopt = fopt;
            end
        end
        save(['res_miplib_',modelfilelst{idx},'_',num2str(maxtime)],'R');
    end
end
