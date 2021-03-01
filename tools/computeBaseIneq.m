function [w,s] = computeBaseIneq(x0,A,b)
    % find active set
    idx_act = abs(A*x0-b)<1e-8;
    A_act = A(idx_act,:);
    
    % get linear part
    w = sum(A_act);
    s = w*x0;
    if abs(s-round(s))<1e-8
        idxz = w==0;
        mm=min(abs(w(~idxz)));
        w=w/mm;
        s=w*x0;
        if abs(s-round(s))<1e-8
            w=w/(2*s);
            s=0.5;
        end
        %ff=max(abs(w))+1;
        %w=w/ff;
        %s=s/ff;
    end
end