% Generalized Mixed-Integer Rounding cut (GMIR)
% X = {(x,y)\in Z+^n \times R+^p: Ax + By <= b}
% Here n>= 1, p >= 0 and b \in R^m.
% Given lambda \in R^m such that lambda'*b \notin Z
%***************************
% Input : A, B, b, lambda
% Output: [lhs, rhs] cut efficients
%***************************
function [lhs,rhs] = GMIR(A,B,b,lambda)
    
    if abs(b'*lambda - round(b'*lambda)) < 1e-8
        error('Wrong Usage');
    end
    
    lambda= lambda(:); % vectorize lambda
    
    lambdam = min(lambda,0); % lambda^-, negative parts of lambda
    
    LA = lambda'*A; % a row vector
    LMA = -lambdam'*A; % a row vector
    LA_hat = LA - floor(LA);
    
    Lb = lambda'*b;
    LMb = -lambdam'*b;
    Lb_hat = Lb - floor(Lb);
    
    part1 = floor(LA) + (max(LA_hat-Lb_hat,0) + LMA)/(1-Lb_hat);
    
    if ~isempty(B)
        LB = lambda'*B;
        LMB = -lambdam'*B;
        part2 = ( min(LB,0) + LMB)/(1-Lb_hat);
    else
        part2=[];
    end
    
    part3 = floor(Lb) + LMb/(1-Lb_hat);
    
    lhs = [part1,part2];
    rhs = part3;
end