function [ W ] = Optimize_W( X,A,v,gamma,nClass )

d1=sum(A);
d2=sum(A,2);
d=d1'+d2;
L=(diag(d)-A-A')./2;
GroupLassoType='LS21';

switch lower(GroupLassoType)
    case lower('JFSSL')
        Y = eig1(L, nClass, 0);
        W = FSSL_subspace(X, Y, gamma);
    case lower('LS21')
        Y = eig1(L, nClass, 0);
        W = LS21_new(X', Y, v,gamma);
    case lower('NDFS') % d^3
        tmp = X * L * X';
        tmp = (tmp + tmp') / 2;
        if exist('W', 'var')
            W = LquadR21_reg(tmp, nClass, gamma, W);
        else
            W = LquadR21_reg(tmp, nClass, gamma);
        end
    case lower('UDFS')
        tmp = X * L * X';
        tmp = (tmp + tmp') / 2;
        W = LquadR21_reg(tmp, nClass, gamma);
    case lower('MCLEASTR')
        Y = eig1(L, nClass, 0);
        W = mcLeastR(X', Y, gamma, struct('rFlag', 1, 'rsL2', 0));
    otherwise
        error('method does not exist!');
end

end

