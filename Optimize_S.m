function [ U ] = Optimize_S(X, Phi, lambda1, Ai,alpha) % update U
lambda = 1/lambda1;
num = size(X,2);
m = size(Ai,3);
PhiX = Phi'*X;
islocal = 1; % only update the similarities of neighbors if islocal=1
dist = L2_distance_1(PhiX,PhiX); % num * num
U = zeros(num);
for i=1:num
    idx = zeros();
    for v = 1:m
        temp = Ai(:,:,v);
        s0 = temp(i,:);
        idx = [idx,find(s0>0)];
    end;
    idxs = unique(idx(2:end));
    if islocal == 1
        idxs0 = idxs;
    else
        idxs0 = 1:num;
    end;
    for v = 1:m
        temp = Ai(:,:,v);
        s1 = temp(i,:);
        si = s1(idxs0);
        di = dist(i,idxs0);
        mw = m*alpha(v);
        lmw = lambda/mw;
        q(v,:) = si-0.5*lmw*di;
    end;
    U(i,idxs0) = SloutionToP19(q,m);
    clear q;
end;
end