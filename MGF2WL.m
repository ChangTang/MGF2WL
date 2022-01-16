function [ W,v,S,alpha,obj ] = MGF2WL( X,Ai,lambda,gamma,c)

[n,~,m]=size(Ai);
[d,~]=size(X);
alpha=ones(1,m)./m;

maxiter=10;
obj = zeros(3,maxiter);
S=sum(Ai,3)./m;
v=ones(1,d)./d;
Serror = 0;
for i=1:maxiter
    disp(['Iter: ',num2str(i)]);
    W  = Optimize_W( X,S,v,gamma,c );
    v=Optimize_v(W); 
    S = Optimize_S(X, W, lambda, Ai,alpha );
    alpha= Optimize_alpha(S,Ai);
    for si=1:m
       Serror = Serror + alpha(si)*gamma*norm(S-Ai(:,:,si),'fro').^2;
    end
    
   d = sum(S,2);D = diag(d);L = D-S;
   traceerror = trace(W'*X*L*X'*W);
   P = inv(diag(v))*W;
   Perror = lambda*norm(P,'fro').^2;
   obj(1,i) = traceerror;
   obj(2,i) = Perror;
   obj(3,i) = Serror;
    
    
end



end

