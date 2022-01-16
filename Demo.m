clear ALL
clc
[ ~,v,~,~,obj ] = MGF2WL( X',Ai,lambda2,lambda1,c);
% X: Original data
% Ai: multiple graphs
% lambda2,lambda1: two paras
% c no. of classes
[~,idx]=sort(v,'descend');
% idx is the feature index sorted in decending order