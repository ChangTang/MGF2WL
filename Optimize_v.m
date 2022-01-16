function [ v ] = Optimize_v( W )

dd=sqrt(sum(W.^2,2));
v=dd./sum(dd);
v=v';
end

