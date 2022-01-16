function [ alpha ] = Optimize_alpha(U, Ai) % update U
    % update alpha
    m = size(Ai,3);
    for v = 1:m
        US = U - Ai(:,:,v);
        distUS = norm(US, 'fro')^2;
        if distUS == 0
            distUS = eps;
        end;
        alpha(v) = 0.5/sqrt(distUS);
    end;
end