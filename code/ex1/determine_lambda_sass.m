function [ lam ] = determine_lambda_sass(A, B, N, K, noise)

S = ones(N,N-K);
S = tril(S,-K);
HTH = B' * ((A*A') \ (B));
lam = max(S'*HTH*noise);


end

