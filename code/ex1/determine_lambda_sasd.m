function [ lam ] = determine_lambda_sasd( H, N, K, noise )

S = ones(N,N-K);
S = tril(S,-K);
HTH = H'*H;
lam = max(S'*HTH*noise);


end

