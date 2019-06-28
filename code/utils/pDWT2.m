function X = pDWT2(x,R,M,p,Nfft)
% X = pDWT2(x,R,M,K,Nfft)
% Parseval windowed wavelet transform with 100x(1-1/M) percent overlap
% (M-times overcomplete) with window shape parameter, K
% Each block is multiplied by a cosine window to power K.
% Input:
% x - 1D signal
% R - block length
% M - over-sampling rate
% K - window shape parameter
% Nfft - length of FFT (Nfft >= R)
% Note: R should be an integer multiple of M.
% Note: K should be less than M (1 <= K < M).
%

if ~isrow(x)
    x = x';
end

[h,g] = compute_wavelet_filter('Daubechies',p);
% disp(['h filter = [' num2str(h) ']']);
% disp(['g filter = [' num2str(g) ']']);

x = [zeros(1,R) x zeros(1,R)];          % to deal with first and last block
Y = buffer(x,  R, R*(M-1)/M, 'nodelay');
X = zeros(size(Y));
[~,n] = size(X);

Jmax = log2(Nfft)-1; Jmin = 0;
for i=1:n
    X(:,i) = Y(:,i);
    for j=Jmax:-1:Jmin
        a1 = X(1:2^(j+1),i);
        a = subsampling(cconvol(a1,h));
        d = subsampling(cconvol(a1,g));
        X(1:2^(j+1),i) = cat(1, a, d );
    end

end

end

