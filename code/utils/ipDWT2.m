function y = ipDWT2(X,R,M,p,N)
% y = ipDWT2(X,M,K)
% Inverse windowed wavelet transform.
% Input:
% X - WDWT produced by 'pDWT2'
% R - block length
% M - over-sampling rate
% K - window shape parameter
% N - length of signal
% Note: Inversion of first and last block of signal is not implemented


[h,g] = compute_wavelet_filter('Daubechies',p);

[Nfft, Nc] = size(X);                   % get size
Z = zeros(size(X));

Jmax = log2(Nfft)-1; Jmin = 0;
for i=1:Nc
    f1 = X(:,i);
    for j=Jmin:Jmax
        a = f1(1:2^j);
        d = f1(2^j+1:2^(j+1));
        a = cconvol(upsampling(a,1),reverse(h),1);
        d = cconvol(upsampling(d,1),reverse(g),1);
        f1(1:2^(j+1)) = a + d;
    end
    Z(:,i) = f1;
end


y = zeros(1,R/M*(Nc+M-1));
i = 0;
for k = 1:Nc
    y(i + (1:R)) = y(i + (1:R)) + Z(:,k).';
    i = i + R/M;
end
y = (1/M).*y(R+(1:N));

end