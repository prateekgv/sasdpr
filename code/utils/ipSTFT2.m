function y = ipSTFT2(X,R,M,K,N)
% y = ipSTFT2(X,M,K)
% Inverse Short-time Fourier Transform (inverse of pSTFT2).
% Input:
% X - STFT produced by 'pSTFT2'
% R - block length
% M - over-sampling rate
% K - window shape parameter
% N - length of signal
% Note: Inversion of first and last block of signal is not implemented

[Nfft, Nc] = size(X);                   % get size

n = (1:R) - 0.5;
win  = sin(pi*n/R).^K;                  % make cosine window

NC = sqrt(sum(win.^2) * M * Nfft/R);    % normalization constant

Y = ifft(X);                            % inverse FFT of each column of X
Y = Y(1:R,:);                           % truncate down to block length
Y = bsxfun(@times, Y, win');

y = zeros(1,R/M*(Nc+M-1));
i = 0;
for k = 1:Nc
    y(i + (1:R)) = y(i + (1:R)) + Y(:,k).';
    i = i + R/M;
end

y = NC * y(R+(1:N));
y = y*(2/M);

if K == 1
    A = 1;
elseif K == 2
    A = 4/3;
elseif K == 3
    A = 8/5;
elseif K == 4
    A = 2^6/(5*7);
elseif K == 5
    A = 128/63;
elseif K == 6
    A = 2^9 / (3*7*11);
elseif K == 7
    A = 2^10 / (3*11*13);
else
    disp('Implmented only for K = 1,...,7')
end

y = A*y;
    