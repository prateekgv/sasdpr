function X = pSTFT2(x,R,M,K,Nfft)
% X = pSTFT2(x,R,M,K,Nfft)
% Parseval short-time Fourier Transform with 100x(1-1/M) percent overlap
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
% % Example:
% [s,fs] = wavread('sp1.wav');
% N = 20000;
% x = s(1:N)';
% R = 501;
% M = 3;
% K = 2;
% Nfft = 512;
% X = pSTFT2(x,R,M,K,Nfft);
% y = ipSTFT2(X,R,M,K,N);
% max(abs(x - y))       % verify perfect reconstruction
%
% sum(abs(x(:)).^2)     % signal energy
% sum(abs(X(:)).^2)     % STFT energy (should be equal)

x = x(:).';                             % ensure x is row vector
n = (1:R) - 0.5;
win  = sin(pi*n/R).^K;                  % cosine window
NC = sqrt(sum(win.^2) * M * Nfft/R);    % normalization constant
x = [zeros(1,R) x zeros(1,R)];          % to deal with first and last block
X = buffer(x,  R, R*(M-1)/M, 'nodelay');
X = bsxfun(@times, win', X);
X = fft(X, Nfft)/NC;                    % FFT applied to each column of X

