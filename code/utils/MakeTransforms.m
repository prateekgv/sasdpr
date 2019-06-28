
function [A, AH, normA] = MakeTransforms(type, N, params)

% [A, AH] = MakeTransforms(type, N, params)
% with A AH = I
%
% INPUT
%    N : signal length
%    type : type of transform can be 'DFT' or 'STFT'
%    params : parameters for transform
%
%    For DFT, params = [Nfft]
%       Nfft = number of DFT coefficients, Nfft >= N
%       Nfft/N is over-sampling factor
%
%   For STFT, params = [R, Nfft]
%       R = length of frame
%       Nfft = number of DFT coefficients, Nfft >= R
%
%    For DCT, params = [Ndct]
%       Ndct = number of DCT coefficients, Ndct >= N
%       Ndct/N is over-sampling factor
%
%   For STDCT, params = [R, Ndct]
%       R = length of frame
%       Ndct = number of DCT coefficients, Ndct >= R
%
% OUTPUT
%    A, AH : function handles for transform
%       AH is the conjugate transpose of A
%       Note: input to AH must be a row vector of length N

% DFT : Discrete Fourier transform with zero-padding
% STFT : Short-time Fourier transform
% DCT : Discrete cosine transform with zero-padding
% STDCT : Short-time DCT
% DWT2 : Windowed discrete wavelet transform (similar to STFT but for wavelets)
% Ivan Selesnick, 2015
% 
% Modified by Prateek Gundannavar, 2019.
% email: prateek@wustl.edu

switch type
    case 'DFT'

        % over-sampled DFT (oversampling factor = Nfft/N)
        Nfft = params(1);
        truncate = @(x, N) x(1:N);
        AH = @(x) fft(x, Nfft)/sqrt(Nfft);
        A = @(X) truncate(ifft(X), N)*sqrt(Nfft);
        
        normA = sqrt(N/Nfft);

    case 'STFT'
        
        R = params(1);
        M = params(2);
        K = params(3);
        Nfft = params(4);
        AH = @(x) pSTFT2(x, R, M, K, Nfft);
        A = @(x) ipSTFT2(x,R,M,K,N);       

        normA = sqrt(R/(M*Nfft));

    case 'DCT'

        % over-sampled DFT (oversampling factor = Nfft/N)
        Ndct = params(1);
        truncate = @(x, N) x(1:N);
        AH = @(x) dct(x, Ndct);
        A = @(X) truncate(idct(X), N);
        
        normA = sqrt(N/Ndct);
        
    case 'STDCT'
        
        R = params(1);
        Ndct = params(2);
        AH = @(x) stdct(x, R, Ndct);
        A = @(x) istdct(x, R, N);       

        normA = sqrt(R/(2*Ndct));
        
    case 'DWT2'
        
        R = params(1);
        M = params(2);
        p = params(3);
        Nfft = params(4);
        AH = @(x) pDWT2(x,R,M,p,Nfft);
        A = @(x) ipDWT2(x,R,M,p,N);
        
        normA = sqrt(R/(M*Nfft));
        
        
end

% x = randn(1, N);         % assuming row vector ! ! !

% % Verify perfect reconstruction property of transform 1 (A*AH = I)
% c = AH(x);
% err = x - A(c);
% re = max(abs(err(:)));          % should be zero
% fprintf('Reconstruction error = %f\n', re)
% 
% % Verify Parseval energy identity
% E = sum(abs(x(:)).^2);
% E1 = sum(abs(c(:)).^2);         % should be unity
% fprintf('Energy ratio = %f\n', E/E1)

