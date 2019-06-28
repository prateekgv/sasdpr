%% function [sig_s,sig_d,bin_vec] = detect_roi(tk_er,fs,threshold,type,method)
%
% Inputs:
%   tk_er   - output of the Teager-Kaiser energy operator
%   fs      - sampling rate
%   threshold - threshold for detecting patterns
%   type    - spindle or kcomplex
%   method  - detoks or sapr
%
% Outputs:
%   sig_s   - start of the pattern of interest
%   sig_d   - duration of the pattern of interest
%   bin_vec - binary vector with the detected pattern
%
% Prateek Gundannavar, prateek@wustl.edu, 2018
% Revised 2019
%% ________________________________________________________________________
%%

function [sig_s, sig_d, bin_vec] = detect_roi( tk_er, fs, threshold, type, method )

n = length(tk_er);
y = zeros(1,n);
y(tk_er >= threshold) = 1;
x = 1:n; th = 0.1;

if y(1) == 1 && y(end) == 1
    y = [~y(1), y, ~y(end)];
    x = 1:n+1;
    x_steps = x(abs(diff(y)) > th);
    x_steps = reshape(x_steps, 2,length(x_steps)/2)';
    x_steps(:,2) = x_steps(:,2)-1;
elseif y(1) == 1
    y = [~y(1), y];
    x_steps = x(abs(diff(y)) > th);
    x_steps = reshape(x_steps, 2,length(x_steps)/2)';
    x_steps(:,2) = x_steps(:,2)-1;
elseif y(end) == 1
    y = [y, ~y(end)];
    x_steps = x(abs(diff(y)) > th);
    x_steps = reshape(x_steps, 2,length(x_steps)/2)';
    x_steps(:,1) = x_steps(:,1)+1;
else
    x_steps = x(abs(diff(y)) > th);
    x_steps = reshape(x_steps, 2,length(x_steps)/2)';
end

sig_s = x_steps(:,1);
sig_d = x_steps(:,2)-x_steps(:,1)+1;

% Remove those kcomplex regions which are < 0.5 sec and > 3.0 sec
if (strcmp(type,'kcomplex')) && (strcmp(method,'sapr'))
    ind = find(sig_d./fs < 0.5 | sig_d./fs > 2.25);
elseif (strcmp(type,'kcomplex')) && (strcmp(method,'detoks'))
    ind = find(sig_d./fs < 0.5 | sig_d./fs > 3.0);
end

% Remove those spindles regions which are < 0.5 sec and > 3.0 sec
if (strcmp(type,'spindle')) && (strcmp(method,'sapr'))
    ind = find(sig_d./fs < 0.5 | sig_d./fs > 3.0);
elseif (strcmp(type,'spindle')) && (strcmp(method,'detoks'))
    ind = find(sig_d./fs < 0.5 | sig_d./fs > 3.0);
end


% remove those regions that above duration threshold
if ~isempty(ind)
    sig_d(ind) = [];
    sig_s(ind) = [];
end

% remove slow-wave activity for K-complexes
tk_vec = zeros(n,1);
if (strcmp(type,'kcomplex')) && (strcmp(method,'sapr'))
    [sig_s,sig_d,~] = remove_sws_effect(fs,sig_s,sig_d,tk_er);
end

% compute the binary vector
bin_vec = zeros(1,n);
if ~isempty(sig_s) && ~isempty(sig_d)
    for ii = 1:size(sig_s,1)
        bin_vec(sig_s(ii):sig_s(ii)+sig_d(ii)) = 1;
    end
end


end


function [sig_s,sig_d,tk_vec] = remove_sws_effect(fs,sig_s,sig_d,tk_er)

% create a new tk_vector
n = length(tk_er);
tk_vec = zeros(n,1);
t = linspace(0,n/fs,n);

if (~isempty(sig_s)) && (~isempty(sig_d))
    for ii = 1:size(sig_s,1)
        tk_vec(sig_s(ii):sig_s(ii)+sig_d(ii)-1) = ...
            tk_er(sig_s(ii):sig_s(ii)+sig_d(ii)-1);
    end
    % remove those peaks which appear within 1.5 seconds duration
    [~,locs] = findpeaks(tk_vec,t);
    rem_locs = find(diff(locs) <= 1.5) + 1;
    if (~isempty(rem_locs))
        locs(rem_locs) = [];
    end
    
    % select only those locations
    flag = zeros(size(sig_d,1),1);
    if (~isempty(locs))
        for ii = 1:length(locs)
            for jj = 1:size(sig_s,1)
                if (locs(ii) >= sig_s(jj)/fs) && (locs(ii) <= ...
                        (sig_s(jj)+sig_d(jj))/fs)
                    flag(jj,1) = 1;
                    break;
                end
            end
        end
        ind = find(flag == 0);
        if ~isempty(ind)
            sig_d(ind) = [];
            sig_s(ind) = [];
        end
    end
    
end

end
