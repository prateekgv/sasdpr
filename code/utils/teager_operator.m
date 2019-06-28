%==========================================================================
% function [tk_er] = teager_operator(sig)
%==========================================================================
% @description: Compute the Teager-Keiser energy operator of the input
%               signal. In this case, the input signal is the oscillatory
%               signal.
% @author     : Prateek Gundannavar
% @date       : 06/12/18
%
% @input
%             - sig         Signal of interest
% @output
%             - tk_er       Output of the Teager-Keiser energy operator
%==========================================================================
function [tk_er] = teager_operator(sig)

sig=sig(:);

squ1   = sig(2:length(sig)-1).^2;
oddi1  = sig(1:length(sig)-2);
eveni1 = sig(3:length(sig));
tk_er  = squ1 - (oddi1.*eveni1);

%make it the same length
tk_er  = [tk_er(1); tk_er; tk_er(length(sig)-2)]; 

end