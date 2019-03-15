function s = lfp_tfa_time2sample( t, t0, ts )
%lfp_tfa_time2sample - compute the sample number given the recorded time of 
%sample, recorded time of reference sample and sampling time. 
%
% USAGE:
%	s = lfp_tfa_time2sample( t, t0, ts )
%
% INPUTS:
%		t       - recorded time of sample (s)
%       t0      - recorded time of reference sample (eg. state onset) (s)
%       ts      - sampling time (s)
% OUTPUTS:
%		s       - sample number
%
    s = round((t - t0)/ts);

end

