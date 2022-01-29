function [ lfp_tfr_normalized ] = lfp_tfa_baseline_normalization( raw_tfs, cfg_baseline )
% lfp_tfa_baseline_normalization - Baseline normlaization of an LFP time frequency
% spectrogram
% USAGE:
%	[ lfp_tfr_normalized ] = lfp_tfa_baseline_normalization( raw_tfs, cfg_baseline )
%
% INPUTS:
%		raw_tfs             - raw lfp power spectrogram (1xMxN array, M =
%		number of freq bins, N = number of time bins)
%       cfg_baseline        - structure containing baseline information 
%            method     : basline normalization method. 
%            Can be 'subtraction', 'division', 'zscore' or 'relchange'
%            'subtraction'
%               P_n(t,f) = P(t,f) - mean
%            'division'
%               P_n(t,f) = P(t,f) / mean
%            'zscore'
%               P_n(t,f) = (P(t,f) - mean) / std
%            'relchange'
%               P_n(t,f) = (P(t,f) - mean) / mean
%            mean       : mean power in each frequency during the baseline period
%            std        : standard deviation of baseline power in each
%            frequency during the baseline period. Required if method is
%            'zscore'
% OUTPUTS:
%		lfp_tfr_normalized  - baseline normalized LFP power spectrogram
%
% See also lfp_tfa_get_state_tfs, lfp_tfa_plot_site_average_tfr, lfp_tfa_process_LFP

    base_mean = cfg_baseline.mean;
    base=repmat(base_mean, [size(raw_tfs, 1) 1 size(raw_tfs, 3)]);
    if strcmp(cfg_baseline.method , 'zscore')
       base_std =  cfg_baseline.std;
    end
    if strcmp(cfg_baseline.method , 'subtraction')
        lfp_tfr_normalized = raw_tfs - base;
    elseif strcmp(cfg_baseline.method , 'division')
        lfp_tfr_normalized = raw_tfs ./ base;
    elseif strcmp(cfg_baseline.method , 'relchange')
        lfp_tfr_normalized = (raw_tfs - base) ./ base;
    elseif strcmp(cfg_baseline.method , 'zscore')
        lfp_tfr_normalized = (raw_tfs - base) ./ repmat(base_std, [size(raw_tfs, 1) 1 size(raw_tfs, 3)]);
        
    end

end