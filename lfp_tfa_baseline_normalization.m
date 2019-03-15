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
%            method     : basline normalization method ('subtraction', 'division', 'zscore', 'relchange')%%- state whose onset to which 
%            mean       : mean power in each frequency bin during the baseline period
%            std        : stddev of baseline power in each frequency bin
%                         baseline.period - start and end of baseline
%                         period wrt refstate onset (seconds)
%                         baseline.method - method for baseline
%                         normalization (not implemented now)
% OUTPUTS:
%		lfp_tfr_normalized  - baseline normalized LFP power spectrogram
%
% See also lfp_tfa_compute_plot_tfr

    base_mean = cfg_baseline.mean;
    if strcmp(cfg_baseline.method , 'zscore')
       base_std =  cfg_baseline.std;
    end
    if strcmp(cfg_baseline.method , 'subtraction')
        lfp_tfr_normalized = raw_tfs - repmat(base_mean, ...
            [1 1 size(raw_tfs, 3)]);
    elseif strcmp(cfg_baseline.method , 'division')
        lfp_tfr_normalized = raw_tfs ./ repmat(base_mean, ...
            [1 1 size(raw_tfs, 3)]);
    elseif strcmp(cfg_baseline.method , 'relchange')
        lfp_tfr_normalized = (raw_tfs - repmat(base_mean, ...
            [1 1 size(raw_tfs, 3)])) ./ repmat(base_mean, ...
            [1 1 size(raw_tfs, 3)]);
    elseif strcmp(cfg_baseline.method , 'zscore')
        lfp_tfr_normalized = (raw_tfs - repmat(base_mean, ...
            [1 1 size(raw_tfs, 3)])) ./ repmat(base_std, ...
            [1 1 size(raw_tfs, 3)]);
    end

end