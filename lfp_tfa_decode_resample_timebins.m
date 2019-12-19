function [ resampled_lfp, timebins ] = lfp_tfa_decode_resample_timebins( raw_lfp, orig_timepoints, nsamples_timebin )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

% find number of bins possible
nbins = floor(length(orig_timepoints) / nsamples_timebin);

% find the new timestamps
timebins_idx = round(nsamples_timebin/2):nsamples_timebin:length(orig_timepoints)-round(nsamples_timebin/2);
timebins = orig_timepoints(timebins_idx);

% resampling by averaging
resampled_lfp = zeros(size(raw_lfp, 1), nbins);
for k = 1:nbins
    bin_samples_idx = (k-1)*nsamples_timebin + 1:k*nsamples_timebin;
    resampled_lfp(:,k) = nanmean(raw_lfp(:, bin_samples_idx), 2);
end

