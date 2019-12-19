function [ timefreqbinned_tfs, bin_timestamps, bin_freq ] = ...
    lfp_tfa_decode_resample_timefreqbins( lfp_tfs, orig_timepoints, orig_freqpoints, nsamples_timebin, nsamples_freqbin)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

% find number of bins possible
ntimebins = floor(length(orig_timepoints) / nsamples_timebin);
nfreqbins = floor(length(orig_freqpoints) / nsamples_freqbin);

% find the bin timestamps
timebins_idx = round(nsamples_timebin/2):nsamples_timebin:length(orig_timepoints)-round(nsamples_timebin/2);
bin_timestamps = orig_timepoints(timebins_idx);

% resampling by averaging
timebinned_tfs = zeros(size(lfp_tfs, 1), size(lfp_tfs, 2), ntimebins);
for k = 1:ntimebins
    bin_samples_idx = (k-1)*nsamples_timebin + 1:k*nsamples_timebin;
    timebinned_tfs(:,:,k) = nanmean(lfp_tfs(:, :, bin_samples_idx), 3);
end

% find the bin frequencies
freqbins_idx = round(nsamples_freqbin/2):nsamples_freqbin:length(orig_freqpoints)-round(nsamples_freqbin/2);
bin_freq = orig_freqpoints(freqbins_idx);

% resampling by averaging
timefreqbinned_tfs = zeros(size(timebinned_tfs, 1), nfreqbins, size(timebinned_tfs, 3));
for f = 1:nfreqbins
    bin_samples_idx = (f-1)*nsamples_freqbin + 1:f*nsamples_freqbin;
    timefreqbinned_tfs(:,f,:) = nanmean(timebinned_tfs(:, bin_samples_idx, :), 2);
end

