function [ site_lfp ] = lfp_tfa_get_ECG_peaks( site_lfp, block_ECG )
%lfp_tfa_compute_site_tfr Summary of this function goes here
%   Detailed explanation goes here

%loop through each run
for b = (unique([site_lfp.trials.block]))
    fprintf('Extracting ECG for block %g\n-----------------------\n', b);
    % concatenate all trials for this run
    concat_time = []; % to concatenate sample time
    block_LFP = [];
    trials_time = [];
    trials_idx = find([site_lfp.trials.block] == b);
    for t = 1:length(trials_idx) % ignore first trial
        
        ts = site_lfp.trials(trials_idx(t)).tsample;
        if ~isempty(concat_time) 
            trial_timestamps = ...
                linspace(ts, ts*length(site_lfp.trials(trials_idx(t)).time), ...
                length(site_lfp.trials(trials_idx(t)).time)) + concat_time(end);
        else
            trial_timestamps = ...
                linspace(0, ts*(length(site_lfp.trials(trials_idx(t)).time)-1), ...
                length(site_lfp.trials(trials_idx(t)).time));
        end
        
        
%         trial_timestamps = ...
%             linspace(0, (length(site_lfp.trials(trials_idx(t)).time) - 1)*ts, ...
%                 length(site_lfp.trials(trials_idx(t)).time));

        concat_time = [concat_time trial_timestamps];
        trials_time = [trials_time; [trial_timestamps(1), trial_timestamps(end)]];
        
        % initialize ECG spikes
        %site_lfp.trials(trials_idx(t)).ECG_spikes = false(size(trial_timestamps));
    end
    
    % get ECG timestamps for this block
    ECG_timestamps = block_ECG.out(b).Rpeak_t;
    ECG_R2Rt = block_ECG.out(b).R2R_t;
    ECG_R2Rvalid = block_ECG.out(b).R2R_valid;
    ECG_R2Rvalid_bpm = block_ECG.out(b).R2R_valid_bpm;
    
    ECG_peaksamples = round(ECG_timestamps/ts) + 1;
    ECG_R2Rsamples = round(ECG_R2Rt/ts) + 1;
    trials_samples = round(trials_time / ts) + 1;
    
    % ECG spikes based on ECG timestamps
    ECG_spikes = false(size(concat_time));
    ECG_spikes(ECG_peaksamples) = true;
    ECG_b2bt = single(nan(size(concat_time)));
    ECG_b2bt(ECG_R2Rsamples) = ECG_R2Rvalid;
    ECG_bpm = single(nan(size(concat_time)));
    ECG_bpm(ECG_R2Rsamples) = ECG_R2Rvalid_bpm;
    
    % fill missing values
    nanx = isnan(ECG_b2bt);
    t    = 1:numel(ECG_b2bt);
    ECG_b2bt(nanx) = interp1(t(~nanx), ECG_b2bt(~nanx), t(nanx));
    nanx = isnan(ECG_bpm);
    t    = 1:numel(ECG_bpm);
    ECG_bpm(nanx) = interp1(t(~nanx), ECG_bpm(~nanx), t(nanx));
               
    % now divide into trials
    for t = 1:length(trials_idx)
        
        trial_ECG_spikes = ECG_spikes(...
            concat_time >= trials_time(trials_idx == trials_idx(t), 1) & ...
            concat_time <= trials_time(trials_idx == trials_idx(t), 2));
        trial_ECG_bpm = ECG_bpm(...
            concat_time >= trials_time(trials_idx == trials_idx(t), 1) & ...
            concat_time <= trials_time(trials_idx == trials_idx(t), 2));
        trial_ECG_b2btime = ECG_b2bt(...
            concat_time >= trials_time(trials_idx == trials_idx(t), 1) & ...
            concat_time <= trials_time(trials_idx == trials_idx(t), 2));
        site_lfp.trials(trials_idx(t)).ECG_spikes = trial_ECG_spikes;
        site_lfp.trials(trials_idx(t)).ECG_bpm = trial_ECG_bpm;
        site_lfp.trials(trials_idx(t)).ECG_b2btime = trial_ECG_b2btime;
    end
        
end



end
