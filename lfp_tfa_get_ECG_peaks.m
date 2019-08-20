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
    ECG_timestamps = block_ECG(b).Rpeak_t;
    % ECG spikes based on ECG timestamps
    ECG_spikes = false(size(concat_time));
    for k = 1:length(ECG_timestamps)
        if min(abs(concat_time - ECG_timestamps(k))) < ts
            time_lgcl_idx = abs(concat_time - ECG_timestamps(k)) == ...
                min(abs(concat_time - ECG_timestamps(k)));
            ECG_spikes(time_lgcl_idx) = true;            
        end
    end
            
    % now divide into trials
    for t = 1:length(trials_idx)
        trial_ECG_spikes = ECG_spikes(...
            concat_time >= trials_time(trials_idx == trials_idx(t), 1) & ...
            concat_time <= trials_time(trials_idx == trials_idx(t), 2));
        site_lfp.trials(trials_idx(t)).ECG_spikes = trial_ECG_spikes;
    end
        
end



end
