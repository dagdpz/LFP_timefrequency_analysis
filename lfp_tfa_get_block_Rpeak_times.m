function [ session_ecg ] = lfp_tfa_get_block_Rpeak_times( session_ecg, block_Rpeak, nrblock )
%lfp_tfa_compute_site_tfr Summary of this function goes here
%   Detailed explanation goes here

    fprintf('Extracting Rpeaks for block %g\n-----------------------\n', nrblock);
    
    % find trials for this block
    trials_idx = find([session_ecg.trials.block] == nrblock);
    if numel(trials_idx) == 0
        fprintf('No ECG data found for block %g\n', nrblock);
        return;
    end
    
    % get ECG timestamps for this block
    ECG_timestamps = block_Rpeak.Rpeak_t;
    if isempty(ECG_timestamps)
        fprintf('No ECG data found for block %g\n', nrblock);
        session_ecg.trials(trials_idx) = [];
        return;
    end
        
    % concatenate all trials for this run
    block_ecg_timestamps = []; % to concatenate sample time
%     block_LFP = [];
    trials_time = [];
    
    trials_time = vertcat(session_ecg.trials(...
        trials_idx).trialperiod);
    ts = session_ecg.trials(trials_idx(1)).tsample;
    block_ecg_timestamps = ...
        (0:ts:round(trials_time(end)/ts)*ts);
    ECG_peaksamples = round(ECG_timestamps/ts) + 1;
    ECG_midsamples = round(mean([ECG_peaksamples(1:end-1); ...
        ECG_peaksamples(2:end)]));
    trials_samples = round(trials_time / ts) + 1;
    
%     for t = 1:length(trials_idx) % ignore first trial
%         
%         ts = session_ecg.trials(trials_idx(t)).tsample;
% %         if ~isempty(concat_time) 
% %             trial_timestamps = ...
% %                 linspace(ts, ts*length(session_ecg.trials(trials_idx(t)).time), ...
% %                 length(session_ecg.trials(trials_idx(t)).time)) + concat_time(end);
% %         else
% %             trial_timestamps = ...
% %                 linspace(0, ts*(length(session_ecg.trials(trials_idx(t)).time)-1), ...
% %                 length(session_ecg.trials(trials_idx(t)).time));
% %         end
%         trial_time = session_ecg.trials(trials_idx(t)).trial_period;        
%         trial_timestamps = trial_time(1):ts:trial_time(end);      
%         
% %         trial_timestamps = ...
% %             linspace(0, (length(site_lfp.trials(trials_idx(t)).time) - 1)*ts, ...
% %                 length(site_lfp.trials(trials_idx(t)).time));
% 
%         block_ecg_timestamps = [block_ecg_timestamps trial_timestamps];
%         trials_time = [trials_time; trial_time];
%         
%         % initialize ECG spikes
%         %site_lfp.trials(trials_idx(t)).ECG_spikes = false(size(trial_timestamps));
%     end
        
    b2btime = diff(ECG_timestamps);
    % detect outliers
    b2btime(b2btime > mean(b2btime) + 3*std(b2btime) | ...
        b2btime < mean(b2btime) - 3*std(b2btime)) = nan;
    beat_idx = 2:length(ECG_timestamps);
    %beat_rate = (beat_idx ./ ECG_timestamps(2:end))*60;
    beat_rate = (1 ./ b2btime)*60;
    % ECG spikes based on ECG timestamps
    ECG_spikes = false(size(block_ecg_timestamps));
    ECG_spikes(ECG_peaksamples) = true;
    ECG_b2bt = single(nan(size(block_ecg_timestamps)));
    ECG_b2bt(ECG_midsamples) = b2btime;
    ECG_bpm = single(nan(size(block_ecg_timestamps)));
    ECG_bpm(ECG_midsamples) = beat_rate;
    % get ECG peak to peak interval and beat rate
    %beat_time_idx = false(1, length(block_ecg_timestamps));
%     for k = 1:length(ECG_timestamps)
%         if min(abs(block_ecg_timestamps - ECG_timestamps(k))) < ts
% %            time_idx = find(abs(concat_time - ECG_timestamps(k)) < ts*0.5);
%             ECG_spikes(abs(block_ecg_timestamps - ...
%                 ECG_timestamps(k)) == min(abs(block_ecg_timestamps - ...
%                 ECG_timestamps(k))))= true;
%             if k < length(ECG_timestamps)
%                 mid_time_idx = find(abs(block_ecg_timestamps - ...
%                     0.5*(ECG_timestamps(k) + ECG_timestamps(k+1))) == ...
%                     min(abs(block_ecg_timestamps - ...
%                     0.5*(ECG_timestamps(k) + ECG_timestamps(k+1)))));
%                 ECG_b2bt(mid_time_idx) = b2btime(k+1);
%                 ECG_bpm(mid_time_idx) = beat_rate(k+1);
%             end
%         end
%         
%     end
    
%     ECG_spikes(beat_time_idx) = true; 
%     ECG_b2bt(beat_time_idx)= b2btime;
%     ECG_bpm(beat_time_idx) = beat_rate;
    
    %beatidx = cumsum(ECG_spike    
    % fill missing values
    nanx = isnan(ECG_b2bt);
    t    = 1:numel(ECG_b2bt);
    ECG_b2bt(nanx) = interp1(t(~nanx), ECG_b2bt(~nanx), t(nanx));
    nanx = isnan(ECG_bpm);
    t    = 1:numel(ECG_bpm);
    ECG_bpm(nanx) = interp1(t(~nanx), ECG_bpm(~nanx), t(nanx));
%     [ECG_b2bt, ~] = fillmissing(ECG_b2bt,'linear','SamplePoints',block_ecg_timestamps);
%     [ECG_bpm, ~] = fillmissing(ECG_bpm,'linear','SamplePoints',block_ecg_timestamps);
    
%     for k = 1:length(b2btime)
%         ECG_beatrate = k/(ECG_timestamps(k));
%         ECG_b2bt(beatidx ==k) = b2btime(k);
%         ECG_bpm(beatidx ==k) = ECG_beatrate*60;
%     end
        
        
            
    % now divide into trials
    for t = 1:length(trials_idx)
        
%         trial_ECG_spikes = ECG_spikes(...
%             block_ecg_timestamps >= trials_time(trials_idx == trials_idx(t), 1) & ...
%             block_ecg_timestamps <= trials_time(trials_idx == trials_idx(t), 2));
%         trial_ECG_bpm = ECG_bpm(...
%             block_ecg_timestamps >= trials_time(trials_idx == trials_idx(t), 1) & ...
%             block_ecg_timestamps <= trials_time(trials_idx == trials_idx(t), 2));
%         trial_ECG_b2btime = ECG_b2bt(...
%             block_ecg_timestamps >= trials_time(trials_idx == trials_idx(t), 1) & ...
%             block_ecg_timestamps <= trials_time(trials_idx == trials_idx(t), 2));
        trial_ECG_spikes = ECG_spikes(trials_samples(t, 1):trials_samples(t, 2));
        trial_ECG_bpm = ECG_bpm(trials_samples(t, 1):trials_samples(t, 2));
        trial_ECG_b2btime = ECG_b2bt(trials_samples(t, 1):trials_samples(t, 2));
        session_ecg.trials(trials_idx(t)).ECG_spikes = trial_ECG_spikes;
        session_ecg.trials(trials_idx(t)).nRpeaks = sum(trial_ECG_spikes);
        session_ecg.trials(trials_idx(t)).ECG_bpm = trial_ECG_bpm;
        session_ecg.trials(trials_idx(t)).ECG_b2btime = trial_ECG_b2btime;
    end


end
