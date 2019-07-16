function combined_evoked_lfp = lfp_tfa_get_combined_evoked_lfp( trials_lfp, state )
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

start_state = state{2}(1);
end_state = state{2}(2);
width = state{3};
nwindows = state{4};
spacing = state{5};

combined_evoked_lfp.time = {}; % timestamps
combined_evoked_lfp.lfp = {}; % evoked LFP response

for t = 1:length(trials_lfp)

    states                = trials_lfp(t).states;
    state_start_onset_t   = states([states(:).id] == ...
        start_state).onset_t;
    end_state_onset_t     = states([states(:).id] == ...
        end_state).onset_t;
    
    % get the LFP samples and timestamps between start and end states
    lfp_data = trials_lfp(t).lfp_data(...
        (trials_lfp(t).time >= state_start_onset_t & ...
        trials_lfp(t).time <= end_state_onset_t));
    lfp_time = trials_lfp(t).time(...
        (trials_lfp(t).time >= state_start_onset_t & ...
        trials_lfp(t).time <= end_state_onset_t));
    % lfp sample indices
    lfp_sample_idx = 1:length(lfp_time);
    % sample time
    lfp_ts = 1/trials_lfp(t).fsample;
    
    % number of samples in each window
    nsamples_window = round(width/lfp_ts);
    
    % check if LFP data is sufficiently long to extract the required number
    % of windows of specified length
    if nwindows*(nsamples_window) > length(lfp_time)
        fprintf('LFP is too short to extract the required windows. \n');
        return
    end
    
    % now get the windows to combine
    if strcmp(spacing, 'uniform') % uniformly spaced windows
        window_spacing = round(length(lfp_time)/nwindows);
        window_mid_idx = (window_spacing/2):window_spacing:...
            length(lfp_time) - (window_spacing/2);
        
        % loop through each window
        for w = 1:length(window_mid_idx)
            % evoked LFP for this state
            combined_evoked_lfp.lfp = [combined_evoked_lfp.lfp, ...
                lfp_data(window_mid_idx(w) - round(nsamples_window/2):...
                window_mid_idx(w) + round(nsamples_window/2))];
            % timestamps
            combined_evoked_lfp.lfp_time = ...
                lfp_time(window_mid_idx(w) - round(nsamples_window/2):...
                window_mid_idx(w) + round(nsamples_window/2));
            % set mid-timestamp to zero
            combined_evoked_lfp.lfp_time = combined_evoked_lfp.time - ...
                combined_evoked_lfp.time(round(length(combined_evoked_lfp.time)/2));
        end
        
    end
    
    if strcmp(spacing, 'random') % randomly spaced windows
        window_spacing = round(length(lfp_time)/nwindows);
        window_start_idx = 1:window_spacing:...
            length(lfp_time) - nsamples_window;
        
        % loop through each window
        for w = 1:length(window_start_idx)
            window_sample_idx = lfp_sample_idx(window_start_idx(w) + round(nsamples_window/2):min(window_start_idx(w) + ...
                window_spacing - round(nsamples_window/2) - 1, length(lfp_time)));
            % pick a random sample as window middle
            window_mid_idx = randsample(window_sample_idx, 1);
            %if p_window > rand % a window occurs
                % evoked LFP for this window
                combined_evoked_lfp.lfp = [combined_evoked_lfp.lfp, ...
                    lfp_data(max(1, window_mid_idx - round(nsamples_window/2)):...
                    min(length(lfp_time), window_mid_idx + floor(nsamples_window/2)))];
                % timestamps
                combined_evoked_lfp.lfp_time = ...
                    lfp_time(max(1, window_mid_idx - round(nsamples_window/2)):...
                        min(length(lfp_time), window_mid_idx + round(nsamples_window/2)));
                % set mid-timestamp to zero
                combined_evoked_lfp.lfp_time = combined_evoked_lfp.lfp_time - ...
                    combined_evoked_lfp.lfp_time(round(length(combined_evoked_lfp.lfp_time)/2));
            %end
            %if wnd, continue, end;
        end
    end
end


if ~isempty(combined_evoked_lfp.lfp)

    % crop each lfp to same number of samples
    nsamples = min(cellfun('length', combined_evoked_lfp.lfp));
    for k = 1:length(combined_evoked_lfp.lfp)
        combined_evoked_lfp.lfp{k} = combined_evoked_lfp.lfp{k}(1:nsamples);
    end
    combined_evoked_lfp.lfp_time = combined_evoked_lfp.lfp_time(1:nsamples);
    
    % evoked LFP average
    arr_state_lfp = vertcat(combined_evoked_lfp.lfp{:});
    combined_evoked_lfp.mean = nanmean(arr_state_lfp, 1);
    combined_evoked_lfp.std = nanstd(arr_state_lfp, 0, 1);
    
end

end

