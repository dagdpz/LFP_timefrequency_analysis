function [ session_ecg ] = lfp_tfa_get_block_Rpeak_times( session_ecg, block_Rpeak, nrblock, plottrials, results_folder )
%lfp_tfa_compute_site_tfr Summary of this function goes here
%   Detailed explanation goes here

    if nargin < 4
        plottrials = 0;
    end

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
    ECG_R2Rt = block_Rpeak.R2R_t;
    ECG_R2Rsamples = block_Rpeak.R2R_sample;
    ECG_R2Rvalid = block_Rpeak.R2R_valid;
    ECG_R2Rvalid_bpm = block_Rpeak.R2R_valid_bpm;
    % consecutive valid idx
    ECG_idx_valid_consec = block_Rpeak.idx_valid_R2R_consec;
        
    % concatenate all trials for this run
    block_ecg_timestamps = []; % to concatenate sample time
%     block_LFP = [];
    trials_time = [];
    
    trials_time = vertcat(session_ecg.trials(...
        trials_idx).trialperiod);
    ts = session_ecg.trials(trials_idx(1)).tsample;
    block_ecg_timestamps = ...
        (0:ts:round(trials_time(end)/ts)*ts);
    block_ecg_validsamples = false(size(block_ecg_timestamps));
    
    consec_start_idx = find(diff([0 ECG_idx_valid_consec]) > 1);
    for i = 1:length(consec_start_idx)
        consec_chunk_start_idx = ECG_idx_valid_consec(consec_start_idx(i)) - 1;
        if i < length(consec_start_idx)
            consec_chunk_end_idx =  ...
                ECG_idx_valid_consec(consec_start_idx(i+1)-1);
        else
            consec_chunk_end_idx = ECG_idx_valid_consec(end);
        end
        block_ecg_validsamples(ECG_R2Rsamples(consec_chunk_start_idx)...
            :ECG_R2Rsamples(consec_chunk_end_idx)) = true;        
    end
    
    ECG_peaksamples = round(ECG_timestamps/ts) + 1;
    %ECG_R2Rsamples = round(ECG_R2Rt/ts) + 1;
%     ECG_midsamples = round(mean([ECG_R2Rsamples(1:end-1); ...
%         ECG_R2Rsamples(2:end)]));
%     ECG_midsamples = [round((ECG_peaksamples(1) + ECG_peaksamples(2))/2), ...
%         ECG_midsamples];
    trials_samples = round(trials_time / ts) + 1;
    
    % ECG spikes based on ECG timestamps
    ECG_spikes = false(size(block_ecg_timestamps));
    ECG_spikes(ECG_peaksamples) = true;
    ECG_b2bt = single(nan(size(block_ecg_timestamps)));
    ECG_b2bt(ECG_R2Rsamples) = ECG_R2Rvalid;
    ECG_bpm = single(nan(size(block_ecg_timestamps)));
    ECG_bpm(ECG_R2Rsamples) = ECG_R2Rvalid_bpm;
    
    % fill missing values
    nanx = isnan(ECG_b2bt);
    t    = 1:numel(ECG_b2bt);
    ECG_b2bt(nanx) = interp1(t(~nanx), ECG_b2bt(~nanx), t(nanx));
    nanx = isnan(ECG_bpm);
    t    = 1:numel(ECG_bpm);
    ECG_bpm(nanx) = interp1(t(~nanx), ECG_bpm(~nanx), t(nanx));
        
    % remove invalid intervals
     ECG_b2bt(~block_ecg_validsamples) = nan;
     ECG_bpm(~block_ecg_validsamples) = nan;
            
    % now divide into trials
    for t = 1:length(trials_idx)
        
        trial_ECG_spikes = ECG_spikes(trials_samples(t, 1):trials_samples(t, 2));
        trial_ECG_bpm = ECG_bpm(trials_samples(t, 1):trials_samples(t, 2));
        trial_ECG_b2btime = ECG_b2bt(trials_samples(t, 1):trials_samples(t, 2));
        trial_ECG_valid = block_ecg_validsamples(...
            trials_samples(t, 1):trials_samples(t, 2));
        session_ecg.trials(trials_idx(t)).ECG_spikes = trial_ECG_spikes;
        session_ecg.trials(trials_idx(t)).nRpeaks = sum(trial_ECG_spikes);
        session_ecg.trials(trials_idx(t)).ECG_bpm = trial_ECG_bpm;
        session_ecg.trials(trials_idx(t)).ECG_b2btime = trial_ECG_b2btime;
        session_ecg.trials(trials_idx(t)).ECG_valid = trial_ECG_valid;
        
        % plot individual trials
        if plottrials
            trial = session_ecg.trials(trials_idx(t));
            h = figure(1); clf; hold on;
            
%             h1 = uicontrol('Position', [5 5 200 40], 'String', 'Continue', ...
%                       'Callback', 'uiresume(gcbf)');
%             h2 = uicontrol('Position', [300 5 200 40], 'String', 'Continue', ...
%                       'Callback', 'plottrials = 0; close;');
            % plot raw ECG
            ax1 = subplot(4,2,[3 4]);
            plot(trial.time, trial.ecg_data);
            box on;
            set(gca, 'Xlim', [trial.time(1), ...
                trial.time(end)]);
            set(gca, 'Ylim', ylim);
            ylabel('ECG amplitude');
%             title([strrep(session_ecg.session, '_', ' '), ...
%                 ', Block ', num2str(nrblock) ...
%                 ', Trial ' num2str(trials_idx(t)) ...
%                 ', Completed = ' num2str(trial.completed)]);
            title('Raw ECG');
            % mark state onsets
            for state = (trial.states)
                if ~isnan(state.onset_t) || ~isempty( state.onset_t)
                    line([state.onset_t state.onset_t], ax1.YLim, ...
                        'Color', 'k', 'LineStyle', '--');
                    text(double(state.onset_t), ax1.YLim(1) + ((ax1.YLim(2) - ax1.YLim(1))*(find([trial.states.id] == state.id)/length(trial.states))), num2str(state.id));
                end
            end
            
            % second axis
%             pos = get(ax1,'position');   % get the position vector
%             pos1=pos(2);              % save the original bottom position
%             pos(2)=pos(2)+0.06; %pos(4)=pos(4)-0.01;  % raise bottom/reduce height->same overall upper position
%             set(ax1,'position',pos)   % and resize first axes
%             pos(2)=pos1; pos(4)=0.00000000001; % reset bottom to original and small height
%             ax1(2)=axes('position',pos,'color','none');
%             trial_period = linspace(trial.trialperiod(1), trial.trialperiod(end), length(trial.time));
%             plot(trial_period, zeros(length(trial_period)), 'k');
%             ylabel('Block timestamp (s)');
%             set(ax1(2), 'Xlim', [trial_period(1), ...
%                 trial_period(end)]);
                                    
            % plot spikes
            subplot(4,2,[5 6]);
            plot(trial.time, trial.ECG_spikes);
            box on;
            set(gca, 'Xlim', [trial.time(1), trial.time(end)]);
            title('ECG R peaks');
            
            % plot ECG R2Rt valid
            subplot(4,2,[7 8]);
            plot(trial.time, trial.ECG_b2btime);
            box on;
            set(gca, 'Xlim', [trial.time(1), trial.time(end)]);
            title('ECG R2R interval')
            ylabel('R2R time (s)');
            xlabel('Time (s)');
            
            % Eye position
            ax_eyepos = subplot(422); box on; hold on;
            title('Eye position');
            xlabel('xpos'); ylabel('ypos');
            plot(real(trial.fix_pos), imag(trial.fix_pos), 's', ...
                'Color', [0.5, 0.5, 0.5], 'MarkerSize', 10, 'MarkerFaceColor', [0.5, 0.5, 0.5]);
            plot(real(trial.eye_pos), imag(trial.eye_pos), 'ro', ...
                'MarkerSize', 10, 'MarkerFaceColor', 'r');
            set(ax_eyepos, 'XLim', [-30 30]);
            set(ax_eyepos, 'YLim', [-30 30]);
            
            % text
            ax_text = subplot(421);
            annotation('textbox', get(ax_text, 'OuterPosition'), ...
                'String', sprintf('Session: %s\nBlock: %s\nTrial: %s\nCompleted: %s\n', ...    
                strrep(session_ecg.session, '_', ' '), ...
                num2str(nrblock), ...
                num2str(trials_idx(t)), ...
                num2str(trial.completed)), 'FontSize', 14, ...
                'EdgeColor', 'none', 'HorizontalAlignment', 'center', ...
                'VerticalAlignment', 'middle');
            delete(ax_text);
            
%             uiwait(gcf);
            
            if t == 1
                disp('Press any key to continue! To abort, press Ctrl+C');
            end
            pause;
            %savefig(h, fullfile(results_folder, 'trials', sprintf('Block_%g_Trial_%g', nrblock, trials_idx(t))));
        end
    end


end    
