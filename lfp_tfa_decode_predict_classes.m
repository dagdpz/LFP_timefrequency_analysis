function lfp_decode_accuracy = lfp_tfa_decode_predict_classes( lfp_decode, lfp_decode_cfg, analyses)
%lfp_tfa_decode_predict_classes - function to decode session-wise LFP data
%   lfp_decode = lfp_tfa_decode_predict_classes( lfp_decode, lfp_tfa_cfg, analyses)
%   INPUTS: 
%       lfp_decode - struct containing session-wise LFP raw data and TFS,
%       output of lfp_tfa_decode_get_conditions_lfp
%       lfp_tfa_cfg - struct containing settings for lfp decoding, see
%       settings/LFP_Decoding_Linus_8sessions/lfp_decoding_settings_Instr_control_IH_ISvsCS.m
%       for example
%       Required fields: 
%           decode.n_cvfolds - number of folds for cross-validation
%           decode.n_tbins_wnd - number of timebins to be considered in
%           the moving window sample
%           root_results_fldr - root folder to save resulting variables and
%           figures
%           decode.classes - struct containing information about the classes to decode
%           decode.nsamples_tfs_tbin - number of samples to be
%           considered for time binning the LFP TFS
%           decode.nsamples_tfs_fbin - number of samples to be
%           considered for frequency binning the LFP TFS
%           decode.nsamples_lfp_tbin - number of samples to be
%           considered for time binning the raw LFP
%   OUTPUTS:
%       lfp_decode - same as input struct 'lfp_decode', but with additional
%       fields containing train and test accuracy (session-wise and average
%       across sessions)
%
%   Required: lfp_tfa_decode_get_class_condition, cvpartition (Statistics
%   and Machine Learning Toolbox), svmtrain and svmpredict (LIBSVM toolbox
%   v3.24), lfp_tfa_decode_plot_accuracy
 
close all;

n_cvfolds = lfp_decode_cfg.decode.n_cvfolds;
wnd_type = lfp_decode_cfg.decode.window_type;
n_tbins_wnd = lfp_decode_cfg.decode.n_tbins_wnd;

results_folder = fullfile(lfp_decode_cfg.root_results_fldr, 'LFP Decoding');
if ~exist(results_folder, 'dir')
    try
        mkdir(results_folder);
    catch e
        warning('Cannot create results folder. \nReason: %s\n\n', e.message());
    end
end

lfp_decode_accuracy = struct();

for a = 1:length(analyses)
    analysis = analyses{a};
    lfp_decode_accuracy(a).analysis = analysis;
    lfp_data = lfp_decode.(analysis);
    lfp_decode_accuracy(a).sessions_avg = struct();
    nepochs = size(lfp_data.session(1).trial, 1);
    for ep = 1:nepochs
        lfp_decode_accuracy(a).sessions_avg.train_accuracy{ep} = [];
        lfp_decode_accuracy(a).sessions_avg.test_accuracy{ep} = [];
    end
    % session-wise decoding
    for i = 1:length(lfp_data.session)
        fprintf('Session %g\n', i);
        % initialize class for each trial
        lfp_data.session(i).class_idx = nan(1, size(lfp_data.session(i).trial, 2));
        classes_target_idx = cell(1, length(lfp_decode_cfg.decode.classes));
        class_trials_idx = cell(1, length(lfp_decode_cfg.decode.classes));
        % assign a class id for each trial
        lfp_decode_accuracy(a).session(i).class_trials_idx = cell(1, ...
            length(lfp_decode_cfg.decode.classes));
        lfp_decode_accuracy(a).session(i).used_sites_idx = cell(1, ...
            length(lfp_decode_cfg.decode.classes));
        for cl = 1:length(lfp_decode_cfg.decode.classes)
            class_condition_idx = lfp_tfa_decode_get_class_condition(...
                lfp_data.session(i).conditions, lfp_decode_cfg.decode.classes(cl));
            class_trials_idx{cl} = ismember(lfp_data.session(i).condition_idx, ...
                class_condition_idx);
            lfp_data.session(i).class_idx(class_trials_idx{cl}) = cl;
            classes_target_idx{cl} = ismember(lfp_data.session(i).targets, ...
                lfp_decode_cfg.decode.classes(cl).target);
            lfp_decode_accuracy(a).session(i).class_trials_idx{cl} = ...
                find(class_trials_idx{cl});
            lfp_decode_accuracy(a).session(i).used_sites_idx{cl} = ...
                find(classes_target_idx{cl});
        end            
        
        nsites = sum(classes_target_idx{1});
        % initialize structs to store session-wise decoding accuracies        
        lfp_decode_accuracy(a).session(i).train_accuracy = cell(nsites, ...
            size(lfp_data.session(i).trial, 1));
        lfp_decode_accuracy(a).session(i).test_accuracy = cell(nsites, ...
            size(lfp_data.session(i).trial, 1));
        lfp_decode_accuracy(a).session(i).timebins = cell(1, size(...
            lfp_data.session(i).trial, 1));
        lfp_decode_accuracy(a).session(i).epoch_name = cell(1, ...
            size(lfp_data.session(i).trial, 1));
        
        % loop through each epoch to analyse
        for ep = 1:size(lfp_data.session(i).trial, 1)
            fprintf('Epoch %s\n', lfp_data.session(i).epoch_name{ep});       
            lfp_decode_accuracy(a).session(i).epoch_name{ep} = ...
                lfp_data.session(i).epoch_name{ep};
            classes_trials_data = lfp_data.session(i).trial(ep, ...
                ~isnan(lfp_data.session(i).class_idx));
            classes_trials_label = lfp_data.session(i).class_idx(...
                ~isnan(lfp_data.session(i).class_idx));
            trial_timesamples = lfp_data.session(i).time{ep};
            if isfield(lfp_data.session, 'freq')
                trial_freqsamples = lfp_data.session(i).freq{ep};
            end                

            % random shuffling of trials individually for each site
            % move this to settings file
            lfp_decode_cfg.decode.shuffle_site_trials = false;

            if lfp_decode_cfg.decode.shuffle_site_trials
                shuffled_class_trials_data = cell(1, length(lfp_decode_cfg.decode.classes));
                shuffled_class_trials_label = [];
                for cl = 1:length(lfp_decode_cfg.decode.classes)
                    class_trials_idx = find(classes_trials_label == cl);
                    shuffled_class_trials_label = [...
                        shuffled_class_trials_label, repmat(...
                        cl, [1, length(class_trials_idx)])];
                    shuffled_class_trials_data{cl} = ...
                        cell(1, length(class_trials_idx));
                    array_class_trials_data = ...
                        permute(cat(4, ...
                        classes_trials_data{class_trials_idx}), [4 1 2 3]);
                    shuffled_array_class_trials_data = array_class_trials_data;
                    for site = 1:nsites
                        shuffled_class_trials_idx = ...
                            randperm(size(array_class_trials_data, 1));
                        for t = 1:length(shuffled_class_trials_idx)
                            shuffled_array_class_trials_data(t, site, :, :) ...
                                = array_class_trials_data(...
                                shuffled_class_trials_idx(t), site, :, :);
                        end
                    end
                    for t = 1:length(class_trials_idx)
                        shuffled_class_trials_data{cl}{t} = ...
                            permute(shuffled_array_class_trials_data(t, :, :, :), ...
                            [2 3 4 1]);
                    end
                end

                shuffled_class_trials_data = [shuffled_class_trials_data{:}];
                % shuffle the classes
                trials_idx = 1:length(shuffled_class_trials_data);
                shuffled_trials_idx = trials_idx(randperm(length(trials_idx)));
                shuffled_class_trials_data = ...
                    shuffled_class_trials_data(shuffled_trials_idx);
                shuffled_class_trials_label = ...
                    shuffled_class_trials_label(shuffled_trials_idx);

                classes_trials_data = shuffled_class_trials_data;
                classes_trials_label = shuffled_class_trials_label;
            end
            
            % Increment number of sites
            for s = 1:nsites
                sites_trials_data = cell(1, length(classes_trials_data));
                fprintf('Site %g\n', s); 
                for t = 1:length(classes_trials_data)
                    if strcmp(analysis, 'lfp_tfs')                    
                        sites_trials_data{t} = classes_trials_data{t}(1:s, :, :); 
                        if lfp_decode_cfg.decode.timebin_lfp_tfs || ...
                                lfp_decode_cfg.decode.freqbin_lfp_tfs
                            [sites_trials_data{t}, trial_timebins, ~] = ...
                                lfp_tfa_decode_resample_timefreqbins(...
                                sites_trials_data{t}, trial_timesamples, trial_freqsamples, ...
                                lfp_decode_cfg.decode.nsamples_tfs_tbin, ...
                                lfp_decode_cfg.decode.nsamples_tfs_fbin);
                        else
                            trial_timebins = trial_timesamples;
                        end
                    elseif strcmp(analysis, 'raw_lfp')
                        sites_trials_data{t} = classes_trials_data{t}(1:s, :); 
                        if lfp_decode_cfg.decode.timebin_lfp
                            [sites_trials_data{t}, trial_timebins] = ...
                                lfp_tfa_decode_resample_timebins(...
                                sites_trials_data{t}, trial_timesamples, ...
                                lfp_decode_cfg.decode.nsamples_lfp_tbin);
                        else
                            trial_timebins = trial_timesamples;
                        end
                    end
                end

                bin_train_accuracy = nan(length(trial_timebins), n_cvfolds);
                bin_test_accuracy = nan(length(trial_timebins), n_cvfolds);

                % loop through each time bin
                % get the valid timebins based on moving window
                if ~n_tbins_wnd || strcmp(wnd_type, 'growing')
                    valid_timebins = 1:length(trial_timebins);
                else
                    valid_timebins = n_tbins_wnd:length(trial_timebins);
                end

                for b = valid_timebins
                    fprintf('Time bin %g\n', b);
                    all_trials_concat = [];
                    % loop through each trial
                    for t = 1:length(sites_trials_data)
                        if strcmp(wnd_type, 'moving')
                            if ~n_tbins_wnd
                                tbins_wnd = b;
                            else
                                tbins_wnd = (b - n_tbins_wnd + 1):(b);
                            end
                        elseif strcmp(wnd_type, 'growing')
                            tbins_wnd = 1:b;
                        end
                        % get data for all sites for time bin until now
                        if strcmp(analysis, 'lfp_tfs')
                            %trial_data = class_trials{t}(:, :, tbins_wnd);
                            trial_data = nanmean(sites_trials_data{t}(:, :, tbins_wnd), 3);
                            %trial_data = reshape(trial_data, 1, []); 
                        elseif strcmp(analysis, 'raw_lfp')
                            %trial_data = class_trials{t}(:, :, tbins_wnd);
                            trial_data = nanmean(sites_trials_data{t}(:, tbins_wnd), 2);                        
                        end
                        trial_data = reshape(trial_data, 1, []);
                        all_trials_concat = cat(1, all_trials_concat, trial_data);
                    end

                    for c = 1:n_cvfolds
                        fprintf('CV fold: %g\n', c);
                        cv = cvpartition(classes_trials_label,'HoldOut',0.5);
                        idx = cv.test;
                        test_idx = idx;
                        train_idx = ~idx;

                        train_data = all_trials_concat(train_idx, :);
                        train_labels = classes_trials_label(train_idx)';
                        test_data = all_trials_concat(test_idx, :);
                        test_labels = classes_trials_label(test_idx)';

                        % normalization to [0, 1] using train data only
%                         train_min = min(train_data);
%                         train_max = max(train_data);
%                         train_data_norm = (train_data - ...
%                             repmat(train_min, size(train_data, 1), 1)) ./ ...
%                             repmat(train_max - train_min, size(train_data, 1), 1);
%                         test_data_norm = (test_data - ...
%                             repmat(train_min, size(test_data, 1), 1)) ./ ...
%                             repmat(train_max - train_min, size(test_data, 1), 1);
                        [train_data_norm, test_data_norm] = ...
                            lfp_tfa_decode_normalize_data(train_data, test_data, 'minmax');


                        % using fitcsvm and predict from MATLAB Statistics
                        % and Machine Learning toolbox 
        %                 SVMModel = fitcsvm(train_data,train_labels,'KernelFunction','linear',...
        %                     'Standardize',true);

        %                 [label,~] = predict(SVMModel,train_data);
        %                 accuracy = (sum(label==train_labels') / length(train_labels));
        %                 bin_train_accuracy(b, c) = accuracy;            
        % 
        %                 [label,~] = predict(SVMModel,test_data);
        %                 accuracy = (sum(label==test_labels') / length(test_labels));
        %                 bin_test_accuracy(b, c) = accuracy;

                        % using svmtrain and svmpredict from libsvm library
                        SVMModel = svmtrain(train_labels, train_data_norm,'-s 0 -t 0');

                        [~,accuracy,~] = svmpredict(train_labels, train_data_norm, SVMModel);
                        bin_train_accuracy(b, c) = accuracy(1)/100;            

                        [~,accuracy,~] = svmpredict(test_labels, test_data_norm, SVMModel);
                        bin_test_accuracy(b, c) = accuracy(1)/100;

                    end

                    fprintf('Mean train score concat: %g\n', nanmean(bin_train_accuracy(b,:)));
                    fprintf('Mean test score concat: %g\n', nanmean(bin_test_accuracy(b,:)));

                end

                lfp_decode_accuracy(a).session(i).train_accuracy{s, ep} = bin_train_accuracy;
                lfp_decode_accuracy(a).session(i).test_accuracy{s, ep} = bin_test_accuracy;
                lfp_decode_accuracy(a).session(i).timebins{ep} = trial_timebins;
                if s == nsites
                    lfp_decode_accuracy(a).sessions_avg.train_accuracy{ep} = cat(2, ...
                        lfp_decode_accuracy(a).sessions_avg.train_accuracy{ep}, nanmean(bin_train_accuracy, 2));
                    lfp_decode_accuracy(a).sessions_avg.test_accuracy{ep} = cat(2, ...
                        lfp_decode_accuracy(a).sessions_avg.test_accuracy{ep}, nanmean(bin_test_accuracy, 2));
                    if i == 1
                        lfp_decode_accuracy(a).sessions_avg.timebins{ep} = trial_timebins;
                        lfp_decode_accuracy(a).sessions_avg.epoch_name{ep} = lfp_data.session(i).epoch_name{ep};
                    end
                end

            end           
            
        end
        
        % plot session-wise accuracy
        figtitle = sprintf('Session %g - %s vs. %s (ntrain = %g, nfold = %g)', i, ...
            lfp_decode_cfg.decode.classes(1).label, ...
            lfp_decode_cfg.decode.classes(2).label, ...
            length(train_labels), n_cvfolds);
        lfp_tfa_decode_plot_accuracy(lfp_decode_accuracy(a).session(i), figtitle, results_folder);

    end
    
    % plot average across sessions
    nsessions = size(lfp_data.sessions_avg.train_accuracy{1}, 2);
    figtitle = sprintf('%s vs. %s (nsessions = %g)', ...
        lfp_decode_cfg.decode.classes(1).label, ...
        lfp_decode_cfg.decode.classes(2).label, ...
        nsessions);
    lfp_tfa_decode_plot_accuracy(lfp_decode_accuracy(a).sessions_avg, figtitle, results_folder);

    
    % save back to main struct
    lfp_decode.(analysis) = lfp_data;

end

% save file
try
    save(fullfile(results_folder, 'lfp_decode_accuracy.mat'), 'lfp_decode_accuracy');
catch e
    warning('Cannot save results. Reason: %s\n\n', e.message());
end

