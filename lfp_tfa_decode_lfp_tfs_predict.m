function lfp_decode = lfp_tfa_decode_lfp_tfs_predict( lfp_decode, lfp_tfa_cfg, analyses)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

n_cvfolds = lfp_tfa_cfg.decode.n_cvfolds;
n_tbins_wnd = lfp_tfa_cfg.decode.n_tbins_wnd;

results_folder = fullfile(lfp_tfa_cfg.root_results_fldr, 'LFP Decoding');
if ~exist(results_folder, 'dir')
    try
        mkdir(results_folder);
    catch e
        warning('Cannot create results folder. \nReason: %s\n\n', e.message());
    end
end

colors = 'br';
for a = 1:length(analyses)
    analysis = analyses{a};
    lfp_data = lfp_decode.(analysis);
    % initialize struct to store average across sessions
%     for ep = 1:size(lfp_data.session(1).trial, 1)
%         lfp_data.sessions_avg.train_accuracy{ep} = [];
%         lfp_data.sessions_avg.test_accuracy{ep} = [];
%     end
    % session-wise decoding
    for i = 1:length(lfp_data.session)
        fprintf('Session %g\n', i);
        lfp_data.session(i).train_accuracy = cell(1, size(lfp_data.session(i).trial, 1));
        lfp_data.session(i).test_accuracy = cell(1, size(lfp_data.session(i).trial, 1));
        lfp_data.session(i).timebins = cell(1, size(lfp_data.session(i).trial, 1));
        for ep = 1:size(lfp_data.session(i).trial, 1)
            fprintf('Epoch %s\n', lfp_data.session(i).epoch_name{ep});
            class_trials = lfp_data.session(i).trial(ep, :);
            class_labels = lfp_data.session(i).classes;
            trial_timesamples = lfp_data.session(i).time{ep};
            if isfield(lfp_data.session, 'freq')
                trial_freqsamples = lfp_data.session(i).freq{ep};
            end

            for t = 1:length(class_trials)
                if strcmp(analysis, 'lfp_tfs')
                    [class_trials{t}, trial_timebins, ~] = ...
                        lfp_tfa_decode_resample_timefreqbins(...
                        class_trials{t}, trial_timesamples, trial_freqsamples, ...
                        lfp_tfa_cfg.decode.nsamples_tfs_tbin, ...
                        lfp_tfa_cfg.decode.nsamples_tfs_fbin);
                elseif strcmp(analysis, 'raw_lfp')
                    [class_trials{t}, trial_timebins] = ...
                        lfp_tfa_decode_resample_timebins(...
                        class_trials{t}, trial_timesamples, ...
                        lfp_tfa_cfg.decode.nsamples_lfp_tbin);
                end
            end

            bin_train_accuracy = nan(length(trial_timebins), n_cvfolds);
            bin_test_accuracy = nan(length(trial_timebins), n_cvfolds);

            % loop through each time bin
            % get the valid timebins based on moving window
            if ~n_tbins_wnd
                valid_timebins = 1:length(trial_timebins);
            elseif mod(n_tbins_wnd, 2) %odd
                valid_timebins = max(1, ceil(n_tbins_wnd/2)):length(trial_timebins) - floor(n_tbins_wnd/2);
            elseif ~mod(n_tbins_wnd, 2) %even
                valid_timebins = max(1, round(n_tbins_wnd/2)):length(trial_timebins) - round(n_tbins_wnd/2);
            end

            for b = valid_timebins
                fprintf('Time bin %g\n', b);
                all_trials_concat = [];
                % loop through each trial
                for t = 1:length(class_trials)
                    if mod(n_tbins_wnd, 2) %odd
                        tbins_wnd = (b - floor(n_tbins_wnd/2)):(b + floor(n_tbins_wnd/2));
                    else % even
                        tbins_wnd = (b - floor(n_tbins_wnd/2) + 1):(b + floor(n_tbins_wnd/2));
                    end
                    % get data for all sites for time bin until now
                    if strcmp(analysis, 'lfp_tfs')                    
                        trial_data = reshape(class_trials{t}(:, :, tbins_wnd), 1, []); 
                    elseif strcmp(analysis, 'raw_lfp')
                        trial_data = reshape(class_trials{t}(:, tbins_wnd), 1, []);
                    end
                    all_trials_concat = cat(1, all_trials_concat, trial_data);
                end

                for c = 1:n_cvfolds
                    %fprintf('CV fold: %g\n', c);
                    cv = cvpartition(length(class_labels),'HoldOut',0.5);
                    idx = cv.test;
                    test_idx = idx;
                    train_idx = ~idx;

                    train_data = all_trials_concat(train_idx, :);
                    train_labels = class_labels(train_idx)';
                    test_data = all_trials_concat(test_idx, :);
                    test_labels = class_labels(test_idx)';       


    %                 SVMModel = fitcsvm(train_data,train_labels,'KernelFunction','linear',...
    %                     'Standardize',true);

    %                 [label,~] = predict(SVMModel,train_data);
    %                 accuracy = (sum(label==train_labels') / length(train_labels));
    %                 bin_train_accuracy(b, c) = accuracy;            
    % 
    %                 [label,~] = predict(SVMModel,test_data);
    %                 accuracy = (sum(label==test_labels') / length(test_labels));
    %                 bin_test_accuracy(b, c) = accuracy;

                    SVMModel = svmtrain(train_labels, train_data,'-s 0 -t 0');

                    [~,accuracy,~] = svmpredict(train_labels, train_data, SVMModel);
                    bin_train_accuracy(b, c) = accuracy(1)/100;            

                    [~,accuracy,~] = svmpredict(test_labels, test_data, SVMModel);
                    bin_test_accuracy(b, c) = accuracy(1)/100;

                end

                fprintf('Mean train score concat: %g\n', nanmean(bin_train_accuracy(b,:)));
                fprintf('Mean test score concat: %g\n', nanmean(bin_test_accuracy(b,:)));

            end

            lfp_data.session(i).train_accuracy{ep} = bin_train_accuracy;
            lfp_data.session(i).test_accuracy{ep} = bin_test_accuracy;
            lfp_data.session(i).timebins{ep} = trial_timebins;
%             lfp_data.sessions_avg.train_accuracy{ep} = cat(2, ...
%                 lfp_data.sessions_avg.train_accuracy{ep}, nanmean(bin_train_accuracy, 2));
%             lfp_data.sessions_avg.test_accuracy{ep} = cat(2, ...
%                 lfp_data.sessions_avg.test_accuracy{ep}, nanmean(bin_test_accuracy, 2));
%             if i == 1
%                 lfp_data.sessions_avg.timebins{ep} = trial_timebins;
%             end
            
        end

        figtitle = sprintf('Session %g - %s vs. %s (ntrials = %g, nfold = %g)', i, ...
            lfp_tfa_cfg.decode.classes(1).label, ...
            lfp_tfa_cfg.decode.classes(2).label, ...
            length(test_labels), n_cvfolds);
        h = figure('name', figtitle);
        subplot(1, length(analyses), a);
        hold on;
        set(gca, 'YLim', [0, 1])
        
        timebin_samples = [0];
        timebin_values = [];
        event_onset_samples = zeros(1, length(lfp_data.session(i).timebins));
        wnd_start_samples = zeros(1, length(lfp_data.session(i).timebins));
        wnd_end_samples = zeros(1, length(lfp_data.session(i).timebins));
        
        for ep = 1:length(lfp_data.session(i).timebins)
            event_onset_sample = find(abs(lfp_data.session(i).timebins{ep}) == ...
                min(abs(lfp_data.session(i).timebins{ep})), 1, 'last');
            event_onset_samples(ep) = timebin_samples(end) + event_onset_sample;
            wnd_start_samples(ep) = timebin_samples(end) + 1;
            wnd_end_samples(ep) = timebin_samples(end) + length(lfp_data.session(i).timebins{ep});
            timebin_samples = ...
                timebin_samples(end) + (1:(length(lfp_data.session(i).timebins{ep}) + 1));
            timebin_values = [timebin_values, lfp_data.session(i).timebins{ep}, nan];
            shadedErrorBar(timebin_samples(1:length(lfp_data.session(i).timebins{ep})), ...
                nanmean(lfp_data.session(i).test_accuracy{ep}, 2), ...
                nanstd(lfp_data.session(i).test_accuracy{ep}, 0, 2), colors(a));
            
            line([event_onset_samples(ep) event_onset_samples(ep)], ylim, 'color', 'k');
            text(event_onset_samples(ep) + 0.5, 0.5, lfp_data.session(i).epoch_name{ep});
            
        end
            
%         shadedErrorBar(trial_timebins, nanmean(lfp_data.session(i).test_accuracy, 2), ...
%             nanstd(lfp_data.session(i).test_accuracy, 0, 2), colors(a));
%     %     line([150 150], ylim, 'color', 'k')
%     %     line([280 280], ylim, 'color', 'k')
% 
        title(figtitle);
        xlabel('Time (s)');
        ylabel('Accuracy');
        xticks = [wnd_start_samples; event_onset_samples; wnd_end_samples];
        xticks = unique(xticks(:));
        set(gca, 'xtick', xticks)
        set(gca, 'xticklabels', round(timebin_values(xticks), 1))
        %grid on;
        %set(gca, 'XLim', [timebin_samples(1), timebin_samples(end)])
        
        if a == length(analyses)
            % save figure
            results_file = fullfile(results_folder, figtitle);
            try
                export_fig(h, results_file, '-png');
                export_fig(h, results_file, '-pdf');
            catch e
                warning('Cannot save figures. Reason: %s\n\n', e.message());
            end
        end        

    end
    
    % save back to main struct
    lfp_decode.(analysis) = lfp_data;

end

% save file
try
    save(fullfile(results_folder, 'lfp_decode_data.mat'), 'lfp_decode');
catch e
    warning('Cannot save results. Reason: %s\n\n', e.message());
end

%rng(lfp_tfa_cfg.random_seed);
% for cl = 1:length(raw_lfp.classes)
%     data_labels = [data_labels; ones(raw_lfp.classes(cl).ntrials, 1) * cl];
% end
% 
% analyses = {'raw_lfp', 'lfp_tfs'};
% colors = ['b', 'r', 'g', 'y'];
% for i = 1:length(analyses)
%     fprintf('%s\n', analyses{i});
%     decode_classes_data = [raw_lfp.classes.(analyses{i})];
%     class_trials = cat(1, decode_classes_data.trial);
% 
%     all_trials = cell(1, size(class_trials, 2));
%     all_trials_pc = cell(1, size(class_trials, 2));
%     all_trials_concat = [];
%     all_trials_concat_pc = [];
%     
%     train_accuracy_epoch = cell(1, size(class_trials, 2));
%     pred_accuracy_epoch = cell(1, size(class_trials, 2));
%     train_accuracy_concat = cell(1, size(class_trials, 2));
%     pred_accuracy_concat = cell(1, size(class_trials, 2));
% 
%     for ep = 1:size(class_trials, 2)        
%                
%         all_trials{ep} = cat(1, class_trials{:, ep});
%         if strcmp(analyses{i}, 'lfp_tfs')
%             all_trials{ep} = reshape(all_trials{ep}, size(all_trials{ep}, 1), []);
%         end
%         all_trials_concat = cat(2, all_trials_concat, all_trials{ep});
%         
% %         [~, score] = pca(all_trials{ep}, 'NumComponents', 100);
% %         all_trials_pc{ep} = score;        
% %         
% %         [~, score] = pca(all_trials_concat, 'NumComponents', 100);
% %         all_trials_concat_pc = score;
% 
%         all_trials_pc{ep} = all_trials{ep};
%         all_trials_concat_pc = all_trials_concat;
%         
%         fprintf('Epoch: %s\n', lfp_tfa_cfg.analyse_epochs{ep, 2});
%         
%         % plot data for classes
%         figure(100);
%         subplot(length(analyses), size(class_trials, 2), (i-1)*size(class_trials, 2) + ep)
%         hold on
%         for cl = 1:size(class_trials, 1)
%             plot(all_trials_pc{ep}(data_labels==cl,:)', 'color', colors(cl));
%         end
%         xlim([1, size(all_trials_pc{ep}, 2)])
%         
%         train_accuracy_epoch{ep} = zeros(1, n_cvfolds);
%         pred_accuracy_epoch{ep} = zeros(1, n_cvfolds);
%         train_accuracy_concat{ep} = zeros(1, n_cvfolds);
%         pred_accuracy_concat{ep} = zeros(1, n_cvfolds);
%         
%         for c = 1:n_cvfolds
%             fprintf('CV fold: %g\n', c);
%             cv = cvpartition(length(data_labels),'HoldOut',0.5);
%             idx = cv.test;
%             test_idx = idx;
%             train_idx = ~idx;
%         
%             train_data = all_trials_pc{ep}(train_idx, :);
%             train_labels = data_labels(train_idx);
%             test_data = all_trials_pc{ep}(test_idx, :);
%             test_labels = data_labels(test_idx);       
% 
% 
%             SVMModel = fitcsvm(train_data,train_labels,'KernelFunction','linear',...
%                 'Standardize',true); % ,'ClassNames',{'negClass','posClass'}       
% 
%             [label,~] = predict(SVMModel,train_data);
%             accuracy = (sum(label==train_labels) / length(train_labels));
%             train_accuracy_epoch{ep}(c) = accuracy;
%             fprintf('Train score epoch: %g\n', accuracy);
% 
%             [label,~] = predict(SVMModel,test_data);
%             accuracy = sum(label==test_labels) / length(test_labels);
%             pred_accuracy_epoch{ep}(c) = accuracy;
%             fprintf('Test score epoch: %g\n', accuracy);
%         
%         
%         
%             % prediction using complete data
%             train_data = all_trials_concat_pc(train_idx, :);
%             train_labels = data_labels(train_idx);
%             test_data = all_trials_concat_pc(test_idx, :);
%             test_labels = data_labels(test_idx);
% 
%             SVMModel = fitcsvm(train_data,train_labels,'KernelFunction','linear',...
%                 'Standardize',true); % ,'ClassNames',{'negClass','posClass'}
% 
%             [label,~] = predict(SVMModel,train_data);
%             accuracy = (sum(label==train_labels) / length(train_labels));
%             train_accuracy_concat{ep}(c) = accuracy;
%             fprintf('Train score concat: %g\n', accuracy);
% 
%             [label,~] = predict(SVMModel,test_data);
%             accuracy = (sum(label==test_labels) / length(test_labels));
%             pred_accuracy_concat{ep}(c) = accuracy;
%             fprintf('Test score concat: %g\n', accuracy);
%         end
% 
%     end  
%     
%     figure(200+i);
%     hold on;
%     stem(cellfun(@mean, train_accuracy_epoch), 'LineWidth', 2);
%     stem(cellfun(@mean, pred_accuracy_epoch), 'LineWidth', 2);
%     plot(cellfun(@mean, train_accuracy_concat), '-o', 'LineWidth', 2);
%     plot(cellfun(@mean, pred_accuracy_concat), '-o', 'LineWidth', 2);
%     legend('Train score - epoch', 'Test score - epoch', 'Train score - concatenated', 'Test score - concatenated') 
%     title(sprintf('%s', analyses{i}));
% 
% end
% 
% 
% 
