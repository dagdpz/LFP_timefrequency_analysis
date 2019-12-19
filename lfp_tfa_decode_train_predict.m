function lfp_decode = lfp_tfa_decode_train_predict( lfp_decode, lfp_tfa_cfg )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

n_cvfolds = lfp_tfa_cfg.decode.n_cvfolds;
data_labels = [];
%rng(lfp_tfa_cfg.random_seed);
for cl = 1:length(lfp_decode.classes)
    data_labels = [data_labels; ones(lfp_decode.classes(cl).ntrials, 1) * cl];
end

analyses = {'raw_lfp', 'lfp_tfs'};
colors = ['b', 'r', 'g', 'y'];
for i = 1:length(analyses)
    fprintf('%s\n', analyses{i});
    decode_classes_data = [lfp_decode.classes.(analyses{i})];
    class_trials = cat(1, decode_classes_data.trial);

    all_trials = cell(1, size(class_trials, 2));
    all_trials_pc = cell(1, size(class_trials, 2));
    all_trials_concat = [];
    all_trials_concat_pc = [];
    
    train_accuracy_epoch = cell(1, size(class_trials, 2));
    pred_accuracy_epoch = cell(1, size(class_trials, 2));
    train_accuracy_concat = cell(1, size(class_trials, 2));
    pred_accuracy_concat = cell(1, size(class_trials, 2));

    for ep = 1:size(class_trials, 2)        
               
        all_trials{ep} = cat(1, class_trials{:, ep});
        if strcmp(analyses{i}, 'lfp_tfs')
            all_trials{ep} = reshape(all_trials{ep}, size(all_trials{ep}, 1), []);
        end
        all_trials_concat = cat(2, all_trials_concat, all_trials{ep});
        
%         [~, score] = pca(all_trials{ep}, 'NumComponents', 100);
%         all_trials_pc{ep} = score;        
%         
%         [~, score] = pca(all_trials_concat, 'NumComponents', 100);
%         all_trials_concat_pc = score;

        all_trials_pc{ep} = all_trials{ep};
        all_trials_concat_pc = all_trials_concat;
        
        fprintf('Epoch: %s\n', lfp_tfa_cfg.analyse_epochs{ep, 2});
        
        % plot data for classes
        figure(100);
        subplot(length(analyses), size(class_trials, 2), (i-1)*size(class_trials, 2) + ep)
        hold on
        for cl = 1:size(class_trials, 1)
            plot(all_trials_pc{ep}(data_labels==cl,:)', 'color', colors(cl));
        end
        xlim([1, size(all_trials_pc{ep}, 2)])
        
        train_accuracy_epoch{ep} = zeros(1, n_cvfolds);
        pred_accuracy_epoch{ep} = zeros(1, n_cvfolds);
        train_accuracy_concat{ep} = zeros(1, n_cvfolds);
        pred_accuracy_concat{ep} = zeros(1, n_cvfolds);
        
        for c = 1:n_cvfolds
            fprintf('CV fold: %g\n', c);
            cv = cvpartition(length(data_labels),'HoldOut',0.5);
            idx = cv.test;
            test_idx = idx;
            train_idx = ~idx;
        
            train_data = all_trials_pc{ep}(train_idx, :);
            train_labels = data_labels(train_idx);
            test_data = all_trials_pc{ep}(test_idx, :);
            test_labels = data_labels(test_idx);       


            SVMModel = fitcsvm(train_data,train_labels,'KernelFunction','linear',...
                'Standardize',true); % ,'ClassNames',{'negClass','posClass'}       

            [label,~] = predict(SVMModel,train_data);
            accuracy = (sum(label==train_labels) / length(train_labels));
            train_accuracy_epoch{ep}(c) = accuracy;
            fprintf('Train score epoch: %g\n', accuracy);

            [label,~] = predict(SVMModel,test_data);
            accuracy = sum(label==test_labels) / length(test_labels);
            pred_accuracy_epoch{ep}(c) = accuracy;
            fprintf('Test score epoch: %g\n', accuracy);
        
        
        
            % prediction using complete data
            train_data = all_trials_concat_pc(train_idx, :);
            train_labels = data_labels(train_idx);
            test_data = all_trials_concat_pc(test_idx, :);
            test_labels = data_labels(test_idx);

            SVMModel = fitcsvm(train_data,train_labels,'KernelFunction','linear',...
                'Standardize',true); % ,'ClassNames',{'negClass','posClass'}

            [label,~] = predict(SVMModel,train_data);
            accuracy = (sum(label==train_labels) / length(train_labels));
            train_accuracy_concat{ep}(c) = accuracy;
            fprintf('Train score concat: %g\n', accuracy);

            [label,~] = predict(SVMModel,test_data);
            accuracy = (sum(label==test_labels) / length(test_labels));
            pred_accuracy_concat{ep}(c) = accuracy;
            fprintf('Test score concat: %g\n', accuracy);
        end

    end  
    
    figure(200+i);
    hold on;
    stem(cellfun(@mean, train_accuracy_epoch), 'LineWidth', 2);
    stem(cellfun(@mean, pred_accuracy_epoch), 'LineWidth', 2);
    plot(cellfun(@mean, train_accuracy_concat), '-o', 'LineWidth', 2);
    plot(cellfun(@mean, pred_accuracy_concat), '-o', 'LineWidth', 2);
    legend('Train score - epoch', 'Test score - epoch', 'Train score - concatenated', 'Test score - concatenated') 
    title(sprintf('%s', analyses{i}));

end

