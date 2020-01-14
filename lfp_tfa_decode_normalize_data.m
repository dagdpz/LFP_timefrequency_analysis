function [ train_data_norm, test_data_norm ] = lfp_tfa_decode_normalize_data( train_data, test_data, normalization )
%UNTITLED14 Summary of this function goes here
%   Detailed explanation goes here

switch(normalization)
    case 'minmax'
        train_min = min(train_data);
        train_max = max(train_data);
        train_data_norm = (train_data - ...
            repmat(train_min, size(train_data, 1), 1)) ./ ...
            repmat(train_max - train_min, size(train_data, 1), 1);
        test_data_norm = (test_data - ...
            repmat(train_min, size(test_data, 1), 1)) ./ ...
            repmat(train_max - train_min, size(test_data, 1), 1);
    case 'zscore'
        train_mean = mean(train_data);
        train_std = std(train_data);
        train_data_norm = (train_data - ...
            repmat(train_mean, size(train_data, 1), 1)) ./ ...
            repmat(train_std, size(train_data, 1), 1);
        test_data_norm = (test_data - ...
            repmat(train_mean, size(test_data, 1), 1)) ./ ...
            repmat(train_std, size(test_data, 1), 1);
end

end

