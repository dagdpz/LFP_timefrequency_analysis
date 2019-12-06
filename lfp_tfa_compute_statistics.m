function [avg, error] = lfp_tfa_compute_statistics(lfp_tfa_data, error_measure)

if nargin < 2  || isempty(error_measure)
    error_measure = 'bootci';
end 

avg = nanmean(lfp_tfa_data, 1);

error = zeros(2, size(avg, 2));

if size(lfp_tfa_data, 1) > 1

    if strcmp(error_measure, 'stddev')
        stddev = nanstd(lfp_tfa_data, 0, 1);
        error = [avg + stddev; ...
            avg - stddev];
    elseif strcmp(error_measure, 'stderr')
        stderr = nanstd(lfp_tfa_data, 0, 1) / size(lfp_tfa_data, 1);
        error = [avg + stderr; ...
            avg - stderr];
    elseif strcmp(error_measure, 'bootci')

        fn_mean = @(x) nanmean(x, 1);
        %rng(random_seed);
        error = bootci(1000, fn_mean, lfp_tfa_data);
    end
end
    