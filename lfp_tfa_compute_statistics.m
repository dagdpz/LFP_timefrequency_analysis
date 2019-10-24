function [avg, error] = lfp_tfa_compute_statistics(lfp_tfa_data, lfp_tfa_cfg)

if nargin > 2 && isfield(lfp_tfa_cfg, 'error_measure')
    error_measure = lfp_tfa_cfg.error_measure;
else
    error_measure = 'bootci';
end 

if strcmp(error_measure, 'stddev')
    avg = nanmean(lfp_tfa_data, 1);
    stddev = nanstd(lfp_tfa_data, 0, 1);
    error = [avg + stddev; ...
        avg - stddev];
elseif strcmp(error_measure, 'stderr')
    avg = nanmean(lfp_tfa_data, 1);
    stderr = nanstd(lfp_tfa_data, 0, 1) / size(lfp_tfa_data, 1);
    error = [avg + stderr; ...
        avg - stderr];
elseif strcmp(error_measure, 'bootci')
    fn_mean = @(x) nanmean(x, 1);
    rng(lfp_tfa_cfg.random_seed);
    avg = mean(bootstrp(1000, fn_mean, lfp_tfa_data), 1);
    rng(lfp_tfa_cfg.random_seed);
    error = bootci(1000, fn_mean, lfp_tfa_data);
end
    
    