function powspctrm_norm = freqBaseLineNorm(ft_TFR, baseline, method)
    powspctrm = ft_TFR.powspctrm;    
    sq_powspctrm = squeeze(ft_TFR.powspctrm);
    powspctrm_norm = zeros(size(powspctrm));
    % find mean power for each freq bin during baseline period
    powspctrm_baseline = sq_powspctrm(:,ft_TFR.time >= baseline(1) & ft_TFR.time <= baseline(2));
    mean_pow_baseline = nanmean(powspctrm_baseline, 2);
    for f = 1:length(ft_TFR.freq)
        powspctrm_norm(1,f,:) = sq_powspctrm(f,:) - mean_pow_baseline(f);
    end
end