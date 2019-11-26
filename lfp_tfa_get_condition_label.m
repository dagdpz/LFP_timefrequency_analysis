function condition_label = lfp_tfa_get_condition_label(cfg_condition, label_length)

typ='';
typs='';
if isfield(cfg_condition, 'type')
    switch cfg_condition.type
        case 1
            typ='Fixation';
            typs='F';
        case 2
            typ='Visually_guided';
            typs='V';
        case 3
            typ='Memory';
            typs='M';
        case 4
            typ='Delay';
            typs='D';
        case 5
            typ='M2S';
            typs='M2';
        case 6
            typ='S2S';
            typs='S2';
    end
end

eff='';
effs='';
if isfield(cfg_condition, 'effector')
    if ~isnan(cfg_condition.effector)
        switch cfg_condition.effector
            case 0
                eff='saccades';
                effs='sac';
            case 1
                eff='free_gaze_reaches';
                effs='fgr';
            case 2
                eff='joint';
                effs='joi';
            case 3
                eff='diss_saccades';
                effs='dsa';
            case 4
                eff='diss_reaches';
                effs='dre';
            case 5
                eff='sen_reaches';
                effs='sen';
            case 6
                eff='central_fix_reaches';
                effs='cfr';
        end
    end
end

if isfield(cfg_condition, 'choice')
    if isnan(cfg_condition.choice) || isinf(cfg_condition.choice)
        ch = '';
        chs = '';
    elseif cfg_condition.choice == 1
        ch='Choice';
        chs='Choi';
    elseif cfg_condition.choice == 0
        ch='Instructed';
        chs='Instr';
    elseif strcmp(cfg_condition.choice, 'diff1')
        ch='(Instructed - Choice)';
        chs = '(Instr - Choi)';      
    end 
else
    ch = '';
    chs = '';
end

if isfield(cfg_condition, 'perturbation')
    if isinf( cfg_condition.perturbation) || isnan(cfg_condition.perturbation)
        pert = '';
        perts = '';
    elseif cfg_condition.perturbation == 0
        pert='Pre-injection';
        perts='Pre';
    elseif cfg_condition.perturbation ==  1
        pert='Post-injection';
        perts='Post';
    elseif strfind(cfg_condition.perturbation,  'diff')
        pert='(Post-injection - Pre-injection)';
        perts = '(Post - Pre)';    

    end
else
    pert = '';
    perts = '';
end

if isfield(cfg_condition, 'success')
    if isinf( cfg_condition.success) || isnan(cfg_condition.success)
        succ = '';
        succs = '';
    elseif cfg_condition.success == 0
        succ='Unsuccessful';
        succs='nosuccess';
    elseif cfg_condition.success ==  1
        succ='Successful';
        succs='success';
    elseif strcmp(cfg_condition.success,  'diff1')
        succ='(Successful - Unsuccessful)';
        succs = '(success - nosuccess)';    

    end
else
    succ = '';
    succs = '';
end

condition_label=[typ, ' ', eff, ' ', ch, ' ', pert, ' ', succ];
condition_label_short=[typs, effs, ' ', chs, ' ', perts, ' ', succs];

if strcmp(label_length, 'short')
    condition_label = condition_label_short;
end    

end