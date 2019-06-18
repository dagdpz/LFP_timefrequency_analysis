function condition_label = lfp_tfa_get_condition_label(cfg_condition, label_length)

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

switch cfg_condition.choice
    case 1
        ch='Choice';
        chs='Choi';
    case 0
        ch='Instructed';
        chs='Instr';
    case 'diff1'
        ch='(Instructed - Choice)';
        chs = '(Instr - Choi)';      
end 

switch cfg_condition.perturbation
    case 0
        pert='Pre-injection';
        perts='Pre';
    case 1
        pert='Post-injection';
        perts='Post';
    case 'diff1'
        pert='(Post-injection - Pre-injection)';
        perts = '(Post - Pre)';      
        
end

condition_label=[typ, ' ', eff, ' ', ch, ' ', pert];
condition_label_short=[typs, effs, ' ', chs, ' ', perts];

if strcmp(label_length, 'short')
    condition_label = condition_label_short;
end    

end