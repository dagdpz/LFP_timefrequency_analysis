function [ cfg_conditions ] = lfp_tfa_sitewise_trial_conditions( site_lfp )
%UNTITLED6 Summary of this function goes here
%   Detailed explanation goes here
    recorded_hemispace = unique([site_lfp.recorded_hemispace]);
    choice = unique([site_lfp.trials.choice_trial]);
    perturbation = unique([site_lfp.trials.perturbation]);
    blocks = unique([site_lfp.trials.block]);
    
    % create conditions
    cfg_conditions = struct();
       
    i = 0;
    for rec_hem = recorded_hemispace        
        for c = choice
            for b = blocks
                i = i + 1;
                cfg_conditions(i).recorded_hemispace = rec_hem;
                cfg_conditions(i).choice = c;
                cfg_conditions(i).block = b;
                cfg_conditions(i).perturbation = ...
                    unique([states_lfp(1).trials([states_lfp(1).trials.block] == b).perturbation]);
                cond_label = [];
                if cfg_conditions(i).recorded_hemispace == 'L'
                    cond_label = [cond_label 'Left_hemisphere_'];
                else
                    cond_label = [cond_label 'Right_hemisphere_'];
                end
                if cfg_conditions(i).choice == 0
                    cond_label = [cond_label 'Instructed_'];
                else
                    cond_label = [cond_label 'Choice_'];
                end
                if cfg_conditions(i).perturbation == 0
                    cond_label = [cond_label 'Control_'];
                else
                    cond_label = [cond_label 'Inactivation_'];
                end
                cond_label = [cond_label, 'Block_', num2str(cfg_conditions(i).block)];
                cfg_conditions(i).label = cond_label;
                                                
            end
        end
    end
    
    
    if isfield(lfp_tfa_cfg, 'add_conditions')
        for c = 1:length(lfp_tfa_cfg.add_conditions)
            if ~isempty(lfp_tfa_cfg.add_conditions(c))
                if ~isempty(lfp_tfa_cfg.add_conditions(c).blocks)
                    for rec_hem = recorded_hemispace        
                        for ch = choice
                            i = i + 1;
                            if strcmp(lfp_tfa_cfg.add_conditions(c).blocks, 'inactivation')
                                cfg_conditions(i).block = blocks(perturbation ~= 0);
                                cfg_conditions(i).perturbation = 1;
                            else
                                cfg_conditions(i).block = lfp_tfa_cfg.add_conditions(c).blocks;
                                cfg_conditions(i).perturbation = perturbation(blocks == lfp_tfa_cfg.add_conditions(c).blocks(1));
                                
                            end                    
                            cfg_conditions(i).choice = ch;
                            if isfield(lfp_tfa_cfg.add_conditions(c), 'perturbation')
                                cfg_conditions(i).perturbation = lfp_tfa_cfg.add_conditions(c).perturbation;
                            end
                            cfg_conditions(i).recorded_hemispace = rec_hem;
                            cond_label = [];
                            if cfg_conditions(i).recorded_hemispace == 'L'
                                cond_label = [cond_label 'Left_hemispace_'];
                            else
                                cond_label = [cond_label 'Right_hemispace_'];
                            end
                            if cfg_conditions(i).choice == 0
                                cond_label = [cond_label 'Instructed_'];
                            else
                                cond_label = [cond_label 'Choice_'];
                            end
                            if cfg_conditions(i).perturbation == 0
                                cond_label = [cond_label 'Control_'];
                            else
                                cond_label = [cond_label 'Inactivation_'];
                            end
                            cond_label = [cond_label, 'Block_', num2str(cfg_conditions(i).block)];
                            cfg_conditions(i).label = cond_label;
                            
                        end
                    end
                end
            end
        end
    end

end

