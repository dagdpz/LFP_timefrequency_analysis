function [ cfg_conditions ] = lfp_tfa_get_trial_conditions( states_lfp, lfp_tfa_cfg )
%UNTITLED6 Summary of this function goes here
%   Detailed explanation goes here
    
    types = lfp_tfa_cfg.compare.types;
    effectors = lfp_tfa_cfg.compare.effectors;
    targets = lfp_tfa_cfg.compare.targets;
    %unique([states_lfp.recorded_hemispace]);
    choice = lfp_tfa_cfg.analyze_choice_trials;
    %unique([states_lfp(1).trials.choice_trial]);
    perturbation = unique([states_lfp(1).trials.perturbation]);
    blocks = unique([states_lfp(1).trials.block]);
    
    npreinj_blocks = sum(perturbation == 0);
    npostinj_blocks = sum(perturbation ~= 0);                     
    
    % create conditions
    cfg_conditions = struct();
    
    i = 0;
%     for rec_hem = recorded_hemispace        
%         for c = choice
%             for b = blocks
%                 i = i + 1;
%                 cfg_conditions(i).recorded_hemispace = rec_hem;
%                 cfg_conditions(i).choice = c;
%                 cfg_conditions(i).block = b;
%                 cfg_conditions(i).perturbation = ...
%                     unique([states_lfp(1).trials([states_lfp(1).trials.block] == b).perturbation]);
%                 cond_label = [];
%                 if cfg_conditions(i).recorded_hemispace == 'L'
%                     cond_label = [cond_label 'Left_hemisphere_'];
%                 else
%                     cond_label = [cond_label 'Right_hemisphere_'];
%                 end
%                 if cfg_conditions(i).choice == 0
%                     cond_label = [cond_label 'Instructed_'];
%                 else
%                     cond_label = [cond_label 'Choice_'];
%                 end
%                 if cfg_conditions(i).perturbation == 0
%                     cond_label = [cond_label 'Control_'];
%                 else
%                     cond_label = [cond_label 'Inactivation_'];
%                 end
%                 cond_label = [cond_label, 'Block_', num2str(cfg_conditions(i).block)];
%                 cfg_conditions(i).label = cond_label;
%                                                 
%             end
%         end
%     end
    
    
    if isfield(lfp_tfa_cfg, 'add_conditions')
        for c = 1:length(lfp_tfa_cfg.add_conditions)
            if ~isempty(lfp_tfa_cfg.add_conditions(c))
                if ~isempty(lfp_tfa_cfg.add_conditions(c).blocks)
                    for rec_hem = targets        
                        for ch = choice
                            i = i + 1;
                            if strcmp(lfp_tfa_cfg.add_conditions(c).blocks, 'control')
                                cfg_conditions(i).block = blocks(perturbation == 0);
                                cfg_conditions(i).perturbation = 0;
                            elseif strcmp(lfp_tfa_cfg.add_conditions(c).blocks, 'inactivation')
                                cfg_conditions(i).block = blocks(perturbation ~= 0);
                                cfg_conditions(i).perturbation = 1;
                            else
                                cfg_conditions(i).block = blocks(lfp_tfa_cfg.add_conditions(c).blocks);
                                cfg_conditions(i).perturbation = sign(perturbation(blocks == blocks(lfp_tfa_cfg.add_conditions(c).blocks(1))));
                                
                            end                    
                            cfg_conditions(i).choice = ch;
                            if isfield(lfp_tfa_cfg.add_conditions(c), 'perturbation')
                                cfg_conditions(i).perturbation = lfp_tfa_cfg.add_conditions(c).perturbation;
                            end
                            cfg_conditions(i).recorded_hemispace = rec_hem;
                            cond_label = [];
                            if cfg_conditions(i).recorded_hemispace == 'L'
                                cond_label = [cond_label 'Left hemisphere_'];
                            else
                                cond_label = [cond_label 'Right hemisphere_'];
                            end
                            if cfg_conditions(i).choice == 0
                                cond_label = [cond_label 'Instructed Trials_'];
                            else
                                cond_label = [cond_label 'Choice Trials_'];
                            end
                            if cfg_conditions(i).perturbation == 0
                                cond_label = [cond_label 'Pre-Injection_'];
                            else
                                cond_label = [cond_label 'Post-Injection_'];
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

