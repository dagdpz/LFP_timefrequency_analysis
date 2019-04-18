function [ cmp_conditions ] = lfp_tfa_compare_conditions( states_lfp, lfp_tfa_cfg )
%UNTITLED6 Summary of this function goes here
%   Detailed explanation goes here
    task_types = lfp_tfa_cfg.compare.types;
    effectors = lfp_tfa_cfg.compare.effectors;
    targets = lfp_tfa_cfg.compare.targets;
    %unique([states_lfp.recorded_hemispace]);
    if lfp_tfa_cfg.compare.choice_trials
        choices = unique([states_lfp(1).trials.choice_trial]);
    else
        choices = lfp_tfa_cfg.compare.choice_trials;
    end
    perturbations = lfp_tfa_cfg.compare.perturbations;
    perturbation_groups = lfp_tfa_cfg.compare.perturbation_groups;
%     reach_hands = lfp_tfa_cfg.analyse.reach_hands;
%     reach_spaces = lfp_tfa_cfg.analyse.reach_spaces;              
    
    % create conditions
    cmp_conditions = struct();
    
    i = 0;
    for target = targets
        target_label = target{1};
        for type = task_types
            type_label = ['Type_' num2str(type)];
            for eff = effectors
                eff_label = ['Eff_' num2str(eff)];            
                for ch = choices
                    if ch == 0
                        ch_label = 'Instr';
                    elseif choice == 1
                        ch_label = 'Choice';
                    else
                        ch_label = [];
                    end
                    for p = perturbations
                        if p == 0
                            p_label = 'Pre';
                        elseif p == 1
                            p_label = 'Post';
                        else
                            p_label = [];
                        end  
%                         for hand = reach_hands
%                             if strcmp(hand_label, 'C') || strcmp(hand_label, 'I')
%                                 hand_label = [hand 'H'];
%                             else
%                                 hand_label = [];
%                             end
%                             for s = reach_spaces
%                                 if strcmp(s_label, 'C') || strcmp(s_label, 'I')
%                                     s_label = [s 'S'];
%                                 else
%                                     s_label = [];
%                                 end
                        i = i + 1;
                        condition_label = [type_label, '_', eff_label, '_', ...
                            target_label, '_', ch_label, '_', p_label];
                        cmp_conditions(i).type = type;
                        cmp_conditions(i).effector = eff;
                        cmp_conditions(i).target = target{1};
                        cmp_conditions(i).choice = ch;
                        cmp_conditions(i).perturbation = p;
                        cmp_conditions(i).perturbation_group = perturbation_groups(perturbations == p);
%                         cmp_conditions(i).reach_hand = hand;
%                         cmp_conditions(i).reach_space = s;
                        cmp_conditions(i).label = condition_label;
%                             end
%                         end
                    end
                end
            end
        end
    end
                        
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
    
%     i = 0;
% 
%         
%     
%     
%     if isfield(lfp_tfa_cfg, 'add_conditions')
%         for c = 1:length(lfp_tfa_cfg.add_conditions)
%             if ~isempty(lfp_tfa_cfg.add_conditions(c))
%                 if ~isempty(lfp_tfa_cfg.add_conditions(c).blocks)
%                     for rec_hem = recorded_hemispace        
%                         for ch = choice
%                             i = i + 1;
%                             if strcmp(lfp_tfa_cfg.add_conditions(c).blocks, 'control')
%                                 cmp_conditions(i).block = blocks(perturbation == 0);
%                                 cmp_conditions(i).perturbation = 0;
%                             elseif strcmp(lfp_tfa_cfg.add_conditions(c).blocks, 'inactivation')
%                                 cmp_conditions(i).block = blocks(perturbation ~= 0);
%                                 cmp_conditions(i).perturbation = 1;
%                             else
%                                 cmp_conditions(i).block = blocks(lfp_tfa_cfg.add_conditions(c).blocks);
%                                 cmp_conditions(i).perturbation = sign(perturbation(blocks == blocks(lfp_tfa_cfg.add_conditions(c).blocks(1))));
%                                 
%                             end                    
%                             cmp_conditions(i).choice = ch;
%                             if isfield(lfp_tfa_cfg.add_conditions(c), 'perturbation')
%                                 cmp_conditions(i).perturbation = lfp_tfa_cfg.add_conditions(c).perturbation;
%                             end
%                             cmp_conditions(i).recorded_hemispace = rec_hem;
%                             cond_label = [];
%                             if cmp_conditions(i).recorded_hemispace == 'L'
%                                 cond_label = [cond_label 'Left hemisphere_'];
%                             else
%                                 cond_label = [cond_label 'Right hemisphere_'];
%                             end
%                             if cmp_conditions(i).choice == 0
%                                 cond_label = [cond_label 'Instructed Trials_'];
%                             else
%                                 cond_label = [cond_label 'Choice Trials_'];
%                             end
%                             if cmp_conditions(i).perturbation == 0
%                                 cond_label = [cond_label 'Pre-Injection_'];
%                             else
%                                 cond_label = [cond_label 'Post-Injection_'];
%                             end
%                             cond_label = [cond_label, 'Block_', num2str(cmp_conditions(i).block)];
%                             cmp_conditions(i).label = cond_label;
%                             
%                         end
%                     end
%                 end
%             end
%         end
%     end

end

