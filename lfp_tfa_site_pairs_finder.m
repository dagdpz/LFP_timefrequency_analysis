function [site_pairs] = lfp_tfa_site_pairs_finder(monkey_session_allsites_path)
% first we load the allsites file
load(monkey_session_allsites_path);

% find the run and block matches in within trials:
match_trialwise = [];
pair_index=0;
for m = 1:length(allsites_lfp)
    match = [];
    for n = 1:length(allsites_lfp)
        if n > m
            blocks_site_1 = [allsites_lfp(m).trials.block];
            blocks_site_2 = [allsites_lfp(n).trials.block];
            trials_site_1 = [allsites_lfp(m).trials.n];
            trials_site_2 = [allsites_lfp(n).trials.n];
            runs_site_1 = [allsites_lfp(m).trials.run];
            runs_site_2 = [allsites_lfp(n).trials.run];
            
            block_run_site_1=[blocks_site_1; runs_site_1; trials_site_1]'; %n trials x 3
            block_run_site_2=[blocks_site_2; runs_site_2; trials_site_2]';
            
            overlap=intersect(block_run_site_1,block_run_site_2,'rows');
            
            
            if isempty(overlap)
                continue
            end
            
            [match_index_site_1,found_overlap_trials_in_this_position_site_1] = ismember(block_run_site_1,overlap,'rows');
            [match_index_site_2,found_overlap_trials_in_this_position_site_2] = ismember(block_run_site_2,overlap,'rows');
            
            %match_index_site_2 = ismember(block_run_site_2,block_run_site_1,'rows'); % are these always in the same order ! (?)
            
            
            %            for k = length(allsites_lfp(m).trials) % so, for each trial(k) in each site(m), it will compare the same trial with a given different site(n)
            %                if ismember(allsites_lfp(m).trials(k).run, [allsites_lfp(n).trials(k).run]) && ismember(allsites_lfp(m).trials(k).block, allsites_lfp(n).trials(k).block)
            %                   match_trialwise = [match_trialwise; 1]; % -> if the run and block values for the given trials match, this will be registeredas a "1".
            %                end
            %            end
            %            if length(match_trialwise) == length(allsites_lfp.trials)
            %                % if the amount of 1s registered in the previous for loop is the same as the number of trials,
            %                % in other words if there is a match for every trial for a given m,n pair;
            %                match = [match; allsites_lfp(n).site_ID]; % -> a list of all site pairs for every m
            %            end
            %            match_trialwise = [];
            
            
            pair_index=pair_index+1;
            
            
            
            site_pair(pair_index).site1 = allsites_lfp(m).site_ID;
            site_pair(pair_index).site2 = allsites_lfp(n).site_ID;
            site_pair(pair_index).site1_idx = m;
            site_pair(pair_index).site2_idx = n;
            site_pair(pair_index).site1_trials = found_overlap_trials_in_this_position_site_1(match_index_site_1);
            site_pair(pair_index).site2_trials = found_overlap_trials_in_this_position_site_2(match_index_site_2);
            
        end
    end
    
    %site_pair(m).pairs = match;
end

site_pairs = site_pair;
