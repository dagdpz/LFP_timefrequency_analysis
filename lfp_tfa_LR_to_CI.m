function pop=lfp_tfa_LR_to_CI(keys,pop)
%% assigning contra and ipsi instead of left/right, dependent on the recorded target (which contains hemisphere information in the end)
% reference is left hemisphere, meaning that right becomes contra and left ipsi
% that is why hemifields, positions and hands need to be inverted for right hemisphere targets.
% for right recording sites LR->CI  (hemifield & positions...?)
% positive space and hand==2 become contra

target=pop.(keys.contra_ipsi_relative_to);
if strcmpi(target(end-1:end),'_R') || strcmpi(target(end-1:end),'_r')
    poptr=num2cell([pop.trials.hemifield]*-1);
    [pop.trials.hemifield]=deal(poptr{:});
    poptr=cellfun(@(x) x.*[-1,1],{pop.trials.position},'uniformoutput',false);
    [pop.trials.position]=deal(poptr{:});
    poptr=cellfun(@(x) x.*[-1,1],{pop.trials.fixation},'uniformoutput',false);
    [pop.trials.fixation]=deal(poptr{:});
    hand1=[pop.trials.reach_hand]==2;
    hand2=[pop.trials.reach_hand]==1;
    [pop.trials(hand1).reach_hand]=deal(1);
    [pop.trials(hand2).reach_hand]=deal(2);
    
    %% fixation index--> temporary solution
    if isfield(pop.trials(1),'fix_index')
        fixindex=[pop.trials.fix_index];
        [pop.trials(fixindex==3).fix_index]=deal(1);
        [pop.trials(fixindex==1).fix_index]=deal(3);
    end
end
end