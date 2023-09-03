
%%%%%%%%%%%%%%%%%%start, forensic_beads_prior_sim%%%%%%%%%%%%%%%%%%%%%%
function forensic_beads_prior_sim;

%Gets conditional probabilities biased by subjective prior


%data_trunc.xlsx, I formated from data_exp_11596-v23_task-mwjx (1).csv,
%which Naina acquired from the Gorilla box for this part of the study flowchart.
%1: event index, 2:participant private id, 3:RT, 4: human probability 5: sequence position (which witness is it?)
%6: suspect gender (1=female), 7:witness gender (1=female),
%8: guilty claim (1=guilty), 9: context (mostly innocent / mostly guilty)
stimuli.raw = xlsread('C:\matlab_files\fiance\forensics_beads_MvF_naina\data_trunc.xlsx');

%In raw, sequence positions 0 have NaNs in place of condition labels for contexts (col 9). Put
%them back in or grpstats will exclude the prior from the average
nan_indices = find( stimuli.raw(:,5) == 0);    %find NaNs
stimuli.raw(nan_indices,9) = stimuli.raw(nan_indices+1,9);  %assign the missing values at pos 0 with the values at pos 1

[stimuli.raw(:,[10 11])] = ...
    get_contexts(stimuli);

stimuli.prior = .5;
[stimuli.raw(:,[12 13])] = ...
    get_model_behaviour(stimuli);


%plot probabilities by sequence positions and context

%get participant averages
cols_to_use = [6 9 5 2 12];  %suspect, context, seq pos, participant
groupvars = { stimuli.raw(:,cols_to_use(1)) stimuli.raw(:,cols_to_use(2)) stimuli.raw(:,cols_to_use(3)) stimuli.raw(:,cols_to_use(4))};   %suspect, context, sequence position, participants
temp = grpstats(stimuli.raw(:,cols_to_use),groupvars,'mean');

%get participant averages
groupvars = {temp(:,1) temp(:,2)  temp(:,3)};   %suspect, context, seq pos
[means meancis] = grpstats(temp(:,[1 2 3 5]),groupvars,{'mean','meanci'});

%now get means and ci's over participants

%Loop through and plot probabilities at sequence positions for different suspects
figure('Color',[1 1 1]);
suspects = unique(means(:,1));
contexts = unique(means(:,2));

for suspect = 1:numel(suspects);
    for context = 1:numel(contexts);
        
        this_data = means(means(:,1) == suspects(suspect) & means(:,2) == contexts(context),4);
        %this_data_ci = meancis(means(:,6) == suspects(suspect),10,1);
        plot(this_data); hold on;
        
    end;    %contexts loop
end;    %suspects loop
ylim([1 100]);
box off;

%%%%%%%%%%%%%%%%%
%and plot adjustments by context too

%get participant averages
cols_to_use = [11 8 2 13]; %context (ad hoc), claim, participant, adjustment
groupvars = { stimuli.raw(:,cols_to_use(1)) stimuli.raw(:,cols_to_use(2)) stimuli.raw(:,cols_to_use(3)) stimuli.raw(:,cols_to_use(4))};   %suspect, context, sequence position, participants
temp = grpstats(stimuli.raw(:,cols_to_use),groupvars,'mean');

%get participant averages
groupvars = {temp(:,1) temp(:,2)};   %suspect, context, seq pos
[means meancis] = grpstats(temp(:,[1 2 4]),groupvars,{'mean','meanci'});

figure('Color',[1 1 1]);
contexts = unique(means(:,1));
claims = unique(means(:,2));

for context = 1:numel(contexts);
    for claim = 1:numel(claims);
        
        this_data = means(means(:,1) == contexts(context),3);
        %this_data_ci = meancis(means(:,6) == suspects(suspect),10,1);
        plot(this_data); hold on;
        
    end;    %contexts loop
end;    %suspects loop
ylim([-15 15]);
xlim([0.5 2.5]);
set(gca,'XTick',[1 2]);
xticklabels({'innocent' 'guilty'});
legend({'innocent context' 'guilty context'});
xlabel('Claim');
box off;

fprintf('');





disp('audi5000');
%%%%%%%%%%%%%%%%%%end, forensic_beads_prior_sim%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%start, get_context%%%%%%%%%%%%%%%%%%%%%%
function contexts = get_contexts(stimuli);


%This next loop takes out each sequence, finds the cumulative sum of guilt
%claims for each sequence position and puts it in 9th col of raw, and finds
%cumulative proportion of guilt claims and puts it in 10th column of raw

raw = stimuli.raw;

seq_starts = find(raw(:,5) == 0);   %first element of each sequence



for sequence=1:size(seq_starts,1); %for every start of a sequence
    
    clear this_sequence_claims* temp_MG;
    
    %extract claims for this sequence (positions 1 through 10)
    this_sequence_claims = raw(seq_starts(sequence,1)+1:seq_starts(sequence,1)+10,8);
    %get cumsum of guilt claims (=1 each)
    this_sequence_claims_cumsum = cumsum( this_sequence_claims);
    this_sequence_claims_cumpro = this_sequence_claims_cumsum./[1:10]';
    %save proportions to array called raw
    contexts(seq_starts(sequence,1):seq_starts(sequence,1)+10,1) = [NaN; NaN; this_sequence_claims_cumpro(1:9,:)]; %CONTEXT DEGREE
    
    
    %%%Now do stuff that depends on seq position. Use clunky if/then to
    %%%find if each position's cum proportion is classified as mostly
    %%%guilty (=1) or innocent (=0) or neither (NaN)(raw col 10). Also, find model prob (raw col 11), human adjustment (raw col 12) and
    %%%model adjustment (raw col 13)
    temp_MG = NaN(11,1);    %This will hold MG/MI classifications, which are NaN by default if unclassifiable
    temp_prob_model = NaN(11,1);
    for position = 1:9; %9, because the last position doesn't give a context to any folling claim / rating so should be left off
        
        %Classify in to MI or MG depending on cumulative proportion of preceding guilt claims. After loop assign whole sequence result to raw col 10
        if this_sequence_claims_cumpro(position,1) == .5;  temp_MG(position+2) = NaN;
        elseif this_sequence_claims_cumpro(position,1) > .5; temp_MG(position+2) = 1;
        else; temp_MG(position+2) = 0;
        end;    %If/then to assign MG/MI values
        
    end;    %loop through this seq positions
    
    contexts(seq_starts(sequence,1):seq_starts(sequence,1)+10,2) = temp_MG; %NEW RAW COL 10, CONTEXT CATEGORY
    
end;    %loop seq starts

%%%%%%%%%%%%%%%%%%end, get_context%%%%%%%%%%%%%%%%%%%%%%







%%%%%%%%%%%%%%%%%%start, get_model_behaviour%%%%%%%%%%%%%%%%%%%%%%
function model_behaviour = get_model_behaviour(stimuli)

raw = stimuli.raw;
prior = stimuli.prior;

% %on which indices is the display screen 0 (prior rating prompt so first rating
seq_start_indices = find(stimuli.raw(:,5)==0);

%For each start index, loop through sequence and get model predictions and
%adjustments for model and human
claim_it = 0;
for seq = 1:size(seq_start_indices);
    
    %     %initialise model probabilities and adjustments
    %     model_behaviour(seq_start_indices(seq),1) = stimuli.prior*100;     %model assumes 50/50 chance
    %     model_behaviour(seq_start_indices(seq),2) = NaN;
    
    
    %Loop through this sequence (assume they are always 10 all the time, as they should be)
    ng = 0; %guilt claim counter for this sequence
    for claim = 1:11;
        
%         claim_it = claim_it + 1;    %counter for index into raw
        
        if claim == 1;
            
            claim_it = seq_start_indices(seq);

            model_behaviour(claim_it,1) = stimuli.prior*100;
            model_behaviour(claim_it,2) = NaN;

        else;
            
            claim_it = claim_it + 1;
            
            
            %get model prediction for every seq position
            pg=[];
            q=0.6;
            
            %get number of guilts (i.e., the number of 1s)
            if raw(claim_it,8) == 1;
                ng = ng + 1;
            end;
%             ng = sum( raw(seq_start_indices(seq)+1:seq_start_indices(seq)+claim-1,8) );
            %get number of draws so far
            nd = claim;
            
            %assign model probability
            %         model_behaviour(seq_start_indices(seq)+claim-1,1) = (1/(1 +
            %         (q/(1-q))^(nd-2*ng)))*100;    %original formula
            %         model_behaviour(seq_start_indices(seq)+claim-1,1) = (1/(1 + ((1-prior)/prior)*(q/(1-q))^(nd-2*ng)))*100;
            new_prior = model_behaviour(claim_it - 1,1)/100;
            model_behaviour(claim_it,1) = (1/(1 + ((1-new_prior)/new_prior)*(q/(1-q))^(nd-2*ng)))*100;
            %assign model adjustment
            model_behaviour(claim_it,2) =  model_behaviour(claim_it-1,1) -  model_behaviour(claim_it,1);
            
        end;    %loop through claims for this sequence
        
    end;    %loop through this sequence (claim)
    
end;    %loop through sequences

fprintf('');



%%%%%%%%%%%%%%%%%%end, simulate%%%%%%%%%%%%%%%%%%%%%%





%%%%%%%%%%%%%%%%%%start, fit_data%%%%%%%%%%%%%%%%%%%%%%
function fit_data(stimuli);
stimuli.num_subs = numel(unique(raw(:,2)));

params = 0;

[params_est ...
    ,  ll ...
    , exitflag, search_out] = ...
    fminsearch(  @(params) f_fitparams( params, stimuli ), params);
%%%%%%%%%%%%%%%%%%end, fit_data%%%%%%%%%%%%%%%%%%%%%%








%%%%%%%%%%%%%%%%%%start, f_fitparams%%%%%%%%%%%%%%%%%%%%%%
function ll = f_fitparams( params, stimuli );

probabilities = []; %To hold the model's behavioural data (probability estimates)
ll = 0; %initialise log likelihood

%on which indices is the display screen 0 (prior rating prompt so first rating
seq_start_indices = find(stimuli.raw(:,5)==0);

%For each start index, loop through sequence and get model predictions and
%adjustments for model and human
for seq = 1:numel(seq_start_indices);
    
    %initialise starting values for sequences
    %I've reoganised, now 9: context degree 10: context category 11: model probability 12: adj humans 13: adj model 14: bias (human - model)
    
    probabilities(seq_start_indices,11) = 50;     %model would assume 50/50 chance (starting value), now fitted parameter
    
    
    
    
end;    %loop through seq starts

%%%%%%%%%%%%%%%%%%end, f_fitparams%%%%%%%%%%%%%%%%%%%%%%




num_subs = numel(unique(raw(:,2)));

%on which indices is the display screen 0 (prior rating prompt so first
%rating
seq_start_indices = find(raw(:,5)==0);

%For each start index, loop through sequence and get model predictions and
%adjustments for model and human
for seq = 1:numel(seq_start_indices);
    
    %initialise starting values for sequences
    %I've reoganised, now 9: context degree 10: context category 11: model probability 12: adj humans 13: adj model 14: bias (human - model)
    
    raw(seq_start_indices,11) = 50;     %model assumes 50/50 chance
    raw(seq_start_indices,12) = NaN;    %the first (prior) human's rating cannot be adjusted by anything previous to it
    raw(seq_start_indices,13) = NaN;    %the first (prior) model's rating cannot be adjusted by anything previous to it
    
    %Loop through this sequence
    for claim = 2:11;
        
        %get model prediction for every seq position
        pg=[];
        q=0.6;
        
        %get number of guilts (i.e., the number of 1s)
        ng = sum( raw(seq_start_indices(seq)+1:seq_start_indices(seq)+claim-1,8) );
        %get number of draws so far
        nd = claim-1;
        %
        %         if switch_majority == 1;
        %
        %             %There's no column that codes for the "true" majority in this
        %             %updated code (it's always based on previous draws) so compute
        %             %it manually now
        %             if sum( raw(seq_start_indices(seq)+1:seq_start_indices(seq)+10,8) ) < 5;    %If there are fewer than five guilts
        %
        %                 ng = nd - ng;
        %
        %             end;
        %
        %         end;    %swicth majority
        
        %assign model probability to column 11
        raw(seq_start_indices(seq)+claim-1,11) = (1/(1 + (q/(1-q))^(nd-2*ng)))*100;
        
        
        %         else
        %
        %             %get model prediction for every seq position
        %             %         if raw(seq_start_indices(seq)+1,9)==1; %Figure out if it's majority G or I. Careful of NaNs on the display screen
        %             ng = sum( raw(seq_start_indices(seq)+1:seq_start_indices(seq)+claim-1,8) );
        %             nd = claim-1;
        % %             nd = (claim-1)-ng;  %how many innocents so far?
        %             raw(seq_start_indices(seq)+claim-1,11) = (1/(1 + (q/(1-q))^(nd-2*ng)))*100;
        % %         else
        % %              nd = sum( raw(seq_start_indices(seq)+1:seq_start_indices(seq)+claim-1,8) );
        % %             ng = claim-1;
        % % %             ng = (claim-1)-nd;
        % %             raw(seq_start_indices(seq)+claim-1,10) = (1-(1/(1 + (q/(1-q))^(nd-2*ng))))*100;;
        % %         end;
        
        
        %4: human probs in 4 9: context degree 10: context category 11: model probability;
        %NEW: 12: adj humans 13: adj model 14: bias (human - model)
        raw(seq_start_indices(seq)+claim-1,12) = raw(seq_start_indices(seq)+claim-1,4) - raw(seq_start_indices(seq)+claim-2,4);
        %model adjustments (original model probs in 10, model adjustments in 13)
        raw(seq_start_indices(seq)+claim-1,13) = raw(seq_start_indices(seq)+claim-1,11) - raw(seq_start_indices(seq)+claim-2,11);
        %human - model adjustments (in 14)
        raw(seq_start_indices(seq)+claim-1,14) = raw(seq_start_indices(seq)+claim-1,4) - raw(seq_start_indices(seq)+claim-1,11);
        
    end;    %loop through this sequence (claim)
end;    %loop through seq starts
