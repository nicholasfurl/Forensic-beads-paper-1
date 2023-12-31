%%%%%%%%%%%%%%%%%%start, forensic_beads_prior_sim%%%%%%%%%%%%%%%%%%%%%%
function forensic_beads_simulate_parameters;

%forensic_beads_simulate_parameters: In order to clean up Git repository
%for public consumption I've renamed this code and removed older versions.
%This code allows one to play with prior and extra guilt and split
%parameters and see effects on sequence position 

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
stimuli.split = 0.6;
stimuli.extra_guilts = 0;
stimuli.formula = 'stimuli.prior'; 
% stimuli.formula = 'model_behaviour(seq_start_indices(seq)+claim-2,1)/100'; 
% stimuli.formula = '1 - model_behaviour(seq_start_indices(seq)+claim-2,1)/100'; 
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
suspects = unique(means(:,1),'stable');
contexts = unique(means(:,2),'stable');

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
groupvars = { stimuli.raw(:,cols_to_use(1)) stimuli.raw(:,cols_to_use(2)) stimuli.raw(:,cols_to_use(3)) };   %context, claim, participant
temp = grpstats(stimuli.raw, groupvars,'mean');

%get participant averages
% groupvars = {temp(:,1) temp(:,2)};   %suspect, context
groupvars = {temp(:,11) temp(:,8)};   %context, claim
[means meancis] = grpstats(temp,groupvars,{'mean','meanci'});

figure('Color',[1 1 1]);
contexts = unique(means(:,11),'stable');
claims = unique(means(:,8),'stable');

for context = 1:numel(contexts);
%     for claim = 1:numel(claims);
        
        this_data = means(means(:,11) == contexts(context),13);
        %this_data_ci = meancis(means(:,6) == suspects(suspect),10,1);
        plot(this_data); hold on;
        
    end;    %contexts loop
% end;    %suspects loop
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

%on which indices is the display screen 0 (prior rating prompt so first rating
seq_start_indices = find(stimuli.raw(:,5)==0);

%For each start index, loop through sequence and get model predictions and
%adjustments for model and human
for seq = 1:numel(seq_start_indices);
    
    %initialise model probabilities and adjustments
    %     model_behaviour(seq_start_indices,1) = stimuli.prior*100;     %model assumes 50/50 chance
%     model_behaviour(seq_start_indices,2) = NaN;
    
    
    %Loop through this sequence
    for claim = 1:11;
        
        %get model prediction for every seq position
        pg=[];
        q=stimuli.split;
        
        %get number of guilts (i.e., the number of 1s)
        ng = sum( raw(seq_start_indices(seq)+1:seq_start_indices(seq)+claim-1,8) ) + stimuli.extra_guilts;
        %get number of draws so far
        nd = claim-1;
        
        %assign model probability
        %         model_behaviour(seq_start_indices(seq)+claim-1,1) = (1/(1 +
        %         (q/(1-q))^(nd-2*ng)))*100;    %original formula
        %         new_prior = 1 - model_behaviour(seq_start_indices(seq)+claim-2,1)/100;
        new_prior = eval(stimuli.formula);
        %         new_split = (q/(1-q))+stimuli.split_inc;
        model_behaviour(seq_start_indices(seq)+claim-1,1) = (1/(1 + ((1-new_prior)/new_prior)*(q/(1-q))^(nd-2*ng)))*100;
        %assign model adjustment
        if claim == 1;
            model_behaviour(seq_start_indices(seq)+claim-1,2) =  NaN;
        else;
            model_behaviour(seq_start_indices(seq)+claim-1,2) =  model_behaviour(seq_start_indices(seq)+claim-1,1) -  model_behaviour(seq_start_indices(seq)+claim-2,1);
        end;
        
    end;    %loop through this sequence (claim)
    
end;    %loop through sequences
%%%%%%%%%%%%%%%%%%end, simulate%%%%%%%%%%%%%%%%%%%%%%

