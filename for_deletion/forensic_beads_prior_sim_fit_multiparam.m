%%%%%%%%%%%%%%%%%%start, forensic_beads_prior_sim%%%%%%%%%%%%%%%%%%%%%%
function forensic_beads_prior_sim_fit;

%This version doesn't let you configure which parameters to use.

%forensic_beads_prior_sim_fit_multiparam.m expands the version of 
%forensic_beads_prior_sim_fit.m that existed on 2/Sept/2023 to introduce a
%number of guilt claims increment parameter.

%forensic_beads_prior_sim_fit.m converts simulation code to model fitting
%code. Beginning with just fitting of the prior in the two suspect conditions.

%forensic_beads_prior_sim_fit.m - simulates conditional probabilities biased by subjective prior 

addpath(genpath('C:\matlab_files\fiance\forensic_beads_pub_repo\Forensic-beads-paper-1\FMINSEARCHBND'))

%data_trunc.xlsx, I formated from data_exp_11596-v23_task-mwjx (1).csv,
%which Naina acquired from the Gorilla box for this part of the study flowchart.
%1: event index, 
%2:participant private id, 
%3:RT, 
%4: human probability 
%5: sequence position (which witness is it?)
%6: suspect gender (1=female), 
%7:witness gender (1=female),
%8: guilty claim (1=guilty), 
%9: context (mostly innocent / mostly guilty)
stimuli.raw = xlsread('C:\matlab_files\fiance\forensics_beads_MvF_naina\data_trunc.xlsx');

%In raw, sequence positions 0 have NaNs in place of condition labels 
%for contexts (col 9) and sometimes suspects (col 6). Put
%them back in or you'll have troubles later
nan_indices = find( stimuli.raw(:,5) == 0);    %find NaNs
stimuli.raw(nan_indices,9) = stimuli.raw(nan_indices+1,9);  %assign the missing values at pos 0 with the values at pos 1
stimuli.raw(nan_indices,6) = stimuli.raw(nan_indices+1,6);  %assign the missing values at pos 0 with the values at pos 1


[stimuli.raw(:,[10 11])] = ...
    get_contexts(stimuli);

%Who are the participants?
participant_list = unique(stimuli.raw(:,2));
num_participants = numel(participant_list);

%initial value of free params
params(1) = .5; %prior, initialised to optimal value (ground truth of paradigm)
params(2) = 0;  %guilt claim increment, intitialised to optimal value
params(3) = 0;  %bias term, intialised to optimal value
params(4) = 1;  %noise term, initialised to optimal value

%data matrix to hold modelling results, with each participant in a row
model_results = [];

%fit for each condition in each participant and save results to new dataset to analyse
for participant = 1:num_participants;
    
    %get data for this participant
    this_ps_data = stimuli.raw(find(stimuli.raw(:,2) == participant_list(participant)),:);
    
    %get condition codes (all participants should always be 0 and 1 but, hey, I like being unnecessarily robust sometimes)
    this_ps_suspects = unique(rmmissing(this_ps_data(:,7)));
    num_suspects = numel(this_ps_suspects);
    
    %Now loop through the detected conditions
    for suspect = 1:num_suspects;
        
        disp(sprintf('fitting participant %d suspect %d', participant, suspect))
        
        %Get the probability rating data to fit for this suspect in this participant
        %We can use the sequence position (5), the claim (8) and the rating itself (4)
        this_ps_suspect_data = this_ps_data(this_ps_data(:,6) == this_ps_suspects(suspect),[5 8 4]);
        
        %pass data, initialised param and function handle to fminsearch
        [params_est(participant, suspect,:) ll(participant, suspect) flag search] = ...
            fminsearchbnd( ...
            @(params) get_model_ll(params, this_ps_suspect_data), ...
            params, ...
            [0 0 0 0], ...  %lower parameter bounds
            [1 Inf Inf Inf] ... %upper parameter bounds
            );
        
        %Now that this participant has been fit, get model performance
        
        %get performance for this model
        %model_results:
        %col 1: sequence position, col 2: claim, col 3: human rating, col
        %4: seq num, col 5: model rating col 6: participant number
        temp = get_model_behaviour(params_est(participant, suspect,:),this_ps_suspect_data);
        model_results = ... 
            [ ...
            model_results; ...
            [temp ...
            this_ps_suspects(suspect)*ones(size(temp,1),1) ...
            participant_list(participant)*ones(size(temp,1),1) ...
            ] ...
            ];

        
    end;    %loop through suspects
       
    save('test_fit_multiparam.m');
    
end;    %loop through participants

fprintf('');


% %from *_old
% stimuli.prior = .5;
% stimuli.split = .6;
% stimuli.extra_guilts = 1;
% stimuli.formula = 'stimuli.prior'; 
% % stimuli.formula = 'model_behaviour(seq_start_indices(seq)+claim-2,1)/100'; 
% % stimuli.formula = '1 - model_behaviour(seq_start_indices(seq)+claim-2,1)/100'; 
% [stimuli.raw(:,[12 13])] = ...
%     get_model_behaviour(stimuli);


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




%%%%%%%%%%%%%%%%%%start, get_model_ll%%%%%%%%%%%%%%%%%%%%%%
function ll = get_model_ll(params,this_ps_suspect_data);

% prior = params(1);

params_get_behav = params([1 2]);
params_bias = params(3);
params_noise = params(4);

%this ps_suspect_data:
%col 1: sequence position, col 2: claim, col 3: human rating, col 4: seq num, col 5: model rating
this_ps_suspect_data = get_model_behaviour(params_get_behav,this_ps_suspect_data);


%compute ll
%ll = -sum(log(this_ps_suspect_data(:,3)/100.*this_ps_suspect_data(:,5)/100));

ll = 0;
for trial = 1:size(this_ps_suspect_data,1);

    %model response noise and bias when gebnerating predictions
%     y_hat = this_ps_suspect_data(trial,5)/100;
    y_hat = params_bias + params_noise*this_ps_suspect_data(trial,5)/100;
    
    %Get "labels" for this trial
    y = this_ps_suspect_data(trial,3)/100;
    
    likelihood = y_hat^y*(1-y_hat)^(1-y);
    ll = ll - likelihood;

%      ll = ll + y*log(y_hat) + (1-y)*log(1-y_hat);

end;    %loop through trials

% ll = -ll;


fprintf('');
%%%%%%%%%%%%%%%%%%end, get_model_ll%%%%%%%%%%%%%%%%%%%%%%




%%%%%%%%%%%%%%%%%%start, get_model_behaviour%%%%%%%%%%%%%%%%%%%%%%
function this_ps_suspect_data = get_model_behaviour(params, this_ps_suspect_data)

prior = params(1);
guilt_claim_inc = params(2);

%on which indices is the display screen 0 (prior rating prompt so first rating)
seq_start_indices = find(this_ps_suspect_data(:,1)==0);

%For each start index, loop through sequence and get model predictions
for seq = 1:numel(seq_start_indices);
    
    %initialise model probabilities
    
    %     this_ps_suspect_data(seq_start_indices,2) = NaN;
    
    %Loop through this sequence
    for claim = 1:11;
        
        %what's the current index into this_ps_suspect_data?
        index = seq_start_indices(seq)+claim-1;
        
        %save sequence number so can loop more easily later
        this_ps_suspect_data(index,4) = seq;
        
        if claim == 1;
            
            this_ps_suspect_data(seq_start_indices(seq),5) = prior*100;
            
        else
            
            %get model prediction for every seq position
            q=.6;
            
            %get number of guilts (i.e., the number of 1s)
            ng = sum( this_ps_suspect_data(seq_start_indices(seq)+1:index,2) ) + guilt_claim_inc;
            
            %get number of draws so far
            nd = claim-1;
            
            %assign model probability
            this_ps_suspect_data(index,5) = (1/(1 + ((1-prior)/prior)*(q/(1-q))^(nd-2*ng)))*100;
                      
        end;    %before first claim (sequence position 0) or a later one?
        
    end;    %loop through this sequence (claim)
    
end;    %loop through sequences
fprintf('');
%%%%%%%%%%%%%%%%%%end, simulate%%%%%%%%%%%%%%%%%%%%%%






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







