%%%%%%%%%%%%%%%%%%start, forensic_beads_parameter_recovery%%%%%%%%%%%%%%%%%%%%%%
function forensic_beads_parameter_recovery;

%forensic_beads_parameter_recovery.m simulates data using stimuli from
%study 1 and a range of prior, response bias, response noise and interaction
%parameters similar to those observed in the human participants.
%Then it fits the model on these data and compares configured parameters
%and correspondiong probabilities against fitted parameters and corresponding
%probabilities.


%forensic_beads_simulate_parameters: In order to clean up Git repository
%for public consumption I've renamed this code and removed older versions.
%This code allows one to play with prior and extra guilt and split
%parameters and see effects on sequence position

%Gets conditional probabilities biased by subjective prior

addpath(genpath('C:\matlab_files\fiance\forensic_beads_pub_repo\Forensic-beads-paper-1\FMINSEARCHBND'))

%That is, which sequences / stimuli to use? (we don't need the behaviour)
study_num_to_analyse = 1;   %Can be 1 for Study 1 (Atheism study) or 2 for Study 2 (gender study).


% %data_trunc.xlsx, I formated from data_exp_11596-v23_task-mwjx (1).csv,
% %which Naina acquired from the Gorilla box for this part of the study flowchart.
% %1: event index, 2:participant private id, 3:RT, 4: human probability 5: sequence position (which witness is it?)
% %6: suspect gender (1=female), 7:witness gender (1=female),
% %8: guilty claim (1=guilty), 9: context (mostly innocent / mostly guilty)
% data.raw = xlsread('C:\matlab_files\fiance\forensics_beads_MvF_naina\data_trunc.xlsx');

%stimuli.raw:
%1: event index,
%2:participant private id,
%3:RT,
%4: human probability
%5: sequence position (which witness is it?)
%6: suspect gender (1=female),
%7:witness gender (1=female),
%8: guilty claim (1=guilty),
%9: context (mostly innocent / mostly guilty)
%10 preceding context (degree)
%11 preceding context (category)
data.raw = get_sub_data(study_num_to_analyse);

%borrowing code from forensic_beads_simulate_parameters.m
%We need col 11, which is preceding context
[data.raw(:,[10 11])] = ...
    get_contexts(data.raw);

%You don's care anymore about the human participant data. You only need
%col 2 "participant" (which here is really just a stimulus set), col 6 the key
%manipulation suspect type, col 5 sequence position. col 8, the claims. You
%can compute the preceding context variable ad hoc shortly. And let's avoid the
%confusions with the column indices this time and table it.
stimuli = data.raw(:,2);
suspect = data.raw(:,6);
seq_pos = data.raw(:,5);
claims = data.raw(:,8);
context = data.raw(:,11);

data = table(stimuli,suspect,seq_pos,claims, context);

%generate the preceding context column

% %Make lists of configured parameters
% config_prior = linspace(0,1,4);
% config_interaction = linspace(0,50,4);
% config_response_bias = linspace(0,1,4);
% config_response_noise = linspace(0,1,4);

% %smaller list, for debugging
config_prior = [.4 .6];
config_interaction = [.5 3];
config_response_bias = [0 1];
config_response_noise = [.75 1];

lower_bounds = [0 0 0 0 ];   %fitting will not try parameters below these values
upper_bounds = [1 50 Inf 1];

%initialise stuff
it = 1;
total_params = numel(config_prior)*numel(config_interaction)*numel(config_response_bias)*numel(config_response_noise);
num_params = numel(lower_bounds);
num_priors = numel(config_prior);
num_interacts = numel(config_interaction);
num_biases = numel(config_response_bias);
num_noises = numel(config_response_noise);
stim = unique(stimuli); %actually, participant numbers, but we use their stimuli
num_stim = numel(stim);
fitted_params = NaN(num_priors, num_interacts, num_biases, num_noises,num_params);
%To initialise in this way, I need to assume that every
%stimulus/participant has the same amount of data (e.g., one 11-rating sequence in
%study 1 and 8 11-rating sequences in study 2. So I'm going to grab the
%number of ratings for just the first sub and initialise based on that. Woe
%to the user who wants to mix and match paradigms.
num_trials_per_sub = sum(data.stimuli==stim(1));
configured_probabilities = NaN(num_priors, num_interacts, num_biases, num_noises,num_stim,num_trials_per_sub);
fitted_probabilities = NaN(num_priors, num_interacts, num_biases, num_noises,num_stim,num_trials_per_sub);

%Run model fitting for each parameter level for each stimulus
for prior = 1:num_priors;
    for interaction = 1:num_interacts
        for bias = 1:num_biases;
            for noise = 1:num_noises;
                
                disp(sprintf('test %d of %d',it,total_params));
                it = it+1;
                
                params(1) = config_prior(prior);
                params(2) = config_interaction(interaction);
                params(3) = config_response_bias(bias);
                params(4) = config_response_noise(noise);
                
                %Study 1 uses a new stimulus for every participant so we
                %can loop through participants to get 599 stimuli. If I do
                %the same for study 2 it should produce the right answer
                %but waste a lot of time (and fitting one prior for both suspects). 
                %Study 2 has the same stim for
                %every participant. Use the other parameter recovery code for
                %study 2, which also fits separate priors in each subject.
                for stimulus = 1:num_stim;
                    
                    this_stim_data = data(data.stimuli==stim(stimulus),:);

                    %save for analysis and comparison to fitted_params later
                    configured_params(prior,interaction,bias,noise,stimulus,:) = params;
                    
                    %Generate behavioural data from this configuration of parameters
                    configured_probabilities(prior, interaction, bias, noise,stimulus,:) = ...
                        get_model_behaviour(params, this_stim_data);
                    
                    %setup dataset for fitting
                    clear data_to_fit;
                    data_to_fit.data = this_stim_data;
                    data_to_fit.configured_probabilities = ...
                        squeeze(configured_probabilities(prior, interaction, bias, noise,stimulus,:));
                    
                    %now fit new parameters to the behaviour generated from these "stimuli"
                    [fitted_params(prior,interaction,bias,noise,:), ...
                        ll_temp, flag search] = ...
                        fminsearchbnd( ...
                        @(params) get_model_ll(params, data_to_fit), ...
                        params, ...
                        lower_bounds, ...  %lower parameter bounds
                        upper_bounds ... %upper parameter bounds
                        );
                    
                    %Generate behavioural data from these fitted parameters
                    fitted_probabilities(prior, interaction, bias, noise,stimulus,:) = ...
                        get_model_behaviour( ...
                        fitted_params(prior,interaction,bias,noise,:), ...
                        data_to_fit.data);
                    
                end;    %stimuli / participants
                
            end;    %priors
        end;    %interactions
    end;    %biases
end;    %noise

fprintf('');

[R_params p_params] = correlate_output(configured_params,fitted_params);

% Create a heatmap of the param correlation matrix
f1 = figure('Color',[1 1 1]);

%parameters
% subplot(1,2,1);
hm_param = heatmap(R_params, 'Colormap', cool, 'ColorbarVisible', 'on',  'XLabel', 'Fitted parameters' , 'YLabel', 'Configured parameters');
hm_param.CellLabelFormat = '%2.2f';
ticklabels = {'Prior', 'Interaction', 'Bias', 'Noise'};
hm_param.XDisplayLabels = ticklabels;
hm_param.YDisplayLabels = ticklabels;
caxis([-1, 1]);

%
% %probabilities
% subplot(1,2,2);
% hm_param = heatmap(R_probs, 'Colormap', cool, 'ColorbarVisible', 'on',  'XLabel', 'Fitted parameters' , 'YLabel', 'Configured parameters');
% hm_param.CellLabelFormat = '%2.2f';
% ticklabels = {'Prior', 'Interaction', 'Bias', 'Noise'};
% hm_param.XDisplayLabels = ticklabels;
% hm_param.YDisplayLabels = ticklabels;
%  caxis([-1, 1]);
%

[R_probs p_probs] = corr( ...
    reshape(configured_probabilities,[],1), ...
    reshape(fitted_probabilities,[],1) ...
    );

disp(sprintf('Correlation between fitted and configured probability ratings: r = %0.4f p = %0.4f',R_probs,p_probs));

save('fb_pr_ws.mat')


disp('audi5000');
%%%%%%%%%%%%%%%%%%end, forensic_beads_parameter_recovery%%%%%%%%%%%%%%%%%%%%%%










%%%%%%%%%%%%%%%%%%start, correlate_output%%%%%%%%%%%%%%%%%%%%%%
function [R P] = correlate_output(A,B);

% Define the dimensions of your arrays
dim1 = size(A, 1);
dim2 = size(A, 2);
dim3 = size(A, 3);
dim4 = size(A, 4);
dim5 = size(A, 5);

% Initialize the correlation matrix R
R = zeros(dim5, dim5);

% Reshape A and B into 2D matrices
A_reshaped = reshape(A, [], dim5);
B_reshaped = reshape(B, [], dim5);

% Calculate the correlation for each pair of elements in A and B
for i = 1:dim5
    for j = 1:dim5
        [R(i, j) P(i, j)] = corr(A_reshaped(:, i), B_reshaped(:, j));
    end
end
%%%%%%%%%%%%%%%%%%end, correlate_output%%%%%%%%%%%%%%%%%%%%%%








%%%%%%%%%%%%%%%%%%start, get_model_ll%%%%%%%%%%%%%%%%%%%%%%
function ll = get_model_ll(params, input_data);

this_ps_suspect_data = input_data.data;

fitted_probabilities = get_model_behaviour(params,this_ps_suspect_data);

ll = 0;
for trial = 1:size(this_ps_suspect_data,1);
    
    %model response noise and bias when generating predictions
    y_hat = fitted_probabilities(trial,end)/100;   %end because model rating must be the last col
    
    %Get "labels" for this trial
    y = input_data.configured_probabilities(trial,1)/100;  %human / participant probability is col 4.
    
    ll_this_trial = y*log(y_hat) + (1-y)*log(1-y_hat);
    if ~isreal(ll_this_trial);
        fprintf('');
    end;
    ll = ll - ll_this_trial;
    
end;    %loop through trials
fprintf('');
%%%%%%%%%%%%%%%%%%end, get_model_ll%%%%%%%%%%%%%%%%%%%%%%













%%%%%%%%%%%%%%%%%%start, get_context%%%%%%%%%%%%%%%%%%%%%%
function contexts = get_contexts(stimuli);


%This next loop takes out each sequence, finds the cumulative sum of guilt
%claims for each sequence position and puts it in 9th col of raw, and finds
%cumulative proportion of guilt claims and puts it in 10th column of raw

raw = stimuli;

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

%%%%%%%%%%%%%%%%%%end, get_contexts%%%%%%%%%%%%%%%%%%%%%








%%%%%%%%%%%%%%%%%%start, get_model_behaviour%%%%%%%%%%%%%%%%%%%%%%
function probabilities = get_model_behaviour(params, this_ps_suspect_data)

prior = params(1);
interaction = params(2);
response_bias = params(3);
response_noise = params(4);

%on which indices is the display screen 0 (prior rating prompt so first rating)
seq_start_indices = find(this_ps_suspect_data.seq_pos==0);

% %get number of column
% cols_this_ps_suspect_data = size(this_ps_suspect_data,2);

%For each start index, loop through sequence and get model predictions
for seq = 1:numel(seq_start_indices);
    
    %initialise model probabilities
    
    %     this_ps_suspect_data(seq_start_indices,2) = NaN;
    
    %Loop through this sequence
    for claim = 1:11;
        
        %what's the current index into this_ps_suspect_data?
        index = seq_start_indices(seq)+claim-1;
        %
        %         %save sequence number so can loop more easily later
        %         this_ps_suspect_data(index,cols_this_ps_suspect_data+1) = seq;
        
        %         if claim == 1
        %
        %             this_ps_suspect_data(seq_start_indices(seq),cols_this_ps_suspect_data+2) = prior*100;
        %
        %         else
        
        %get model prediction for every seq position
        q=.6;
        
        %get number of guilts (i.e., the number of 1s)
        ng = sum( this_ps_suspect_data.claims(seq_start_indices(seq)+1:index) ) ;
        
        %get number of draws so far
        nd = claim-1;
        
        %             %interaction term
        this_claim = sum(this_ps_suspect_data.claims(index));  %guilt or innocent claim now?
        if isnan(this_claim);   %first rating in sequence, before any witnesses
            interactTerm = 0;   %no contribution of interaction yet
        else;
            interactTerm = this_claim*ng;   %claim * preceding context interaction
        end;
        
        %assign model probability
        %             noiseless_p = ( (1/(1 + ((1-prior)/prior)*(q/(1-q))^(nd-2*ng)))*100 ) + (interaction*ng);
        noiseless_p = ( (1/(1 + ((1-prior)/prior)*(q/(1-q))^(nd-2*ng)))*100 ) + (interaction*interactTerm); %free parameter is weight on context*claim interaction
        
        
        %add noise and response bias
        probabilities(index,1) = ...
            response_bias + response_noise*noiseless_p;
        
        %         end;    %before first claim (sequence position 0) or a later one?
        
    end;    %loop through this sequence (claim)
    
end;    %loop through sequences
fprintf('');
%%%%%%%%%%%%%%%%%%end, simulate%%%%%%%%%%%%%%%%%%%%%%











%%%%%%%%%start, get_sub_data%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function data = get_sub_data(study_num_to_analyse);

if study_num_to_analyse == 1;
    
    fnames = {...
        'C:\matlab_files\fiance\forensic_beads_pub_repo\Forensic-beads-paper-1\Study1\01_christi_mostly_innoce.xlsx'...
        'C:\matlab_files\fiance\forensic_beads_pub_repo\Forensic-beads-paper-1\Study1\02_atheist_mostly_innoce.xlsx'...
        'C:\matlab_files\fiance\forensic_beads_pub_repo\Forensic-beads-paper-1\Study1\03_christi_mostly_guilty.xlsx'...
        'C:\matlab_files\fiance\forensic_beads_pub_repo\Forensic-beads-paper-1\Study1\04_atheist_mostly_guilty.xlsx'...
        };
    sus_rel_codes = [1 0 1 0];
    context_codes = [0 0 1 1];
    data = [];
    for file = 1:4;
        
        clear temp;
        temp = xlsread(fnames{file});
        
        data = [...
            data;
            temp(:,6) ...                                   %1: event index
            temp(:,1) ...                                   %2: participant private id
            temp(:,2) ...                                   %3: RT
            temp(:,3) ...                                   %4: probability estimate
            repmat([0:10]',size(temp,1)/11 ,1) ...          %5: sequence position (including prior 0-10)
            sus_rel_codes(file)*ones(size(temp,1),1) ...    %6: suspect religion (0=Atheist)
            temp(:,5) ...                                   %7: witness gender (1=female)
            temp(:,4) ...                                   %8: witness claim (1=guilt)
            context_codes(file)*ones(size(temp,1),1) ...    %9: context (1=mosty guilty)
            ];
        
    end;    %Loop through datafiles (file)
    
    
    
elseif study_num_to_analyse == 2;
    
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
    data = xlsread('C:\matlab_files\fiance\forensic_beads_pub_repo\Forensic-beads-paper-1\Study2\data_trunc.xlsx');
    
end;    %Which study's dataset?


%In raw, sequence positions 0 have NaNs in place of condition labels
%for contexts (col 9) and sometimes suspects (col 6). Put
%them back in or you'll have troubles later
nan_indices = find(data(:,5) == 0);    %find NaNs
data(nan_indices,9) = data(nan_indices+1,9);  %assign the missing values at pos 0 with the values at pos 1
data(nan_indices,6) = data(nan_indices+1,6);  %assign the missing values at pos 0 with the values at pos 1

%Get context operationalised according to preceding claims only
%binarised as guilty innocent (10) or numeric as number of guilts (11)
[data(:,[10 11])] = ...
    get_contexts(data);

%%%%%%%%%end, get_sub_data%%%%%%%%%%%%%%%%%%%%%%%%%%%%






