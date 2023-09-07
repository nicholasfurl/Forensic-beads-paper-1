%%%%%%%%%%%%%%%%%%start, forensic_beads_prior_sim%%%%%%%%%%%%%%%%%%%%%%
function forensic_beads_prior_sim_fit_multiparam_study2;

%I've given up fighting with the code and multiplying its complexity to
%accomodate one model that fits per suspects per participant (study 2) and
%another that works for one suspect per partici[pant (study 1). Let's just
%use forensic_beads_prior_sim_fit_multiparam.m for Study 1 and for Study 2
%if you want to look at separate parameter fits for all parameters in the
%different suspects. Use forensic_beads_prior_sim_fit_multiparam_study.m if
%you want study 2 where suspect parameters vary over participants but
%there's one fit per participant for response bias, response noise and
%disconfirmatory.

%forensic_beads_prior_sim_fit_multiparam_2studies.m expands the version of
%forensic_beads_prior_sim_fit_multiparam.m that existed on 2/Sept/2023 to
%integrate in Study 1 as well as Study 2 (which was the only one
%implemented in the previous programs).

%forensic_beads_prior_sim_fit_multiparam.m expands the version of
%forensic_beads_prior_sim_fit.m that existed on 2/Sept/2023 to introduce a
%number of guilt claims increment parameter.

%forensic_beads_prior_sim_fit.m converts simulation code to model fitting
%code. Beginning with just fitting of the prior in the two suspect conditions.

%forensic_beads_prior_sim_fit.m - simulates conditional probabilities biased by subjective prior

addpath(genpath('C:\matlab_files\fiance\forensic_beads_pub_repo\Forensic-beads-paper-1\FMINSEARCHBND'))

% study_num_to_analyse = 2;   %Can be 1 for Study 1 (Atheism study) or 2 for Study 2 (gender study).

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
stimuli.raw = get_sub_data(2);

%Who are the participants?
participant_list = unique(stimuli.raw(:,2),'stable');
num_participants = numel(participant_list);

%initial value of free params
params(1) = .5; %prior, initialised to optimal value (ground truth of paradigm)
params(2) = .5; %prior, initialised to optimal value (ground truth of paradigm)
params(3) = 0;  %guilt claim increment, intitialised to optimal value
params(4) = 0;  %bias term, intialised to optimal value
params(5) = 1;  %noise term, initialised to optimal value

lower_bounds = [0 0 0 0 0];   %fitting will not try parameters below these values
upper_bounds = [1 1 Inf Inf 1];

%indices into params that designate which are free. Handy way to play
%around with models by changing parameterisation. "Initial" values in
%params become hard coded if not indexed here. Use empty array with ground
%truth settings in params to get ideal observer.
free_params_idx = [3];

num_params = numel(params);

% num_suspects = 2;   %suspect type (atheist/Christian or male/female)

%initialise data matrix to hold modelling results, with each participant in a row
model_fitting_results = [];
model_behaviour_results = [];

%initialise matrix to hold parameters. Ensure it is filled with NaNs first
%so it's flexible enough to handle within (Study 2) or between (Study 1) participant
%manipulations of suspect
% params_est = nan(num_participants,num_suspects,num_params);
% ll = nan(num_participants,num_suspects);

for participant = 1:num_participants;
    
    %get data for this participant
    this_ps_data = stimuli.raw(find(stimuli.raw(:,2) == participant_list(participant)),:);
    
    
    disp(sprintf('fitting participant %d', participant))
    
    %split free and fixed params
    free_params = params(free_params_idx);
    
    free_lower_bounds = lower_bounds(free_params_idx);
    free_upper_bounds = upper_bounds(free_params_idx);
 
    %pass data, initialised param and function handle to fminsearch
    [params_temp, ll_temp, flag search] = ...
        fminsearchbnd( ...
        @(free_params) get_model_ll(free_params, params, free_params_idx, this_ps_data), ...
        free_params, ...
        free_lower_bounds, ...  %lower parameter bounds
        free_upper_bounds ... %upper parameter bounds
        );
 
    %Now that this participant has been fit, get model performance
    %model_fitting_results:
    %col1: participant id, col2: suspect code, col3: ll, cols 4 to end: params
    model_fitting_results = ...
        [model_fitting_results; ...
        [participant_list(participant) ll_temp params_temp]
        ];
    
    %get performance for this model
    %This needs to be divided by suspect, once for each of the
    %suspect-specific parameters,both of which are estimated on this loop
       
    %Accumulate results for plotting. Need to re-loop suspects if Study 2
    this_ps_suspect_codes = unique(rmmissing(this_ps_data(:,6)),'stable');     %get suspect codes present in this participant
    this_ps_num_suspects = numel(this_ps_suspect_codes);              %How many codes for this participant?
    
    for suspect = 1:this_ps_num_suspects;
        
        clear this_suspect_data
        this_ps_suspect_data = this_ps_data(this_ps_data(:,6) == this_ps_suspect_codes(suspect),:); %should be same suspect-sepcific data as used to estiate parameter inside of fitting function above
        params_performance = params;    %So this just initialises the params vector to be used with all the iniatalised (not fitted) values, as some of the paraneters will be fitted and some fixed
        params_performance(free_params_idx) = params_temp;  %replace default parameter list with any updated free parameter values
        params_performance = params_performance([suspect 3:end]); %Narrow down to just this suspect
        temp = get_model_behaviour(params_performance,this_ps_suspect_data);
        model_behaviour_results = [ model_behaviour_results; temp ];
        
    end;    %loop through suspects
    
    %     save('test_fit_multiparam.m');
    
    end;    %loop through participants
    
    %Before plotting, you need to compute human and model adjustments
    %returns model_behaviour_results (each row is a sequence position):
%1: event index,
%2:participant private id,
%3:RT,
%4: human probability
%5: sequence position (which witness is it?)
%6: suspect (1=female, 1 = Christian),
%7:witness gender (1=female),
%8: guilty claim (1=guilty),
%9: context (mostly innocent / mostly guilty)
%10 preceding context (degree)
%11 preceding context (category)
%12 sequence number (per suspect)
%13 model's probability
%14 human's adjustments
%15 model's adjustments
model_behaviour_results = get_adjustments(model_behaviour_results);




%%%%%%%%%%%%%%PARAMETER PLOT BEGIN%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Alright, at this stage, I have estimated parameters in the columns of
%model_fitting_results. Can make plot expandable in case number of
%parameters changes. Will need to make multiple plots / bars for the
%different suspects.

%groupby and mean aggregate by participant
[param_means param_ci] = grpstats(model_fitting_results(:,3:end), {ones(size(model_fitting_results,1),1)},{'mean', 'meanci'});

%plot
h1 = figure('Color',[1 1 1]);

input_struct.fig = h1;
input_struct.sp = [1,1,1];
input_struct.input_means = param_means';
if numel(param_means) == 1;
    input_struct.input_ci = param_ci(1) - param_means;
else
    input_struct.input_ci = squeeze(param_ci(:,:,1))'-param_means';
end;
input_struct.xlabel = 'Parameter';
input_struct.ylabel = 'Parameter value';

b = make_grouped_bar_with_errors(input_struct);
%%%%%%%%%%%%%%PARAMETER PLOT ENDS%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



%%%%%%%%%%%%%%SEQUENCE POSITION PLOT BEGIN%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%I have estimated behaviour at different sequence positions in rows of
%model_behavioural results and codes for claims, contexts and suspects for each row
%Can overlay line plots for different suspects
%in different contexts in one plot. In another, I can plot claims *
%contexts.
%
% %get participant averages
cols_to_use = [6 9 5];  %suspect, context, seq pos
groupvars = { ...
    model_behaviour_results(:,cols_to_use(1)) ...
    model_behaviour_results(:,cols_to_use(2)) ...
    model_behaviour_results(:,cols_to_use(3)) ...
    };   %suspect, context, sequence position
[means meancis] = grpstats(model_behaviour_results,groupvars,{'mean','meanci'});

%now get means and ci's over participants

%Loop through and plot probabilities at sequence positions for different suspects
figure('Color',[1 1 1]);
suspects = unique(means(:,6),'stable');
contexts = unique(means(:,9),'stable');

for suspect = 1:numel(suspects);
    for context = 1:numel(contexts);
        
        this_data = means(means(:,6) == suspects(suspect) & means(:,9) == contexts(context),13);
        %this_data_ci = meancis(means(:,6) == suspects(suspect),10,1);
        plot(this_data); hold on;
        
    end;    %contexts loop
end;    %suspects loop
ylim([1 100]);
box off;
%%%%%%%%%%%%%%SEQUENCE POSITION PLOT END%%%%%%%%%%%%%%%%%%%%%%%%%%%%%





%%%%%%%%%%%%%%DISCONFIRMATORY ADJUSTMENT PLOT BEGIN%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% %collapse over sequence positions
claim_idx = 8;
%context_idx = 9;  %context (mostly innocent / guilty)
context_idx = 11;  %preceding context (category)
suspect_idx = 6;
participant_idx = 2;
adjustment_idx = 15;

% %get participant averages
% cols_to_use = [context_idx claim_idx participant_idx adjustment_idx]; %context (ad hoc), claim, participant, adjustment
% groupvars = { model_behaviour_results(:,cols_to_use(1)) model_behaviour_results(:,cols_to_use(2)) model_behaviour_results(:,cols_to_use(3)) };   %context, claim, participant
% temp = grpstats(model_behaviour_results, groupvars,'mean');
% 
% %get participant averages
% % groupvars = {temp(:,1) temp(:,2)};   %suspect, context
% groupvars = {temp(:,11) temp(:,8)};   %context, claim
% [means meancis] = grpstats(temp,groupvars,{'mean','meanci'});
% 
% figure('Color',[1 1 1]);
% contexts = unique(means(:,11),'stable');
% claims = unique(means(:,8),'stable');
% 
% for context = 1:numel(contexts);
%     
%     this_data = means(means(:,11) == contexts(context),15);
%     %this_data_ci = meancis(means(:,6) == suspects(suspect),10,1);
%     plot(this_data); hold on;
%     
% end;    %contexts loop
% 
% ylim([-15 15]);
% xlim([0.5 2.5]);
% set(gca,'XTick',[1 2]);
% xticklabels({'innocent' 'guilty'});
% legend({'innocent context' 'guilty context'});
% xlabel('Claim');
% box off;





group_vars = { ...
    model_behaviour_results(:,context_idx), ... %claim
    model_behaviour_results(:,claim_idx), ...  %context (mostly innocent / guilty)
    model_behaviour_results(:,suspect_idx), ...  %suspect
    model_behaviour_results(:,participant_idx), ...  %participant
    };
 
%collapse over sequence lengths
mbr_collapse_seqpos = grpstats(model_behaviour_results,group_vars,'mean');
%Now get means with confidence intervals computed over participants
group_vars = { ...
    mbr_collapse_seqpos(:,context_idx), ... %claim
    mbr_collapse_seqpos(:,claim_idx), ...  %preceding context (category)
    mbr_collapse_seqpos(:,suspect_idx), ...  %suspect
    };
[mbr_disconfirm_means, mbr_disconfirm_cis] = ...
    grpstats(mbr_collapse_seqpos,group_vars,{'mean' 'meanci'});

h2 = figure('Color',[1 1 1]);

suspects = unique(mbr_disconfirm_means(:,suspect_idx),'stable');
suspects_num = numel(suspects);

%make suplots for each suspect
for suspect = 1:suspects_num;
    
    this_suspect_data = mbr_disconfirm_means(mbr_disconfirm_means(:,suspect_idx)==suspects(suspect),adjustment_idx);  %means
    this_suspect_ci = mbr_disconfirm_cis(mbr_disconfirm_means(:,suspect_idx)==suspects(suspect),adjustment_idx,1) - this_suspect_data;  %means
    
    this_suspect_data = reshape(this_suspect_data',[2,2]);
    this_suspect_ci = reshape(this_suspect_ci',[2,2]);
    
    input_struct.fig = h2;
    input_struct.sp = [1,2,suspect];
    input_struct.input_means = this_suspect_data;
    input_struct.input_ci = this_suspect_ci;
    input_struct.xlabel = 'Claim Type';
    input_struct.ylabel = 'Adjustment';
    
    b = make_grouped_bar_with_errors(input_struct);

end;
% %%%%%%%%%%%%%%DISCONFIRMATORY ADJUSTMENT PLOT END%%%%%%%%%%%%%%%%%%%%%%%%%%%%



disp('audi5000');
%%%%%%%%%%%%%%%%%%end, forensic_beads_prior_sim_fit_multiparam_2studies%%%%%%%%%%%%%%%%%%%%%%




%%%%%%%%%%%%%%%%%%start, get_model_ll%%%%%%%%%%%%%%%%%%%%%%
function ll = get_model_ll(free_params,params, free_params_idx, this_ps_data);

%replace default parameter list with any updated free parameter values
params(free_params_idx) = free_params;

%This participant might have both suspects (Study 2) or just one (Study 1)
this_ps_suspect_codes = unique(rmmissing(this_ps_data(:,6)),'stable');     %get suspect codes present in this participant
this_ps_num_suspects = numel(this_ps_suspect_codes);              %How many codes for this participant?

ll = 0;
%Now loop through the detected conditions
for suspect = 1:this_ps_num_suspects;
    
    this_ps_suspect_data = this_ps_data(this_ps_data(:,6) == this_ps_suspect_codes(suspect),:);
    
    %alter parameter list to be specific for this suspect
    this_params = [params(suspect) params(3:end)];   %The current free parameter value of the suspect prior tested in this suspect loop iteration and the other current values of parameters
    %this ps_suspect_data:
    %now the same as raw, but adds cols col 12: seq num, col 13: model rating
    this_ps_suspect_data = get_model_behaviour(this_params,this_ps_suspect_data);
    
    for trial = 1:size(this_ps_suspect_data,1);
        
        %model response noise and bias when generating predictions
        y_hat = this_ps_suspect_data(trial,end)/100;   %end because model rating must be the last col
        
        %Get "labels" for this trial
        y = this_ps_suspect_data(trial,4)/100;  %human / participant probability is col 4.
        
        ll_this_trial = y*log(y_hat) + (1-y)*log(1-y_hat);
        if ~isreal(ll_this_trial);
            fprintf('');
        end;
        ll = ll - ll_this_trial;
        
    end;    %loop through trials, this suspect
    
    
end;


fprintf('');
%%%%%%%%%%%%%%%%%%end, get_model_ll%%%%%%%%%%%%%%%%%%%%%%











%%%%%%%%%%%%%%%%%%start, get adjustments%%%%%%%%%%%%%%%%%%%%%%
function model_behaviour_results = get_adjustments(model_behaviour_results);

%on which indices is the display screen 0 (prior rating prompt so first rating)
seq_start_indices = find(model_behaviour_results(:,5)==0);

%In case I need to change cols of input later
cols_this_ps_suspect_data = size(model_behaviour_results,2);

%For each start index, loop through sequence
for seq = 1:numel(seq_start_indices);
    
    %Loop through this sequence
    for claim = 1:11;
        
        
        %what's the current index into this_ps_suspect_data?
        index = seq_start_indices(seq)+claim-1;
        
        if claim == 1
            
            model_behaviour_results( ...
                seq_start_indices(seq), ...
                [cols_this_ps_suspect_data+1 cols_this_ps_suspect_data+2]) ...
                = NaN;
            
        else
            
            %human adjustments
            model_behaviour_results(index,cols_this_ps_suspect_data+1) = ...
                model_behaviour_results(index,4) -  model_behaviour_results(index-1,4);
            
            %model adjustments
            model_behaviour_results(index,cols_this_ps_suspect_data+2) = ...
                model_behaviour_results(index,13) -  model_behaviour_results(index-1,13);
            
        end;    %sequence position 0 (prior) or another? claim == 1?
        
    end;    %claims
    
end;    %sequences

fprintf('');

%%%%%%%%%%%%%%%%%%end, get adjustments%%%%%%%%%%%%%%%%%%%%%%





%%%%%%%%%%%%%%%%%end, make_grouped_bar_with_errors%%%%%%%%%%%%%%%%%%%%%%
function sp_h = make_grouped_bar_with_errors(input_struct);

figure(input_struct.fig);

sp_h = subplot(input_struct.sp(1), input_struct.sp(2), input_struct.sp(3));

b = bar(input_struct.input_means,'grouped');
hold on;

% Calculate the number of groups and number of bars in each group
[ngroups,nbars] = size(input_struct.input_means);

% Get the x coordinate of the bars
x = nan(nbars, ngroups);
for i = 1:nbars
    x(i,:) = b(i).XEndPoints;
end

% Plot the errorbars
errorbar(x',input_struct.input_means,input_struct.input_ci,'k','linestyle','none');

box off;
xlabel(input_struct.xlabel);
ylabel(input_struct.ylabel);
%%%%%%%%%%%%%%%%%end, make_grouped_bar_with_errors%%%%%%%%%%%%%%%%%%%%%%










%%%%%%%%%%%%%%%%%%start, get_model_behaviour%%%%%%%%%%%%%%%%%%%%%%
function this_ps_suspect_data = get_model_behaviour(params, this_ps_suspect_data)

prior = params(1);
guilt_claim_inc = params(2);

%on which indices is the display screen 0 (prior rating prompt so first rating)
seq_start_indices = find(this_ps_suspect_data(:,5)==0);

%sometimes I need to change the columns of this_ps_suspect_data outside of
%this function. But this function adds a column to the right end of
%this_ps_suspect_data. So let's make get_model_behaviour more adaptable.
%It'll add on to whatever number of columns are here
cols_this_ps_suspect_data = size(this_ps_suspect_data,2);

%For each start index, loop through sequence and get model predictions
for seq = 1:numel(seq_start_indices);
    
    %initialise model probabilities
    
    %     this_ps_suspect_data(seq_start_indices,2) = NaN;
    
    %Loop through this sequence
    for claim = 1:11;
        
        %what's the current index into this_ps_suspect_data?
        index = seq_start_indices(seq)+claim-1;
        
        %save sequence number so can loop more easily later
        this_ps_suspect_data(index,cols_this_ps_suspect_data+1) = seq;
        
        %         if claim == 1
        %
        %             this_ps_suspect_data(seq_start_indices(seq),cols_this_ps_suspect_data+2) = prior*100;
        %
        %         else
        
        %get model prediction for every seq position
        q=.6;
        
        %get number of guilts (i.e., the number of 1s)
        ng = sum( this_ps_suspect_data(seq_start_indices(seq)+1:index,8) )   + guilt_claim_inc;
        
        %get number of draws so far
        nd = claim-1;
        
        %assign model probability
        noiseless_p = (1/(1 + ((1-prior)/prior)*(q/(1-q))^(nd-2*ng)))*100;
        
        %add noise and response bias
        this_ps_suspect_data(index,cols_this_ps_suspect_data+2) = ...
            params(3) + params(4)*noiseless_p;
        
        %         end;    %before first claim (sequence position 0) or a later one?
        
    end;    %loop through this sequence (claim)
    
end;    %loop through sequences
fprintf('');
%%%%%%%%%%%%%%%%%%end, simulate%%%%%%%%%%%%%%%%%%%%%%






%%%%%%%%%%%%%%%%%%start, get_context%%%%%%%%%%%%%%%%%%%%%%
function contexts = get_contexts(raw);


%This next loop takes out each sequence, finds the cumulative sum of guilt
%claims for each sequence position and puts it in 9th col of raw, and finds
%cumulative proportion of guilt claims and puts it in 10th column of raw

% raw = stimuli.raw;

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
    
    contexts(seq_starts(sequence,1):seq_starts(sequence,1)+10,2) = temp_MG; %CONTEXT CATEGORY
    
end;    %loop seq starts

%%%%%%%%%%%%%%%%%%end, get_context%%%%%%%%%%%%%%%%%%%%%%




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



