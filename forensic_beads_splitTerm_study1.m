%%%%%%%%%%%%%%%%%%start, forensic_beads_prior_sim%%%%%%%%%%%%%%%%%%%%%%
function forensic_beads_prior_sim_fit_multiparam_2studies_interacterm;

%forensic_beads_prior_sim_fit_multiparam_2studies_interacterm.m adds a new
%parameter to make another stab at modelling the underadjustment to
%confirmatory evidence effect.

%forensic_beads_prior_sim_fit_multiparam_2studies.m
%is the one you can use to either run Study 1 or run Study 2 such that
%you have a separate set of all parameters fitted for each suspect. If you
%want Study 2 where prior varies over suspects and the other parameters are
%one fitted per participant, use
%forensic_beads_prior_sim_fit_multiparam_study2.m

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
addpath(genpath('C:\matlab_files\fiance\forensic_beads_pub_repo\Forensic-beads-paper-1\klabhub-bayesFactor-3d1e8a5'))
addpath(genpath('C:\matlab_files\fiance\forensic_beads_pub_repo\Forensic-beads-paper-1\shaded_plots'));

study_num_to_analyse = 1;   %Can be 1 for Study 1 (Atheism study) or 2 for Study 2 (gender study).

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
stimuli.raw = get_sub_data(study_num_to_analyse);

% %In raw, sequence positions 0 have NaNs in place of condition labels
% %for contexts (col 9) and sometimes suspects (col 6). Put
% %them back in or you'll have troubles later
% nan_indices = find( stimuli.raw(:,5) == 0);    %find NaNs
% stimuli.raw(nan_indices,9) = stimuli.raw(nan_indices+1,9);  %assign the missing values at pos 0 with the values at pos 1
% stimuli.raw(nan_indices,6) = stimuli.raw(nan_indices+1,6);  %assign the missing values at pos 0 with the values at pos 1

%
% [stimuli.raw(:,[10 11])] = ...
%     get_contexts(stimuli);

%Who are the participants?
participant_list = unique(stimuli.raw(:,2));
num_participants = numel(participant_list);

%initial value of free params
params(1) = .5; %prior, initialised to optimal value (ground truth of paradigm)
params(2) = 0;  %guilt claim increment, intitialised to optimal value
params(3) = 0;  %bias term, intialised to optimal value
params(4) = 1;  %noise term, initialised to optimal value
params(5) = 0;  %interaction term (weight on amount of confirmation). Starting value of 0 is the optimal value.
params(6) = .6;    %split term, initialised to optimal value

lower_bounds = [.5 0 0 1 0 .6];   %fitting will not try parameters below these values
upper_bounds = [.5 0 0 1  0 6];

%indices into params that designate which are free. Handy way to play
%around with models by changing parameterisation. "Initial" values in
%params become hard coded if not indexed here.
free_params_idx = [1 3 4 5 6];

num_params = numel(params);

num_suspects = 2;   %suspect type (atheist/Christian or male/female)

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
    
    %This participant might have both suspects (Study 2) or just one (Study 1)
    this_ps_suspect_codes = unique(rmmissing(this_ps_data(:,6)));     %get suspect codes present in this participant
    this_ps_num_suspects = numel(this_ps_suspect_codes);              %How many codes for this participant?
    
    %Now loop through the detected conditions
    for suspect = 1:this_ps_num_suspects;
        
        disp(sprintf('fitting participant %d suspect %d', participant, suspect))
        
        %Get the probability rating data to fit for this suspect in this participant
        this_ps_suspect_data = this_ps_data(this_ps_data(:,6) == this_ps_suspect_codes(suspect),:);
        
        %split free and fixed params
        free_params = params(free_params_idx);
        %         fixed_params = params(setdiff(1:end,free_params_idx));
        
        free_lower_bounds = lower_bounds(free_params_idx);
        free_upper_bounds = upper_bounds(free_params_idx);
        
        %pass data, initialised param and function handle to fminsearch
        %         [params_est(participant, this_ps_suspect_codes(suspect)+1,:) ll(participant, this_ps_suspect_codes(suspect)+1) flag search] = ...
        [new_params, new_ll, flag search] = ...
            fminsearchbnd( ...
            @(free_params) get_model_ll(free_params, params, free_params_idx, this_ps_suspect_data), ...
            free_params, ...
            free_lower_bounds, ...  %lower parameter bounds
            free_upper_bounds ... %upper parameter bounds
            );
        
            %In the simplify version, we take the parameter value list, replace the
    %intial values with the updated ones and then assign that whole thing,
    %free and fixed together right now to model fitting results.
    params_temp = params;  %So this just initialises the params vector to be used with all the initialised (not fitted) values, as some of the parameters will be fitted and some fixed
    params_temp(free_params_idx) = new_params;  %replace default parameter list with any updated free parameter values
    
    %Now that this participant has been fit, get model performance
    %model_fitting_results:
    %col1: participant id, 
    %col2: suspect code, 
    %col3: ll, 
    %cols 4 to end: params
    model_fitting_results = ...
        [model_fitting_results; ...
        [participant_list(participant) this_ps_suspect_codes(suspect) new_ll params_temp]
        ];
        
%         %Accumulate results for plotting
%         
%         %Now that this participant has been fit, get model performance
%         %model_fitting_results:
%         %col1: participant id, col2: suspect code, col3: ll, cols 4 to end: params
%         model_fitting_results = ...
%             [model_fitting_results; ...
%             [participant_list(participant) this_ps_suspect_codes(suspect) ll_temp params_temp]
%             ];
        
        %get performance for this model
%         params_performance = params;
%         params_performance(free_params_idx) = params_temp;
        temp = get_model_behaviour(params_temp,this_ps_suspect_data);
        model_behaviour_results = [ model_behaviour_results; temp ];
        
    end;    %loop through suspects
    
    save('test_fit_multiparam.m');
    
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
model_behaviour_results = get_adjustments(model_behaviour_results); %should tack human and model adjustment cols onto the end




%%%%%%%%%%%%%%PARAMETER PLOT BEGIN%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Alright, at this stage, I have estimated parameters in the columns of
%model_fitting_results. Can make plot expandable in case number of
%parameters changes. Will need to make multiple plots / bars for the
%different suspects.
% 
% %Some issues when there's only one col (e.g., only one paranmeter is estimated), as tends to happen with matlab
% parameter_values = model_fitting_results(:,3:end);  %ll through all params
% num_params_test = size(parameter_values,2);

%indices into model_fitting_results that I'll need (for transparency)
mfr_prior_param_index = 4;  %index into model_fitting results dim 2 for prior param
mfr_suspect_index = 2;
mfr_ll_index = 3;

%groupby and mean aggregate by suspect code:
% [param_means param_ci] = grpstats(parameter_values, {model_fitting_results(:,2)},{'mean', 'meanci'});
[param_means param_ci] = grpstats(model_fitting_results, {model_fitting_results(:,mfr_suspect_index)},{'mean', 'meanci'});


%plot
h1 = figure('Color',[1 1 1]);
input_struct.fig = h1;

%parameter value plot
input_struct.sp = [1,2,1];
%make_grouped_bar_with_errors wants params in rows and suspects in cols
% if num_params_test == 1;
%     input_struct.input_means = param_means;
%     input_struct.input_ci = (param_ci(:,2)-param_means);
% else
    input_struct.input_means = ...
        param_means(:,mfr_prior_param_index:end)';
    input_struct.input_ci = ...
        param_ci(:,mfr_prior_param_index:end,2)' - ...
        input_struct.input_means;
% end;
input_struct.xlabel = 'Parameter';
input_struct.ylabel = 'Mean Parameter value';
input_struct.title = 'Parameter values';
b = make_grouped_bar_with_errors(input_struct);

%log-likelihood plot
input_struct.sp = [1,2,2];

    input_struct.input_means = ...
        param_means(:,mfr_ll_index)';
    input_struct.input_ci = ...
        param_ci(:,mfr_ll_index,2)' - ...
        input_struct.input_means;

input_struct.xlabel = '';
input_struct.ylabel = 'negative log-likelihood';
input_struct.title = 'Model fit';
b = make_grouped_bar_with_errors(input_struct);


% %trad paired t-test
[th tpvals tci tstats] =   ...
    ttest2( ...
    model_fitting_results(model_fitting_results(:,mfr_suspect_index)==0,mfr_prior_param_index), ...
    model_fitting_results(model_fitting_results(:,mfr_suspect_index)==1,mfr_prior_param_index) ...
    );

% %paired Bayesian test
[bf10,bfpvals,bfci,bfstats] = ...
    bf.ttest2( ...
    model_fitting_results(model_fitting_results(:,mfr_suspect_index)==0,mfr_prior_param_index), ...
    model_fitting_results(model_fitting_results(:,mfr_suspect_index)==1,mfr_prior_param_index) ...
    );

disp(sprintf('Suspect effect on prior: trad pval: %0.2f bf10: %7.4f',tpvals,bf10));


%Ok now let's make a fancy kernel density plot for the pub
cmap = [0 0 0;
    .25 .25 .25
    ];

alphas = [0 .75];

h1 = figure('Color',[1 1 1]);

for suspect = 1:2;  %suspect parameters loop
    
    %Kernel density
    [kernel(:,suspect) xi(:,suspect)] = ...
        ksdensity( ...
        model_fitting_results(model_fitting_results(:,mfr_suspect_index)==suspect-1,mfr_prior_param_index) ...
        );
    plot_shaded(xi(:,suspect),kernel(:,suspect),'Alpha',alphas(suspect),'Color',cmap(suspect,:),'LineWidth',2);
    hold on;
    
    %     %binned histogram
    %      plot_histogram_shaded(model_fitting_results(:,2+suspect),'bins',25,'alpha',alphas(suspect),'color',cmap(suspect,:));
    
    
end;    %suspect parameters loop

box off;
xlabel('Estimated prior parameter value');
ylabel('Probability (kernel) density');

fprintf('');

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
% cmap = lines(4);
cmap = [0 0 0;
    0 0 0;
    .5 .5 .5;
    .5 .5 .5
    ];
figure('Color',[1 1 1]);

legend_labels = {'suspect 0, innocent sequence' 'suspect 0, guilty sequence' 'suspect 1, innocent sequence'  'suspect 1, guilty sequence'};

suspects = unique(means(:,cols_to_use(1)));
contexts = unique(means(:,cols_to_use(2)));

cmap_it = 1;
for suspect = 1:numel(suspects);
    for context = 1:numel(contexts);
        
        this_data = means(means(:,cols_to_use(1)) == suspects(suspect) & means(:,cols_to_use(2)) == contexts(context),13);
        this_cis = ...
            meancis(meancis(:,cols_to_use(1)) == suspects(suspect) & meancis(:,cols_to_use(2)) == contexts(context),13,2) - this_data;
        %this_data_ci = meancis(means(:,6) == suspects(suspect),10,1);
        
        errorbar( ...
            this_data, ...
            this_cis, ...
            'Color', cmap(cmap_it,:) ...
            );
        hold on;
        plot( ...
            this_data, ...
            'Marker','o', ...
            'Color', cmap(cmap_it,:) ...
            );
        
        hold on;
        
        text(2,30-cmap_it*6,legend_labels{cmap_it},'Color',cmap(cmap_it,:));
        
        cmap_it = cmap_it + 1;
        
    end;    %contexts loop
end;    %suspects loop


ylim([1 100]);
box off;
%legend({'suspect 0, innocent sequence' 'suspect 1, innocent sequence' 'suspect 0, guilty sequence' 'suspect 1, guilty sequence'});
xlabel('Sequence position');
ylabel('Probability');

fprintf('');
%%%%%%%%%%%%%%SEQUENCE POSITION PLOT END%%%%%%%%%%%%%%%%%%%%%%%%%%%%%





%%%%%%%%%%%%%%DISCONFIRMATORY ADJUSTMENT PLOT BEGIN%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Get variables to average for human and model needed to collapse over sequence lengths
group_vars = { ...
    model_behaviour_results(:,8), ... %claim
    model_behaviour_results(:,11), ...  %preceding context (category)
    model_behaviour_results(:,6), ...  %suspect
    model_behaviour_results(:,2), ...  %participant
    };

%Now collapse over sequence lengths
mbr_collapse_seqpos = grpstats(model_behaviour_results,group_vars,'mean');

%Now get means with confidence intervals to get data per participant
group_vars = { ...
    mbr_collapse_seqpos(:,8), ... %claim
    mbr_collapse_seqpos(:,11), ...  %preceding context (category)
    mbr_collapse_seqpos(:,6), ...  %suspect
    };

%Now, get data per participant
[mbr_disconfirm_means, mbr_disconfirm_cis] = ...
    grpstats(mbr_collapse_seqpos,group_vars,{'mean' 'meanci'});

h2 = figure('Color',[1 1 1]);

suspects = unique(mbr_disconfirm_means(:,6));
suspects_num = numel(suspects);

%make suplots for each suspect
subplot_it = 1;
titles = {'Model suspect 0', 'Human suspect 0' 'Model suspect 1', 'Human suspect 1' };
for suspect = 1:suspects_num;
    
    %setup plot
    input_struct.fig = h2;
    input_struct.xlabel = 'Claim Type';
    input_struct.ylabel = 'Adjustment';
    
    %indices for getting plot data (just to be organised)
    col_index_model_adj = 15;
    col_index_human_adj = 14;
    col_index_claims = 8;
    col_index_context = 11;
    col_index_suspects = 6;
    this_suspect_indices = mbr_disconfirm_means(:,col_index_suspects)==suspects(suspect);
    
    %Get model data
    this_suspect_data = ...
        mbr_disconfirm_means( ...
        this_suspect_indices, ...
        [col_index_model_adj col_index_claims col_index_context]...
        );  %means
    
    this_suspect_ci = ...
        mbr_disconfirm_cis( ...
        this_suspect_indices, ...
        [col_index_model_adj col_index_claims col_index_context],2 ...
        );  %CIs
    
    this_suspect_ci = reshape( ...
        this_suspect_ci(:,1) - this_suspect_data(:,1), ...
        [2,2] ...
        );
    
    this_suspect_data = reshape( ...
        this_suspect_data(:,1), ...
        [2,2] ...
        );
    
    input_struct.input_means = this_suspect_data';
    input_struct.input_ci = this_suspect_ci';
    input_struct.title = titles{subplot_it};
    input_struct.sp = [2,2,subplot_it]; subplot_it = subplot_it + 1;
    b = make_grouped_bar_with_errors(input_struct);
    legend({'innocent context' 'guilty context'});
    
    %human plot
    
    %Get model data
    this_suspect_data = ...
        mbr_disconfirm_means( ...
        this_suspect_indices, ...
        [col_index_human_adj col_index_claims col_index_context]...
        );  %means
    
    this_suspect_ci = ...
        mbr_disconfirm_cis( ...
        this_suspect_indices, ...
        [col_index_human_adj col_index_claims col_index_context],2 ...
        );  %CIs
    
    this_suspect_ci = reshape( ...
        this_suspect_ci(:,1) - this_suspect_data(:,1), ...
        [2,2] ...
        );
    
    this_suspect_data = reshape( ...
        this_suspect_data(:,1), ...
        [2,2] ...
        );
    
    input_struct.input_means = this_suspect_data';
    input_struct.input_ci = this_suspect_ci';
    input_struct.title = titles{subplot_it};
    input_struct.sp = [2,2,subplot_it]; subplot_it = subplot_it + 1;
    b = make_grouped_bar_with_errors(input_struct);
    legend({'innocent context' 'guilty context'});
    
end;
%%%%%%%%%%%%%%DISCONFIRMATORY ADJUSTMENT PLOT END%%%%%%%%%%%%%%%%%%%%%%%%%%%%



disp('audi5000');
%%%%%%%%%%%%%%%%%%end, forensic_beads_prior_sim_fit_multiparam_2studies%%%%%%%%%%%%%%%%%%%%%%




%%%%%%%%%%%%%%%%%%start, get_model_ll%%%%%%%%%%%%%%%%%%%%%%
function loss = get_model_ll(free_params,params, free_params_idx, this_ps_suspect_data);

% prior = params(1);

% params_get_behav = params([1 2]);
% params_bias = params(3);
% params_noise = params(4);

params(free_params_idx) = free_params;

%this ps_suspect_data:
%now the same as raw, but adds cols col 12: seq num, col 13: model rating
this_ps_suspect_data = get_model_behaviour(params,this_ps_suspect_data);

loss = 0;
for trial = 1:size(this_ps_suspect_data,1);
    
    y_hat = this_ps_suspect_data(trial,end);   %end because model rating must be the last col    
    y = this_ps_suspect_data(trial,4);  %human / participant probability is col 4.
    
    if isnan(y_hat);
        continue
    end;
    
    %squared error loss
    loss = loss + (y_hat - y)^2;
    
    %log likelihood loss (I think)
%     ll_this_trial = y*log(y_hat+eps) + (1-y)*log(1-y_hat+eps);
%     
%     if ~isreal(ll_this_trial);
%         disp('imaginary ll');
%         fprintf('');
%     end;
%     
%     ll = ll - ll_this_trial;
    
end;    %loop through trials
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
% x = nan(nbars, ngroups);
if ngroups == 1;
    x = b.XEndPoints;
else
    for i = 1:nbars
        x(i,:) = b(i).XEndPoints;
    end
end;

% Plot the errorbars
errorbar(x',input_struct.input_means,input_struct.input_ci,'k','linestyle','none');

box off;
xlabel(input_struct.xlabel);
ylabel(input_struct.ylabel);

title(input_struct.title);
%%%%%%%%%%%%%%%%%end, make_grouped_bar_with_errors%%%%%%%%%%%%%%%%%%%%%%






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
    contexts(seq_starts(sequence,1):seq_starts(sequence,1)+10,1) = [.5; .5; this_sequence_claims_cumpro(1:9,:)]; %CONTEXT DEGREE
    
    
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









%%%%%%%%%%%%%%%%%%start, get_model_behaviour%%%%%%%%%%%%%%%%%%%%%%
function this_ps_suspect_data = get_model_behaviour(params, this_ps_suspect_data)

prior = params(1);
guilt_claim_inc = params(2);
response_bias = params(3);
response_noise = params(4);
interaction = params(5);
split = params (6);

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
        q=split;
        %         q = interaction;
        
        %get number of guilts (i.e., the number of 1s)
        ng = sum( this_ps_suspect_data(seq_start_indices(seq)+1:index,8) ) ;
        
        %get number of draws so far
        nd = claim-1  + guilt_claim_inc;
        
        %             %interaction term
        if claim == 1;
            
            interactTerm = 0;
            
        else;
            
            I_claim = 0;
            G_claim = 0;
            if this_ps_suspect_data(index,8) == 0; I_claim = 1;  %I want innocence, not guilt claims
            else; G_claim = 1;
            end;
            
            %binarised measure of context
            %this_context = this_ps_suspect_data(index,10) - .5; %-.5 if fully innocent, 0 if ambiguous, .5 is full guilty context
            %binarised measure of preceding context
            this_context = this_ps_suspect_data(index,11);
            
            interactTerm = I_claim*this_context;   %claim * preceding context interaction
%             interactTerm = I_claim*this_context - G_claim*this_context;   %claim * preceding context interaction

            %             interactTerm = this_context;   %claim * preceding context interaction

%             interactTerm = ng;   %claim * preceding context interaction
        end;
 
        %assign model probability
        %             noiseless_p = ( (1/(1 + ((1-prior)/prior)*(q/(1-q))^(nd-2*ng)))*100 ) + (interaction*ng);
        noiseless_p = ( (1/(1 + ((1-prior)/prior)*(q/(1-q))^(nd-2*ng))) ) + (interaction*interactTerm); %free parameter is weight on context*claim interaction
        %         noiseless_p = (1/(1 + ((1-prior)/prior)*(q/(1-q))^(nd-2*ng)))*100; %free parameter is weight on context*claim interaction
        
        
        %add noise and response bias
        %         this_ps_suspect_data(index,cols_this_ps_suspect_data+2) = ...
        %             response_bias + response_noise*noiseless_p;
        
        noise_p =        response_bias + response_noise*noiseless_p;
        
        if noise_p <= 0;
            noise_p = 0;
        elseif noise_p >= 1;
            noise_p = 1;
        end;
        
        this_ps_suspect_data(index,cols_this_ps_suspect_data+2) = noise_p*100;
        
        if isnan(this_ps_suspect_data(index,end));
            fprintf('');
        end;
        
        %         end;    %before first claim (sequence position 0) or a later one?
        
    end;    %loop through this sequence (claim)
    
end;    %loop through sequences
fprintf('');
%%%%%%%%%%%%%%%%%%end, simulate%%%%%%%%%%%%%%%%%%%%%%



