%%%%%%%%%%%%%%%%%%start, forensic_beads_prior_sim%%%%%%%%%%%%%%%%%%%%%%
function forensic_beads_simplify_study2;

%Fits a probability estimation model to human probability estimates with
%prior, split, bias and noise parameters. For each participant, fits split
%noise and bias plus fits a prior to each suspect within participant.

addpath(genpath('C:\matlab_files\fiance\forensic_beads_pub_repo\Forensic-beads-paper-1\FMINSEARCHBND'))
addpath(genpath('C:\matlab_files\fiance\forensic_beads_pub_repo\Forensic-beads-paper-1\klabhub-bayesFactor-3d1e8a5'))
addpath(genpath('C:\matlab_files\fiance\forensic_beads_pub_repo\Forensic-beads-paper-1\shaded_plots'))

%%%%%%
%Set-able stuff!

ideal_observer = 0;    %set tpo one if you want ideal observer; set to something else (like 0) if you want to free the four parameters

%initial value of free params
params(1) = .5; %prior, initialised to optimal value (ground truth of paradigm)
params(2) = .5; %prior, initialised to optimal value (ground truth of paradigm)
params(3) = .7;  %split term, initialised to optimal value
params(4) = 0;  %bias term, intialised to optimal value
params(5) = 1;  %noise term, initialised to optimal value

if ideal_observer == 1;
    
    %If you wantr to simulate ideal observer performance with fixed params
    lower_bounds = [params(1) params(2) params(3) params(4) params(5)];
    upper_bounds =lower_bounds;
    
else;
    
    %fitted parameters constrained to be between these values
    %technically all params will be free, but you can fix some by forcing their
    %range to be one value only.
    lower_bounds = [0 0 0 0 .5];
    upper_bounds = [1 1 1 1  1];
    
end;


%%%%%%%%
%Now get data ....

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
var_names = ...
    {'Event',
    'Pid',
    'RT',
    'HumanProbability',
    'SequencePosition',
    'Suspect',  %1 = prejudice group
    'WitnessGender',    %1=female
    'Claim',    %1=guilty
    'Context',  %1=Mostly Guilty
    'ContextDegree',    %Higher values, more preceding guilt claims
    'ContextCat'   %1=more guilt claims
    };
data = ...   %assign data to struct
    array2table( ...    %convert imported data to table
    get_sub_data, ... %import data
    'VariableNames', var_names ...
    );

data = sortrows(data,{'Pid','Suspect','Context'});

%Who are the participants?
participant_list = unique(data.Pid,'stable');    %Vitally important participant num order is maintained or it'll get mismatched with results arrays later!!!
num_participants = numel(participant_list);

%initialise data matrix to hold modelling results
var_names = {'Pid','Loss','PriorMale','PriorFemale','Split','Bias','Noise'};

%initialise table to hold modelling results
model_fitting_results = array2table(nan(0,numel(var_names)), 'VariableNames',var_names);

%initialise vector to hold models' probabilities
model_behaviour_results = [];

for participant = 1:num_participants;
    
    clear this_ps_data this_ps_suspect_codes this_ps_num_suspects;
    
    %get data for this participant
    this_ps_data = ...
        data( ...
        data.Pid == participant_list(participant),...
        :);
    
    disp(sprintf('fitting participant %d', participant))
    
    %fit params to data
    options = optimset('MaxFunEvals',10000);
    [new_params, new_loss, flag search] = ...
        fminsearchbnd( ...
        @(params) get_model_loss(params, this_ps_data), ...
        params, ...
        lower_bounds, ...  %lower parameter bounds
        upper_bounds, ... %upper parameter bounds
        options ...
        );
    
    %save parameter fitting output (see line before loop for col names)
    model_fitting_results = ...
        [model_fitting_results; ...
        num2cell([participant_list(participant) new_loss new_params])
        ];
    
    
    %get performance for this model
    %This needs to be divided by suspect, once for each of the
    %suspect-specific parameters,both of which are estimated on this loop
    
    %Accumulate results for plotting. Need to re-loop suspects if Study 2
    model_behaviour_results = [ ...
        model_behaviour_results;
        get_behaviour_for_this_ps_sequences(this_ps_data,new_params) ...
        ];

end;    %loop through participants

%Add the accumulated model probabilities as column of data table
data.ModelProbability = model_behaviour_results;

%add adjustment columns to data table
data = [ ...
    data, ...
    array2table(get_adjustments(data),'VariableNames',{'HumanAdjust','ModelAdjust'}) ...
    ];


%Run quick t-test to see if suspect gender has significant effect on prior

%get (within-participant) differences between the two prior parameters
prior_diffs = model_fitting_results.PriorMale - model_fitting_results.PriorFemale;

%trad paired t-test
[th tpvals tci tstats] =   ...   
    ttest(prior_diffs);

%paired Bayesian test
[bf10,bfpvals,bfci,bfstats] = ...
        bf.ttest( prior_diffs);

disp(sprintf('Suspect effect on prior: trad pval: %0.2f bf10: %7.0f',tpvals,bf10));


%Bar plots of parameter values and loss
make_parameter_plot(model_fitting_results);

%Kernel distribution plot of the prior parameter values
make_kernel_dists(model_fitting_results);

%Plot model probabilities by sequence position
plot_sequence_behaviour(data);

%Plot adjustments by suspect and context
plot_adjustments(data);

disp('audi5000');
%%%%%%%%%%%%%%%%%%end, forensic_beads_simplify_study2%%%%%%%%%%%%%%%%%%%%%%






%%%%%%%%%%%%%BEGIN get_behaviour_for_this_ps_sequences%%%%%%%%%%%%%%%%
function sequence_probabilities = get_behaviour_for_this_ps_sequences(this_ps_data,new_params);

%Operates on a single participants' data to get simulated behaviour for
%every sequence without mixing up the order of the suspects or contexts or
%anything else


%on which indices is the display screen 0 (prior rating prompt so first rating)
seq_start_indices = find(this_ps_data.SequencePosition==0);

%initialise output
sequence_probabilities = [];

%For each start index, loop through sequence and get model predictions
for seq = 1:numel(seq_start_indices);
    
    %probably unnecessary
    clear this_suspect_code this_params this_seq_data;

    %What suspect parameter do I need to use for this sequence?
    this_suspect_code = this_ps_data.Suspect(seq_start_indices(seq));
    
    %modify param vector to pick out the prior for this sequences suspect
    this_params = [new_params(this_suspect_code+1) new_params(3:end)];
    
    %Get data for this one sequence
    this_seq_data = this_ps_data( seq_start_indices(seq):seq_start_indices(seq)+10, :);
    
    %Hand just the one sequence to get_model_behaviour, together with
    %suspect-specific parameter vector and then accumulate it with the
    %other sequences for this participant to be returned by function
    sequence_probabilities = [ ...
        sequence_probabilities; ...
        get_model_behaviour(this_params,this_seq_data)*100 ...
        ];
   
end;    %seq: loop through this participant's sequences

%%%%%%%%%%%%%END get_behaviour_for_this_ps_sequences%%%%%%%%%%%%%%%%








%%%%%%%%%%%%%%DISCONFIRMATORY ADJUSTMENT PLOT BEGIN%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function plot_adjustments(data);


%I wouldn't need to use two steps where I group by participant first for
%study 1, but I'd need to if I wanted to do Study 2 and want cis over Ps and not trials.
temp = grpstats(data,{'Claim', 'ContextCat','Suspect','Pid'},{'mean'});

%Now collapse over participant too, getting ci's over P's as we go
means = grpstats(temp,{'Claim', 'ContextCat','Suspect'},{'mean' 'meanci'});

h2 = figure('Color',[1 1 1]);

suspects = unique(means.Suspect);
suspects_num = numel(suspects);

%make suplots for each suspect
subplot_it = 1;
titles = {'Model suspect 0', 'Human suspect 0' 'Model suspect 1', 'Human suspect 1' };
for suspect = 1:suspects_num;

    %indices for getting plot data (just to be organised)
    this_suspect_indices = means.Suspect==suspects(suspect);
    
    
    
    %MODEL PLOT
    
    clear this_suspect_data this_suspect_ci input_struct;
    
    %Get model adjustment means for this suspect
    %should be in order:
    %claim context
    %0 0
    %0 1
    %1 0
    %1 1
    this_suspect_data = ...
        means.mean_mean_ModelAdjust(this_suspect_indices);
    
    this_suspect_ci = ...
        means.meanci_mean_ModelAdjust(this_suspect_indices,2) - this_suspect_data;
    
    %After reshape should be claims in cols and contexts in rows
    this_suspect_data = reshape( ...
        this_suspect_data(:,1), ...
        [2,2] ...
        );
    
    this_suspect_ci = reshape( ...
        this_suspect_ci, ...
        [2,2] ...
        );
    
    %setup plot
    input_struct.fig = h2;
    
    input_struct.sp = [2,2,subplot_it]; 
    input_struct.title = titles{subplot_it};
    subplot_it = subplot_it + 1;
    
    input_struct.input_means = this_suspect_data';
    input_struct.input_ci = this_suspect_ci';
    input_struct.xlabel = 'Claim Type';
    input_struct.ylabel = 'Adjustment';
    input_struct.ylim = [-16 16];
    
    b = make_grouped_bar_with_errors(input_struct);
    legend({'innocent context' 'guilty context'});
 
    %HUMAN PLOT
    
    clear this_suspect_data this_suspect_ci input_struct;
    
    %Get model adjustment means for this suspect
    %should be in order:
    %claim context
    %0 0
    %0 1
    %1 0
    %1 1
    this_suspect_data = ...
        means.mean_mean_HumanAdjust(this_suspect_indices);
    
    this_suspect_ci = ...
        means.meanci_mean_HumanAdjust(this_suspect_indices,2) - this_suspect_data;
    
    %After reshape should be claims in cols and contexts in rows
    this_suspect_data = reshape( ...
        this_suspect_data(:,1), ...
        [2,2] ...
        );
    
    this_suspect_ci = reshape( ...
        this_suspect_ci, ...
        [2,2] ...
        );
    
    %setup plot
    input_struct.fig = h2;
    
    input_struct.sp = [2,2,subplot_it]; 
    input_struct.title = titles{subplot_it};
    subplot_it = subplot_it + 1;
    
    input_struct.input_means = this_suspect_data';
    input_struct.input_ci = this_suspect_ci';
    input_struct.xlabel = 'Claim Type';
    input_struct.ylabel = 'Adjustment';
    input_struct.ylim = [-16 16];
    
    b = make_grouped_bar_with_errors(input_struct);
    legend({'innocent context' 'guilty context'});
   
end;    %Loop through suspects
%%%%%%%%%%%%%%DISCONFIRMATORY ADJUSTMENT PLOT ENDS%%%%%%%%%%%%%%%%%%%%%%%%%%%%










%%%%%%%%%%%%%%SEQUENCE POSITION PLOT BEGIN%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function plot_sequence_behaviour(data);

%boring plot throatclearing stuff
cmap = [0 0 0;
    .5 .5 .5;
    ];

figure('Color',[1 1 1]);

legend_labels = {'Suspect 0' 'Suspect 1'};


%get data to plot

%I wouldn't need to use two steps where I group by participant first for
%study 1, but I'd need to if I wanted to do Study 2.
temp = grpstats(data,{'Suspect', 'Context','SequencePosition','Pid'},{'mean'});

%Now collapse over participant too, getting ci's over P's as we go
means = grpstats(temp,{'Suspect', 'Context','SequencePosition'},{'mean' 'meanci'});
% 
% suspects = unique(means.Suspect,'sorted');
% contexts = unique(means.Context,'sorted');

cmap_it = 1;
for context = [1 2];
    
    subplot(1,2,context);
    
    for suspect = [1 2];
        
        %data for this suspect / context combo
        this_data = ...
            means.mean_mean_ModelProbability( means.Suspect==suspect-1 & means.Context==context-1 );
        
        %cis for this suspect / context combo
        this_cis = ...
            means.meanci_mean_ModelProbability( means.Suspect==suspect-1 & means.Context==context-1, 2 ) - this_data;
        
        errorbar( ...
            this_data, ...
            this_cis, ...
            'Color', cmap(suspect,:) ...
            );
        hold on;
        plot( ...
            this_data, ...
            'Marker','o', ...
            'Color', cmap(suspect,:) ...
            );
        
        hold on;
        
        text(8,60-suspect*6,legend_labels{suspect},'Color',cmap(suspect,:));
        
        cmap_it = cmap_it + 1;
        
    end;    %suspects loop
    
    title(sprintf('Context: %d',context-1));
    ylim([1 100]);
    box off;
    xlabel('Sequence position');
    ylabel('Probability');
    
end;    %contexts loop
%%%%%%%%%%%%%%SEQUENCE POSITION PLOT END%%%%%%%%%%%%%%%%%%%%%%%%%%%%%









%%%%%%%%%%%%%%MAKE_PARAMETER_PLOT BEGIN%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function make_parameter_plot(model_fitting_results);

model_fitting_results.ones = ones(size(model_fitting_results,1),1);

%Every row is a participant already and suspects are in cols, so just get mean and ci for whole column
param_means = grpstats(model_fitting_results, 'ones',{'mean', 'meanci'});

%plot
h1 = figure('Color',[1 1 1]);
input_struct.fig = h1;

%Put parameter values in first subplot
input_struct.sp = [1,2,1];

input_struct.input_means = ...
    [param_means.mean_PriorMale',
    param_means.mean_PriorFemale',
    param_means.mean_Split',
    param_means.mean_Bias',
    param_means.mean_Noise'];

input_struct.input_ci = ...
    [param_means.meanci_PriorMale(:,2)'
    param_means.meanci_PriorFemale(:,2)',
    param_means.meanci_Split(:,2)'
    param_means.meanci_Bias(:,2)'
    param_means.meanci_Noise(:,2)'] - ...
    input_struct.input_means;

input_struct.xlabel = 'Parameter';
input_struct.ylabel = 'Mean Parameter value';
input_struct.title = 'Parameter values';
input_struct.ylim = [0 1];
b = make_grouped_bar_with_errors(input_struct);

clear input_struct; %just in case ...


%loss plot
input_struct.fig = h1;
input_struct.sp = [1,2,2];

input_struct.input_means = ...
    param_means.mean_Loss;

input_struct.input_ci = ...
    param_means.meanci_Loss(:,2) - ...
    input_struct.input_means;

input_struct.xlabel = '';
input_struct.ylabel = 'Squared error (loss)';
input_struct.title = 'Model fit';
input_struct.ylim = [0 150000];

b = make_grouped_bar_with_errors(input_struct);
%%%%%%%%%%%%%%MAKE_PARAMETER_PLOT ENDS%%%%%%%%%%%%%%%%%%%%%%%%%%%%%








%%%%%%%%%%%%%%MAKE KERNEL DISTS BEGIN%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function make_kernel_dists(model_fitting_results);

cmap = [0 0 0;
    .25 .25 .25
    ];


alphas = [0 .75];

h1 = figure('Color',[1 1 1]);

% suspects = unique(model_fitting_results.Suspect);
% suspects_num = numel(suspects);

for suspect = [1, 2];  %suspect parameters loop
    
    %Kernel density
    [kernel(:,suspect) xi(:,suspect)] = ...
        ksdensity( ...
        table2array(model_fitting_results(:,suspect+2)) ...
        );
    
    plot_shaded(xi(:,suspect),kernel(:,suspect),'Alpha',alphas(suspect),'Color',cmap(suspect,:),'LineWidth',2);
    hold on;
    
end;    %suspect parameters loop

box off
xlabel('Estimated prior parameter value');
ylabel('Probability (kernel) density');
%%%%%%%%%%%%%%MAKE KERNEL DISTS ENDS%%%%%%%%%%%%%%%%%%%%%%%%%%%%%











%%%%%%%%%%%%%%%%%%start, get_model_loss%%%%%%%%%%%%%%%%%%%%%%
function loss = get_model_loss(params, this_ps_data);

fitted_probabilities = ...
    get_behaviour_for_this_ps_sequences(this_ps_data,params);


% %Get suspects (Should return two suspects: 0 and 1
% this_ps_suspect_codes = unique(rmmissing(this_ps_data.Suspect));     %get suspect codes present in this participant
% this_ps_num_suspects = numel(this_ps_suspect_codes);              %How many codes for this participant?
% 
% loss = 0;
% %Now loop through suspects
% for suspect = 1:this_ps_num_suspects;
%     
%     %Get the probability rating data to fit for this suspect in this participant
%     this_ps_suspect_data = ...
%         this_ps_data(this_ps_data.Suspect == this_ps_suspect_codes(suspect),:);
%     
%     %alter parameter list to be specific for this suspect
%     this_params = [params(suspect) params(3:end)];   %The current free parameter value of the suspect prior tested in this suspect loop iteration and the other current values of parameters
%     
%     %Behaviour for this suspect, using suspect-specific prior and the other parameters shared between participants
%     %Returns proportions so multiply by 100 to compare against human probabilities
%     fitted_probabilities = ...
%         get_model_behaviour(this_params,this_ps_suspect_data)*100;
    
    %Loop through the trials (claims) for this suspect
    loss = 0;
    for trial = 1:size(this_ps_data,1);
        
        %squared error loss, accumulating over both suspects
        y_hat = fitted_probabilities(trial,end);
        y = this_ps_data.HumanProbability(trial,1);  %human / participant probability is col 4.
        loss = loss + (y_hat - y)^2;
        
    end;    %loop through trials, this suspect
    
% end;    %loop through the two suspects
%%%%%%%%%%%%%%%%%%end, get_model_loss%%%%%%%%%%%%%%%%%%%%%%











%%%%%%%%%%%%%%%%%%start, get adjustments%%%%%%%%%%%%%%%%%%%%%%
function adjustment_cols = get_adjustments(data);

%on which indices is the display screen 0 (prior rating prompt so first rating)
seq_start_indices = find(data.SequencePosition==0);

%initialise (human and model adjustment cols)
adjustment_cols = nan(size(data,1),2);

%For each start index, loop through sequence
for seq = 1:numel(seq_start_indices);
    
    %Loop through this sequence
    for claim = 2:11;   %ignore "first" claim in this sequence, with no preceding claim (so no adjustment)
        
        %what's the current index into data?
        index = seq_start_indices(seq)+claim-1;
        
        %human adjustments
        adjustment_cols(index,1) = ...
            data.HumanProbability(index) -  data.HumanProbability(index-1);
        
        %model adjustments
        adjustment_cols(index,2) = ...
            data.ModelProbability(index) -  data.ModelProbability(index-1);
        
    end;    %claims
    
end;    %sequences
%%%%%%%%%%%%%%%%%%end, get adjustments%%%%%%%%%%%%%%%%%%%%%%







%%%%%%%%%%%%%%%%%begin, make_grouped_bar_with_errors%%%%%%%%%%%%%%%%%%%%%%
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
ylim([input_struct.ylim(1) input_struct.ylim(2)]);

title(input_struct.title);
%%%%%%%%%%%%%%%%%end, make_grouped_bar_with_errors%%%%%%%%%%%%%%%%%%%%%%







%%%%%%%%%%%%%%%%%%start, get_context%%%%%%%%%%%%%%%%%%%%%%
function contexts = get_contexts(raw);

%This next loop takes out each sequence, finds the cumulative sum of guilt
%claims for each sequence position and puts it in 9th col of raw, and finds
%cumulative proportion of guilt claims and puts it in 10th column of raw

degree_half = NaN;  %When computing context degree, do I want undefined contexts (i.e., first and second ratings) to be treated as ambiguous contexts (50/50 guilt/innocent) to be removed from analysis (NaN)?

seq_starts = find(raw(:,5) == 0);   %first element of each sequence

for sequence=1:size(seq_starts,1); %for every start of a sequence
    
    clear this_sequence_claims* temp_MG;
    
    %extract claims for this sequence (positions 1 through 10)
    this_sequence_claims = raw(seq_starts(sequence,1)+1:seq_starts(sequence,1)+10,8);
    %get cumsum of guilt claims (=1 each)
    this_sequence_claims_cumsum = cumsum( this_sequence_claims);
    this_sequence_claims_cumpro = this_sequence_claims_cumsum./[1:10]';
    %save proportions to array called raw
    contexts(seq_starts(sequence,1):seq_starts(sequence,1)+10,1) = ...
        [degree_half; degree_half; this_sequence_claims_cumpro(1:9,:)]; %CONTEXT DEGREE
      
    %%%Now do stuff that depends on seq position. Use clunky if/then to
    %%%find if each position's cum proportion is classified as mostly
    %%%guilty (=1) or innocent (=0) or neither (NaN)(raw col 10). 
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
function data = get_sub_data;

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
function model_probabilities = get_model_behaviour(params, this_ps_suspect_data)

prior = params(1);
split = params(2);
response_bias = params(3);
response_noise = params(4);

%on which indices is the display screen 0 (prior rating prompt so first rating)
seq_start_indices = find(this_ps_suspect_data.SequencePosition==0);

%initialise output
model_probabilities = nan(size(this_ps_suspect_data,1),1);

%For each start index, loop through sequence and get model predictions
for seq = 1:numel(seq_start_indices);
    
    %Loop through this sequence
    for claim = 1:11;
        
        %what's the current index into this_ps_suspect_data?
        index = seq_start_indices(seq)+claim-1;
        
        q=split;
        
        %get number of guilts (i.e., the number of 1s)
        ng = sum( this_ps_suspect_data.Claim(seq_start_indices(seq)+1:index) ) ;
        
        %get number of draws so far
        nd = claim-1;
        
        %condition probability
        noiseless_p = (1/(1 + ((1-prior)/prior)*(q/(1-q))^(nd-2*ng)));
        
        %add noise and response bias
        noise_p =        response_bias + response_noise*noiseless_p;
        
        if noise_p <= 0;
            noise_p = 0;
        elseif noise_p >= 1;
            noise_p = 1;
        end;
        
        model_probabilities(index,1) = noise_p;
        
    end;    %loop through this sequence (claim)
    
end;    %loop through sequences
%%%%%%%%%%%%%%%%%%end, get_model_behaviour%%%%%%%%%%%%%%%%%%%%%%

