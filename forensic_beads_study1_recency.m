%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function forensic_beads_study1_recency;

%forensic_beads_study1_recency tried to model oddbaall effect by
%considering only most recent samples rather than all sequences as elapsed,
%with the recency weighting a ree parameter. Based on consultation with
%Matteo.

%forensic_beads_splitTerm_study1_simplify:
%Fits a probability estimation model to human probability estimates with
%prior, split, bias and noise parameters. Fits separate four-param models
%to each participant.

%, is forensic_beads_splitTerm_study1
%But I've done some rewriting and gone through line by line looking for bugs, improvements
%to explain why it doesn't match param recovery

addpath(genpath('C:\matlab_files\fiance\forensic_beads_pub_repo\Forensic-beads-paper-1\FMINSEARCHBND'))
addpath(genpath('C:\matlab_files\fiance\forensic_beads_pub_repo\Forensic-beads-paper-1\klabhub-bayesFactor-3d1e8a5'))
addpath(genpath('C:\matlab_files\fiance\forensic_beads_pub_repo\Forensic-beads-paper-1\shaded_plots'));

%Setable stuff!

%initial value of free params
params(1) = .5; %prior, initialised to optimal value (ground truth of paradigm)
params(2) = .6;    %split term, initialised to the ground truth / optimal value
params(3) = 0;  %decision bias, intialised to optimal value
params(4) = 1;  %decision noise, initialised to optimal value
params(5) = 1;  %learning rate, initaialied to optimal value
 
%fitted parameters constrained to be between these values
%Can be used to run simulations with fixed parameters too, instead of model
%fitting: technically all params will be free, but you can fix some or all by forcing their
%range to be one value only.
lower_bounds = [.5 .6 0 1 1];
upper_bounds = [.5 .6 0 1 1];

%If you want to simulate ideal observer performance with fixed params
% lower_bounds = [params(1) params(2) params(3) params(4)];
% upper_bounds = lower_bounds;

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
% .... And We will later on add:
%12 model probability
%13 human's adjustments
%14 model's adjustments
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

%Who are the participants?
participant_list = unique(data.Pid,'stable');   %Vitally important participant num order is maintained or it'll get mismatched with results arrays later!!!
num_participants = numel(participant_list);

%initialise table to hold modelling results
var_names = {'Pid','Suspect','Loss','Prior','Split','Bias','Noise','Learn'};
model_fitting_results = array2table(nan(0,numel(var_names)), 'VariableNames',var_names);

%initialise vector to hold model's probabilities
model_behaviour_results = [];

for participant = 1:num_participants;
    
    clear this_ps_data new_params new_loss;
    
    %get data for this participant
    this_ps_data = ...
        data( ...
        data.Pid == participant_list(participant),...
        :);
 
        disp(sprintf('fitting participant %d', participant));
        
        %fit params to data.
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
        %We can just fit all data for each participant and take suspect
        %from first row of participant data here because in study 1, we
        %know every participant sees one suspect.
        model_fitting_results = ...
            [model_fitting_results; ...
            num2cell([participant_list(participant) this_ps_data.Suspect(1) new_loss new_params])
            ];
        
        %get performance for this model
        model_behaviour_results = [ ...
            model_behaviour_results; ...
            get_model_behaviour(new_params,this_ps_data)*100 ...
            ];
        
%     end;    %loop through suspects
    
end;    %loop through participants

save('fb_fit_study1_simpler.m');

%Add the accumulated model probabilities as column of data table
data.ModelProbability = model_behaviour_results;

%add adjustment columns to data table
data = [ ...
    data, ...
    array2table(get_adjustments(data),'VariableNames',{'HumanAdjust','ModelAdjust'}) ...
    ];

%Stats tests comparing effects of suspect on prior parameter

suspect_0 = model_fitting_results.Prior(model_fitting_results.Suspect==0);
suspect_1 = model_fitting_results.Prior(model_fitting_results.Suspect==1);

%trad t-test
[th tpvals tci tstats] =  ttest2(suspect_0, suspect_1);

%paired Bayesian test
[bf10,bfpvals,bfci,bfstats] = bf.ttest2( suspect_0, suspect_1);

disp(sprintf('Suspect effect on prior: trad pval: %0.2f bf10: %7.4f',tpvals,bf10));

%Bar plots of parameter values and loss
make_parameter_plot(model_fitting_results);

%Kernel distribution plot of the prior parameter values
make_kernel_dists(model_fitting_results);

%Plot model probabilities by sequence position
plot_sequence_behaviour(data);

%Plot adjustments by suspect and context
plot_adjustments(data);

disp('audi5000');
%%%%%%%%%%%%%%%%%%end, forensic_beads_prior_sim_fit_multiparam_2studies%%%%%%%%%%%%%%%%%%%%%%














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
    input_struct.ylim = [-15 15];
    
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
    input_struct.ylim = [-15 15];
    
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

% %I wouldn't need to use two steps where I group by participant first for
% %study 1, but I'd need to if I wanted to do Study 2.
temp = grpstats(data,{'Suspect', 'Context','SequencePosition','Pid'},{'mean'});

%Now collapse over participant too, getting ci's over P's as we go
means = grpstats(temp,{'Suspect', 'Context','SequencePosition'},{'mean' 'meanci'});

suspects = unique(means.Suspect);
contexts = unique(means.Context);

cmap_it = 1;
for context = 1:numel(contexts);
    
    subplot(1,2,context);
    
    for suspect = 1:numel(suspects);
        
        %data for this suspect / context combo
        this_data = ...
            means.mean_mean_ModelProbability( means.Suspect==suspects(suspect) & means.Context==contexts(context) );
        
        %cis for this suspect / context combo
        this_cis = ...
            means.meanci_mean_ModelProbability( means.Suspect==suspects(suspect) & means.Context==contexts(context), 2 ) - this_data;
        
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
        
        text(2,30-suspect*6,legend_labels{suspect},'Color',cmap(suspect,:));
        
        cmap_it = cmap_it + 1;
        
    end;    %suspects loop
    
    title(sprintf('Context: %d',contexts(context)));
    ylim([1 100]);
    box off;
    xlabel('Sequence position');
    ylabel('Probability');
    
end;    %contexts loop

%%%%%%%%%%%%%%SEQUENCE POSITION PLOT ENDS%%%%%%%%%%%%%%%%%%%%%%%%%%%%%








%%%%%%%%%%%%%%MAKE KERNEL DISTS BEGIN%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function make_kernel_dists(model_fitting_results);

cmap = [0 0 0;
    .25 .25 .25
    ];


alphas = [0 .75];

h1 = figure('Color',[1 1 1]);

suspects = unique(model_fitting_results.Suspect);
suspects_num = numel(suspects);

for suspect = 1:suspects_num;  %suspect parameters loop
    
    %Kernel density
    [kernel(:,suspect) xi(:,suspect)] = ...
        ksdensity( ...
        model_fitting_results.Prior(model_fitting_results.Suspect==suspects(suspect)) ...
        );
    
    plot_shaded(xi(:,suspect),kernel(:,suspect),'Alpha',alphas(suspect),'Color',cmap(suspect,:),'LineWidth',2);
    hold on;
    
end;    %suspect parameters loop

box off
xlabel('Estimated prior parameter value');
ylabel('Probability (kernel) density');
%%%%%%%%%%%%%%MAKE KERNEL DISTS ENDS%%%%%%%%%%%%%%%%%%%%%%%%%%%%%





%%%%%%%%%%%%%%MAKE_PARAMETER_PLOT BEGIN%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function make_parameter_plot(model_fitting_results);

%average by suspect code:
param_means = grpstats(model_fitting_results, 'Suspect',{'mean', 'meanci'});

%plot
h1 = figure('Color',[1 1 1]);
input_struct.fig = h1;

%Put parameter values in first subplot
input_struct.sp = [1,2,1];

input_struct.input_means = ...
    [param_means.mean_Prior',
    param_means.mean_Split',
    param_means.mean_Bias',
    param_means.mean_Noise'];

input_struct.input_ci = ...
    [param_means.meanci_Prior(:,2)'
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
input_struct.ylim = [0 5000];

b = make_grouped_bar_with_errors(input_struct);

%While I'm at it, get the grand standard deviation of the residuals for use
%in parameter recovery

%average all rows (sequences / one sequence per participant)
model_fitting_results.residual = sqrt( model_fitting_results.Loss );
model_fitting_results.Ones = ones(size(model_fitting_results,1),1);
stds = grpstats(model_fitting_results, 'Ones',{'std'});

disp(sprintf('std of all residuals is: %4.4f',stds.std_residual));


%%%%%%%%%%%%%%MAKE_PARAMETER_PLOT ENDS%%%%%%%%%%%%%%%%%%%%%%%%%%%%%







%%%%%%%%%%%%%%%%%%start, get_model_loss%%%%%%%%%%%%%%%%%%%%%%
function loss = get_model_loss(params, this_ps_suspect_data);

fitted_probabilities = get_model_behaviour(params,this_ps_suspect_data)*100;

loss = 0;
for trial = 1:size(this_ps_suspect_data,1);
    
    %squared error loss
    y_hat = fitted_probabilities(trial,end);
    y = this_ps_suspect_data.HumanProbability(trial,1);  %human / participant probability is col 4.
    loss = loss + (y_hat - y)^2;
    
end;    %loop through trials
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
ylim([input_struct.ylim(1) input_struct.ylim(2)]);

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
    contexts(seq_starts(sequence,1):seq_starts(sequence,1)+10,1) = [NaN; NaN; this_sequence_claims_cumpro(1:9,:)]; %CONTEXT DEGREE
    
    
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




