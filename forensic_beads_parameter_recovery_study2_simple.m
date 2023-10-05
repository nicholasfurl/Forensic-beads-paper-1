%%%%%%%%%%%%%%%%%%start, forensic_beads_parameter_recovery%%%%%%%%%%%%%%%%%%%%%%
function forensic_beads_parameter_recovery_study2_simple;

%A simpler, rewritten, debugged and double checked version of parameter recovery for
%study 2, fitting different priors to male and female sequences
%within-participant but one split, bias and noise parameter per participant.

addpath(genpath('C:\matlab_files\fiance\forensic_beads_pub_repo\Forensic-beads-paper-1\FMINSEARCHBND'));
addpath(genpath('C:\matlab_files\fiance\parameter_recovery\beta_fixed_code\Model_fitting_hybrid_study\plotSpread'));

%There's a version with fewer config'ed param levels that runs much faster
%and is suitable for debugging. Set small_list to 1 to run that. If you
%want a more thorough range of cofig'ed params, use a different nbumber
%(e.g., 0)
small_list = 0;

%1: event index,
%2:participant private id,
%3:RT,
%4: human probability
%5: sequence position (which witness is it?)
%6: suspect gender (1=female),
%7:witness gender (1=female),
%8: guilty claim (1=guilty),
%9: context (mostly innocent / mostly guilty)
data.raw = get_sub_data;

%get relevant data
stimuli = data.raw(:,2);
suspect = data.raw(:,6);
seq_pos = data.raw(:,5);
claims = data.raw(:,8);

%assemple into table
data = table(stimuli,suspect,seq_pos,claims);

%Just to ensure we have no more issues with mismatched sorting, ensure data
%are sorted by stimulus, suspect, seq_pos (preserving claim order) and then
%always loop in these orders
data = sortrows(data,{'stimuli','suspect'});


if small_list == 1;
    
    %smaller list, for debugging
    config_prior1 = [.4 .6];
    config_prior2 = [.4 .6];
    config_split = [.6 .8];
    config_response_bias = [.4 .8];
    config_response_noise = [.4 .8];
    
else;
    
    %Make lists of configured parametersp
    config_prior1 = linspace(0,1,4);
    config_prior2 = linspace(0,1,4);
    config_split = linspace(.6,.9,4);
    config_response_bias = linspace(0,1,4);
    config_response_noise = linspace(0,1,4);
    
end;    %small list or no?

lower_bounds = [0 0 .5 0 0];   %fitting will not try parameters outside these values
upper_bounds = [1 1  1 1 1];

%initialise stuff
num_params = numel(lower_bounds);   %Number of parameter variables (not combinations of their values)
num_priors1 = numel(config_prior1);
num_priors2 = numel(config_prior2);
num_splits = numel(config_split);
num_biases = numel(config_response_bias);
num_noises = numel(config_response_noise);
total_params = num_priors1*num_priors2*num_splits*num_biases*num_noises; %Number of combinations of parameter values to be tested (Not the number of parameter variables

%So there is technically only one stimulus set, that all 104 participants viewed
%But if I'm going to add random noise to probability estimates, it
%effectively creates distrinct "participants" or runs of each stim set. So
%let's keep the num_stim (counts participants as stimulus sets)
stim = unique(stimuli); %actually, participant numbers, but we treat a participant here as a stimulus set for which behaviour is simulated
num_stim = numel(stim);

%For study 2, this value should be 11 claims*8 sequences = 88 trials
num_trials_per_sub = sum(data.stimuli==stim(1));

%initialise output arrays
configured_params = NaN(num_priors1, num_priors2,num_splits, num_biases, num_noises,num_stim,num_params);
fitted_params = NaN(num_priors1, num_priors2, num_splits, num_biases, num_noises,num_stim,num_params);
configured_probabilities = NaN(num_priors1, num_priors2, num_splits, num_biases, num_noises,num_stim,num_trials_per_sub);
fitted_probabilities = NaN(num_priors1, num_priors2, num_splits, num_biases, num_noises,num_stim,num_trials_per_sub);


%Run model fitting for each parameter level for each stimulus
it = 1; %only used for disp to std out
for prior1 = 1:num_priors1;
    for prior2 = 1:num_priors2;
        for split = 1:num_splits;
            for bias = 1:num_biases;
                for noise = 1:num_noises;
                    
                    disp(sprintf('test %d of %d',it,total_params)); %only use of it
                    it = it+1;
                    
                    params(1) = config_prior1(prior1);
                    params(2) = config_prior2(prior2);
                    params(3) = config_split(split);
                    params(4) = config_response_bias(bias);
                    params(5) = config_response_noise(noise);
                    
                    %Same stim set but will have different noise term added
                    %to responses for each of the 104 runs
                    for stimulus = 1:num_stim;
                        
                        %Data for this "participant" or stimulus set. Will
                        %be the same set for every "participant", albeit
                        %reordered a little. Returns Pid (labeld "stimuli",
                        %suspect codes, seq_pos and claims for all 88
                        %trials across the 8 sequences
                        this_stim_data = data(data.stimuli==stim(stimulus),:);
                        
                        %save for analysis and comparison to fitted_params later
                        configured_params(prior1,prior2,split,bias,noise,stimulus,:) = params;
                        
                        %Need to generate behaviour separately for each
                        %sequence using the correct prior param for each
                        %get model behaviour will be called for each
                        %sequence inside of the get_behaviour_for_this_ps_sequences function call
                        configured_probabilities(prior1, prior2, split, bias, noise,stimulus,:) = ...
                            get_behaviour_for_this_ps_sequences(this_stim_data,params);
                        
                        
                        %setup dataset for fitting
                        clear data_to_fit;
                        data_to_fit.data = this_stim_data;
                        data_to_fit.configured_probabilities = ...
                            squeeze(configured_probabilities(prior1,prior2, split, bias, noise,stimulus,:));
                        
                        %now fit new parameters to the behaviour generated from these "stimuli"
                        options = optimset('MaxFunEvals',1500);
                        [new_params, ...
                            loss_temp, flag search] = ...
                            fminsearchbnd( ...
                            @(params) get_model_loss(params, data_to_fit), ...
                            params, ...
                            lower_bounds, ...  %lower parameter bounds
                            upper_bounds, ... %upper parameter bounds
                            options ...
                            );
                       
                    fitted_params(prior1,prior2,split,bias,noise,stimulus,:) = new_params;
                    
                    %Generate behavioural data from these fitted parameters
                    %get_behaviour_for* already *100 to output of get_model_behaviour
                    fitted_probabilities(prior1,prior2, split, bias, noise,stimulus,:) = ...
                         get_behaviour_for_this_ps_sequences(this_stim_data,new_params);
                        
                    end;    %stimuli / participants
                    
                end;    %priors1
            end;    %priors2
        end;    %splits
    end;    %biases
end;    %noise

save('fb_pr_s2_winSubs_ws.mat')

[R_params p_params] = correlate_output(configured_params,fitted_params);

% Create a heatmap of the param correlation matrix
f1 = figure('Color',[1 1 1]);

%parameters
% subplot(1,2,1);
hm_param = heatmap(R_params, 'Colormap', cool, 'ColorbarVisible', 'on',  'XLabel', 'Fitted parameters' , 'YLabel', 'Configured parameters');
hm_param.CellLabelFormat = '%2.2f';
ticklabels = {'Prior male', 'Prior female', 'Interaction', 'Bias', 'Noise'};
% hm_param.XDisplayLabels = ticklabels;
% hm_param.YDisplayLabels = ticklabels;
caxis([-1, 1]);


[R_probs p_probs] = corr( ...
    reshape(configured_probabilities,[],1), ...
    reshape(fitted_probabilities,[],1) ...
    );

disp(sprintf('Correlation between fitted and configured probability ratings: r = %0.4f p = %0.4f',R_probs,p_probs));


%Make plots of recovered parameter values
dims = size(fitted_params);
f2 = figure('Color',[1 1 1]);
commands = { ...
    'squeeze(fitted_params(param_level,:,:,:,:,:,param));' ...
    'squeeze(fitted_params(:,param_level,:,:,:,:,param));' ...
    'squeeze(fitted_params(:,:,param_level,:,:,:,param));' ...
    'squeeze(fitted_params(:,:,:,param_level,:,:,param));' ...
    'squeeze(fitted_params(:,:,:,:,param_level,:,param));' ...
    };
x_labels = { ...
    'config_prior1' ...
    'config_prior2' ...
    'config_split' ...
    'config_response_bias' ...
    'config_response_noise' ...
    };
for param=1:num_params
    
    subplot(1,num_params,param);
    
    these_levels = unique(eval(x_labels{param}));
    
    for param_level = 1:numel(these_levels);
        
        clear this_data;
        eval(['this_data = ' commands{param}]);
        
        this_data = reshape(this_data,[],1);
        x_data = eval([x_labels{param} '(param_level);']);
        
        sw = .05;
        handles = plotSpread(this_data ...
            ,'xValues',x_data ...
            ,'distributionColors',[.75 .75 .75] ...
            ,'distributionMarkers','.' ...
            , 'spreadWidth', sw ...
            );
        
        f_a = .1;
        these_levels_barwidth = mean(diff(these_levels))/2;
        if isnan(these_levels_barwidth); these_levels_barwidth=1; end;
        bar(x_data,x_data ...
            ,'FaceColor',[0 0 0] ...
            ,'FaceAlpha',f_a ...
            ,'EdgeColor',[0 0 0] ...
            ,'BarWidth', these_levels_barwidth ...
            );
        
    end;    %param level
    set(gca,'XTick',these_levels)
    ylabel(ticklabels{param});
    xlim([min(these_levels)-these_levels_barwidth*2 max(these_levels)+these_levels_barwidth*2])
    ylim([0 1]);
end;    %loop through the different params to make a plot of each

disp('audi5000');
%%%%%%%%%%%%%%%%%%end, forensic_beads_parameter_recovery%%%%%%%%%%%%%%%%%%%%%%









%%%%%%%%%%%%%BEGIN get_behaviour_for_this_ps_sequences%%%%%%%%%%%%%%%%
function sequence_probabilities = get_behaviour_for_this_ps_sequences(this_ps_data, params);

%Operates on a single participants' data to get simulated behaviour for
%every sequence without mixing up the order of the suspects or contexts or
%anything else


%on which indices is the display screen 0 (prior rating prompt so first rating)
seq_start_indices = find(this_ps_data.seq_pos==0);

%initialise output
sequence_probabilities = [];

%For each start index, loop through sequence and get model predictions
for seq = 1:numel(seq_start_indices);
    
    %probably unnecessary
    clear this_suspect_code this_params this_seq_data;
    
    %What suspect parameter do I need to use for this sequence?
    this_suspect_code = this_ps_data.suspect(seq_start_indices(seq));
    
    %modify param vector to pick out the prior for this sequences suspect
    this_params = [params(this_suspect_code+1) params(3:end)];
    
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













%%%%%%%%%%%%%%%%%%start, correlate_output%%%%%%%%%%%%%%%%%%%%%%
function [R P] = correlate_output(A,B);

% Define the dimensions of your arrays
dim1 = size(A, 1);
dim2 = size(A, 2);
dim3 = size(A, 3);
dim4 = size(A, 4);
dim5 = size(A, 5);
dim6 = size(A, 6);
dim7 = size(A, 7);

% Initialize the correlation matrix R
R = zeros(dim7, dim7);

% Reshape A and B into 2D matrices
A_reshaped = reshape(A, [], dim7);
B_reshaped = reshape(B, [], dim7);

% Calculate the correlation for each pair of elements in A and B
for i = 1:dim7
    for j = 1:dim7
        [R(i, j) P(i, j)] = corr(A_reshaped(:, i), B_reshaped(:, j));
    end
end
%%%%%%%%%%%%%%%%%%end, correlate_output%%%%%%%%%%%%%%%%%%%%%%











%%%%%%%%%%%%%%%%%%start, get_model_loss%%%%%%%%%%%%%%%%%%%%%%
function loss = get_model_loss(params, this_ps_data);

fitted_probabilities = ...
    get_behaviour_for_this_ps_sequences(this_ps_data.data,params);

%Loop through the trials (claims) for this suspect
loss = 0;
for trial = 1:size(this_ps_data.data,1);
    
    %squared error loss
    y_hat = fitted_probabilities(trial,end);
    y = this_ps_data.configured_probabilities(trial,1);  %human / participant probability is col 4.
    loss = loss + (y_hat - y)^2;
    
end;    %loop through trials, this suspect
%%%%%%%%%%%%%%%%%%end, get_model_loss%%%%%%%%%%%%%%%%%%%%%%




%%%%%%%%%%%%%%%%%%start, get_model_behaviour%%%%%%%%%%%%%%%%%%%%%%
function model_probabilities = get_model_behaviour(params, this_ps_suspect_data)

prior = params(1);
split = params(2);
response_bias = params(3);
response_noise = params(4);

%on which indices is the display screen 0 (prior rating prompt so first rating)
seq_start_indices = find(this_ps_suspect_data.seq_pos==0);

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
        ng = sum( this_ps_suspect_data.claims(seq_start_indices(seq)+1:index) ) ;
        
        %get number of draws so far
        nd = claim-1;
        
        %condition probability
        noiseless_p = (1/(1 + ((1-prior)/prior)*(q/(1-q))^(nd-2*ng)));
        
        %add noise and response bias
        noise_p =        response_bias + response_noise*noiseless_p;
        
        %add some Gaussian noise, using std of the residuals of model fitting
        std_resid = 73.77/100;
        noise_p = noise_p + randn(1,1)*std_resid;
        
        if noise_p <= 0;
            noise_p = 0;
        elseif noise_p >= 1;
            noise_p = 1;
        end;
        
        model_probabilities(index,1) = noise_p;
        
    end;    %loop through this sequence (claim)
    
end;    %loop through sequences
%%%%%%%%%%%%%%%%%%end, get_model_behaviour%%%%%%%%%%%%%%%%%%%%%%














%%%%%%%%%start, get_sub_data%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function data = get_sub_data(study_num_to_analyse);

data = xlsread('C:\matlab_files\fiance\forensic_beads_pub_repo\Forensic-beads-paper-1\Study2\data_trunc.xlsx');

%In raw, sequence positions 0 have NaNs in place of condition labels
%for contexts (col 9) and sometimes suspects (col 6). Put
%them back in or you'll have troubles later
nan_indices = find(data(:,5) == 0);    %find NaNs
data(nan_indices,9) = data(nan_indices+1,9);  %assign the missing values at pos 0 with the values at pos 1
data(nan_indices,6) = data(nan_indices+1,6);  %assign the missing values at pos 0 with the values at pos 1
%%%%%%%%%end, get_sub_data%%%%%%%%%%%%%%%%%%%%%%%%%%%%






