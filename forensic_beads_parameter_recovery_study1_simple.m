%%%%%%%%%%%%%%%%%%start, forensic_beads_parameter_recovery%%%%%%%%%%%%%%%%%%%%%%
function forensic_beads_parameter_recovery_study1_simple;

%forensic_beads_parameter_recovery_study1_simple.m: Im just going through
%v2 to check for bugs by simplifying / rewriting some passages. This will
%only work for Study 1 now.

%Applies stimuli from study 1 to generate behaviour for ranges of the four
%parameters and then fits model to the simulated behaviour and correlates
%configured and recovered parameters and probability / behaviour.

addpath(genpath('C:\matlab_files\fiance\forensic_beads_pub_repo\Forensic-beads-paper-1\FMINSEARCHBND'))
addpath(genpath('C:\matlab_files\fiance\parameter_recovery\beta_fixed_code\Model_fitting_hybrid_study\plotSpread'));

%1: event index,
%2:participant private id,
%3:RT,
%4: human probability
%5: sequence position (which witness is it?)
%6: suspect gender (1=female),
%7:witness gender (1=female),
%8: guilty claim (1=guilty),
%9: context (mostly innocent / mostly guilty)
data.raw = get_sub_data; %keep same as code that fits model to human behaviour

%Don't need human participant data
%We're just fitting every sequence separately here, so don't need
stimuli = data.raw(:,2);
seq_pos = data.raw(:,5);
claims = data.raw(:,8);

%claims are the sequences, seq position tells me when sequences start and end.
data = table(stimuli,seq_pos,claims);

%Make lists of configured parametersp
config_prior = linspace(0,1,4);
config_split = linspace(.6,.9,4);
config_response_bias = linspace(0,1,4);
config_response_noise = linspace(0,1,4);

%smaller list, for debugging
% config_prior = [.4 .6];
% config_split = [.6 .8];
% config_response_bias = [.4 .8];
% config_response_noise = [.4 .8];

lower_bounds = [0 .5 0 0];   %fitting will not try parameters outside these values
upper_bounds = [1  1 1 1];

%initialise stuff
it = 1; %Below I'm going to nest a loop for every parameter, so this will keep track of which iteration I'm on overall.

num_params = numel(lower_bounds);
num_priors = numel(config_prior);
num_splits = numel(config_split);
num_biases = numel(config_response_bias);
num_noises = numel(config_response_noise);
total_params = num_priors*num_splits*num_biases*num_noises;

stim = unique(stimuli); %actually, participant numbers, but we treat a participant here as a stimulus set for which behaviour is simulated
num_stim = numel(stim);

%Trial means a witness.
%in fact, for Study 1, there is one sequence
%So this should be 11 (pre-witness rating + 10 sequence positions)
%Note this assumes all participants have equal trial nums (True for Study 1)
num_trials_per_sub = sum(data.stimuli==stim(1));


%initialise output arrays
configured_params = NaN(num_priors, num_splits, num_biases, num_noises,num_stim,num_params);
fitted_params = NaN(num_priors, num_splits, num_biases, num_noises,num_stim,num_params);
configured_probabilities = NaN(num_priors, num_splits, num_biases, num_noises,num_stim,num_trials_per_sub);
fitted_probabilities = NaN(num_priors, num_splits, num_biases, num_noises,num_stim,num_trials_per_sub);


%Run model fitting for each parameter level for each stimulus
for prior = 1:num_priors;
    for split = 1:num_splits
        for bias = 1:num_biases;
            for noise = 1:num_noises;
                
                disp(sprintf('test %d of %d',it,total_params));
                it = it+1;
                
                params(1) = config_prior(prior);
                params(2) = config_split(split);
                params(3) = config_response_bias(bias);
                params(4) = config_response_noise(noise);
                
                %Get four parameters for this stimulus set / participant (Study 1, this should be only one sequence)
                for stimulus = 1:num_stim;
                    
                    %Get the data for the stimuli available to this participant
                    this_stim_data = data(data.stimuli==stim(stimulus),:);
                    
                    %save for analysis and comparison to fitted_params later
                    configured_params(prior,split,bias,noise,stimulus,:) = params;
                    
                    %Generate behavioural data from this configuration of parameters
                    %get_model_behaviour returns proportions so *100
                    configured_probabilities(prior, split, bias, noise,stimulus,:) = ...
                        get_model_behaviour(params, this_stim_data)*100;
                    
                    %setup dataset for fitting
                    clear data_to_fit;
                    data_to_fit.data = this_stim_data;
                    data_to_fit.configured_probabilities = ...
                        squeeze(configured_probabilities(prior, split, bias, noise,stimulus,:));
                    
                    %now fit new parameters to the behaviour generated from these "stimuli"
                    options = optimset('MaxFunEvals',1500);
                    [fitted_params(prior,split,bias,noise,stimulus,:), ...
                        loss_temp, flag search] = ...
                        fminsearchbnd( ...
                        @(params) get_model_loss(params, data_to_fit), ...
                        params, ...
                        lower_bounds, ...  %lower parameter bounds
                        upper_bounds, ... %upper parameter bounds
                        options ...
                        );
                    
                    %Generate behavioural data from these fitted parameters
                    fitted_probabilities(prior, split, bias, noise,stimulus,:) = ...
                        get_model_behaviour( ...
                        fitted_params(prior,split,bias,noise,stimulus,:), ...
                        data_to_fit.data)*100;
                    
                end;    %stimuli / participants
                
            end;    %priors
        end    %splits
    end;    %biases
end;    %noise

save('fb_s1_pr_simple.mat')

[R_params p_params] = correlate_output(configured_params,fitted_params);

% Create a heatmap of the param correlation matrix
f1 = figure('Color',[1 1 1]);

%parameters
% subplot(1,2,1);
hm_param = heatmap(R_params, 'Colormap', cool, 'ColorbarVisible', 'on',  'XLabel', 'Fitted parameters' , 'YLabel', 'Configured parameters');
hm_param.CellLabelFormat = '%2.2f';
ticklabels = {'Prior', 'Split', 'Bias', 'Noise'};
hm_param.XDisplayLabels = ticklabels;
hm_param.YDisplayLabels = ticklabels;
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
    'squeeze(fitted_params(param_level,:,:,:,:,param));' ...
    'squeeze(fitted_params(:,param_level,:,:,:,param));' ...
    'squeeze(fitted_params(:,:,param_level,:,:,param));' ...
    'squeeze(fitted_params(:,:,:,param_level,:,param));' ...
    };
x_labels = { ...
    'config_prior' ...
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
            ,'distributionColors',[.5 .5 .5] ...
            ,'distributionMarkers','.' ...
            , 'spreadWidth', sw ...
            );
        
        f_a = .1;
        these_levels_barwidth = mean(diff(these_levels))/2;
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










%%%%%%%%%%%%%%%%%%start, correlate_output%%%%%%%%%%%%%%%%%%%%%%
function [R P] = correlate_output(A,B);

% Define the dimensions of your arrays
dim1 = size(A, 1);
dim2 = size(A, 2);
dim3 = size(A, 3);
dim4 = size(A, 4);
dim5 = size(A, 5);
dim6 = size(A, 6);

% Initialize the correlation matrix R
R = zeros(dim6, dim6);

% Reshape A and B into 2D matrices
A_reshaped = reshape(A, [], dim6);
B_reshaped = reshape(B, [], dim6);

% Calculate the correlation for each pair of elements in A and B
for i = 1:dim6
    for j = 1:dim6
        [R(i, j) P(i, j)] = corr(A_reshaped(:, i), B_reshaped(:, j));
    end
end
%%%%%%%%%%%%%%%%%%end, correlate_output%%%%%%%%%%%%%%%%%%%%%%








%%%%%%%%%%%%%%%%%%start, get_model_loss%%%%%%%%%%%%%%%%%%%%%%
function loss = get_model_loss(params, input_data);

%model probabilities returned as proportions, so *100
fitted_probabilities = ...
    get_model_behaviour(params,input_data.data)*100;

loss = 0;
for trial = 1:size(input_data.data,1);
    
    %squared error loss
    y_hat = fitted_probabilities(trial,1);
    y = input_data.configured_probabilities(trial,1);
    loss = loss + (y_hat - y)^2;
    
end;    %loop through trials
%%%%%%%%%%%%%%%%%%end, get_model_ll%%%%%%%%%%%%%%%%%%%%%%
















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
        
        %get number of guilts so far (i.e., the number of 1s)
        ng = sum( this_ps_suspect_data.claims(seq_start_indices(seq)+1:index) ) ;
        
        %get number of draws so far
        nd = claim-1;
        
        %assign model probability
        noiseless_p = ( (1/(1 + ((1-prior)/prior)*(q/(1-q))^(nd-2*ng))) ); 

        %add noise and response bias parameters
        noise_p =        response_bias + response_noise*noiseless_p;
        
        %add some Gaussian noise, using std of the residuals of model fitting
        std_resid = 39.27;
        noise_p = noise_p + randn(1,1)*std_resid;
        
        %Make sure probability stays between 0 and 1
        if noise_p <= 0;
            noise_p = 0;
        elseif noise_p >= 1;
            noise_p = 1;
        end;
        
        %accumulate data for posterity
        %Probabilities will need to be multiplied by 100 after return from this function
        model_probabilities(index,1) = noise_p;

    end;    %loop through this sequence (claim)
end;    %loop through sequences
%%%%%%%%%%%%%%%%%%end, simulate%%%%%%%%%%%%%%%%%%%%%%











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
%%%%%%%%%end, get_sub_data%%%%%%%%%%%%%%%%%%%%%%%%%%%%






