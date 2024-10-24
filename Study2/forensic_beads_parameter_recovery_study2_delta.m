%%%%%%%%%%%%%%%%%%start, forensic_beads_parameter_recovery%%%%%%%%%%%%%%%%%%%%%%
function forensic_beads_parameter_recovery_study2_delta;

%forensic_beads_parameter_recovery_study2_delta: I modified
% forensic_beads_parameter_recovery_study2_simple to do parameter recovery
% with the delta model.

%A simpler, rewritten, debugged and double checked version of parameter recovery for
%study 2, fitting different priors to male and female sequences
%within-participant but one split, bias and noise parameter per participant.

addpath(genpath('C:\matlab_files\fiance\forensic_beads_pub_repo\Forensic-beads-paper-1\FMINSEARCHBND'));
addpath(genpath('C:\matlab_files\fiance\parameter_recovery\beta_fixed_code\Model_fitting_hybrid_study\plotSpread'));

%Settable things

%1 if parameters derived from empirical study
%2 if small test grid search
%3 if large grid search
config_list = 1;

%Same as used in Study 2
stimuli = [ ...
        NaN     1     1     0     1     1     1     0     1     1     0
        NaN     0     1     1     0     1     1     0     1     1     1
        NaN     1     0     1     0     1     1     0     1     1     1
        NaN     0     0     1     0     1     0     0     0     1     0
        NaN     0     1     0     0     0     0     1     0     1     0
        NaN     0     1     1     1     1     0     0     1     1     1
        NaN     0     0     1     0     1     0     0     1     0     0
        NaN     0     0     0     0     1     0     1     1     0     0
        ];
suspects = [0 0 1 1 0 1 1 0];    %male = 0, female = 1
num_stim = size(stimuli,1);


%List of pre-configured parameters

%1 if paranmeters derived from empirical study
%2 if small test grid search
%3 if large grid search
if config_list == 1;

    %cols: prior male, prior female, alpha (learning rate), beta (inverse temperature)
    parameters = [

         0.42109      0.41899    0.0014651       28.574
      0.50557       0.5084   0.00067853       45.996
      0.50204      0.50154   0.00040513       143.08
      0.50464      0.53516     0.076806       2.1028
      0.50554      0.50644    0.0012988       26.946
       0.4756      0.27096      0.16921       3.2426
      0.50014      0.49914   0.00018585       320.59
      0.20225      0.14763     0.012558       4.6323
      0.50007      0.50009   0.00025639       1414.3
      0.43311       0.4113      0.51541   8.8631e-13
      0.49993      0.50002   6.5438e-05       887.87
      0.34418      0.91714      0.19219   3.1915e-13
      0.41595      0.39632      0.22303      0.93723
   6.7546e-11      0.16787     0.024134       1.4174
      0.49992      0.50003    0.0001516       1417.7
      0.50005       0.4996   0.00033709       747.21
      0.51149      0.44885     0.049758       10.871
      0.50033      0.50016   0.00029028       1416.5
      0.72963      0.72923   6.7102e-12       1.2185
      0.66403      0.94621    9.527e-12       3.8864
       0.5028      0.50104   0.00029517       253.65
      0.59366      0.54209       0.1421       1.7856
      0.51068      0.48628      0.26499       2.2589
      0.50096      0.50021   0.00028878       665.99
      0.50125      0.50103   0.00059909       377.31
      0.49994      0.49996   0.00017653       1410.4
      0.50023      0.50012   0.00010526       950.36
      0.50001      0.50004    0.0001682       1417.3
      0.32813      0.18766      0.33196        3.181
      0.50086      0.50005   0.00023058       1415.5
      0.33088      0.15048     0.057236       1.7799
      0.49993      0.49986   0.00015776         1199
      0.94698   6.8057e-11       0.4353      0.37379
      0.73895       0.4471     0.072222       2.4032
      0.50074      0.49726    0.0013699       193.87
      0.87588      0.66235      0.58053   3.8584e-13
      0.50032      0.50007     0.000331       889.42
   4.8706e-10     0.026614     0.001279       4.1896
       0.4267      0.42804     0.058573       11.228
      0.61334      0.55629   4.8882e-11       1.6846
      0.50134      0.49969   0.00036649       490.74
       0.6001      0.62227   1.3643e-11      0.47262
      0.49995      0.49998        1.001       1.0947
      0.51008      0.39905     0.093843       2.7239
      0.50875      0.50139    0.0039971       74.026
      0.50391      0.44817      0.70756      0.50277
    6.006e-11     0.018113     0.039536       1.3609
      0.15189   1.8169e-11      0.16356       2.5633
      0.50008      0.49986   0.00015822         1169
      0.89083      0.86399   4.7182e-13      0.72975
       0.5072       0.4886     0.032946       11.284
        0.523      0.49432     0.010261       20.141
      0.51436      0.50201   0.00065224       27.623
      0.55508      0.56153      0.06356       6.7819
      0.50028      0.49755   0.00021243       147.79
      0.57289      0.49093      0.01533       6.2234
      0.83369      0.99903   3.4297e-13       4.3417
      0.73096      0.38507      0.55883   5.8595e-14
      0.50023      0.50126   0.00018193         1415
      0.49997      0.50001   0.00016805       404.26
            1      0.88384      0.74815      0.40519
            1      0.58466     0.085119      0.25452
      0.54987      0.53158      0.11972       1.7101
      0.19645     0.090255      0.15411        2.651
      0.67013      0.62452      0.04917       2.6438
      0.50228      0.50058   0.00034129       129.33
      0.87447      0.90861      0.58521   4.2305e-15
      0.82453      0.85995   9.2372e-11      0.73793
    0.0038757       0.0577   1.0015e-11       2.5999
      0.58377      0.46567     0.087278       2.6491
            1      0.63764      0.87172      0.51662
      0.50312      0.50155    0.0014849       93.148
      0.59431   1.4525e-09     0.024021       1.2867
      0.51251      0.50482      0.03711       5.4504
       0.5012      0.50067   0.00028137       199.79
      0.99983      0.86607   2.7787e-12      0.62925
       0.2027      0.91749      0.52562   1.3792e-14
      0.50435      0.50216     0.010079       24.263
      0.55803      0.49351     0.026843        3.518
      0.67745       0.6747    0.0099137       1.5074
    2.613e-10     0.023052    0.0038845       1.9093
      0.50062      0.50023   0.00012462       1416.4
      0.53101      0.48917      0.08122       2.7153
      0.72869      0.18931   2.6665e-09      0.85338
   3.0688e-10    0.0065048     0.074862       2.5847
      0.51251      0.48597     0.032359       12.152
      0.79913      0.78725   1.6655e-11       1.3339
     0.073653   7.5989e-13      0.14655       1.7448
      0.50014      0.50002   0.00026351       915.66
      0.56027      0.22281      0.31539       1.7825
      0.50003      0.49999   0.00026687       1071.4
      0.49996      0.49965   0.00039695       607.93
      0.54002      0.48808     0.081071         5.18
      0.50016      0.49999   0.00024883       1417.1
      0.64108      0.32457      0.25137      0.45992
      0.99978      0.89189   6.4509e-12       1.5466
      0.51835      0.40803     0.088189       3.5958
      0.55313      0.49306     0.094031       6.8617
       0.8813      0.68867      0.56965   2.0765e-13
      0.50186      0.49956    0.0010797        190.7
            1      0.73142      0.64415       3.3496
      0.91611       0.7649       0.5225   5.5397e-12
      0.54909      0.42809      0.99452       1.0472
       0.5067       0.5213      0.95118      0.13877

    ];

elseif config_list == 2;    %small grid

    %smaller list, for debugging
    malepriors = [.4 .6];
    femalepriors = [.4 .6];
    alpha = [.3 .7];
    beta = [1 100];
    
    [M, F, A, B] = ndgrid(malepriors, femalepriors, alpha, beta);

    parameters = [M(:), F(:), A(:), B(:)];


elseif config_list == 3;    %big grid

    %smaller list, for debugging
    malepriors = [0:.2:1];
    femalepriors = [0:.2:1];
    alpha = [eps:.2:1];
    beta = [1:50:1000];
    
    [M, F, A, B] = ndgrid(malepriors, femalepriors, alpha, beta);

    parameters = [M(:), F(:), A(:), B(:)];

end;

%male prior, female prior, alpha, beta
params = [0.5 0.5 0.5 12];
lower_bounds = [0 0 eps 0];   %fitting will not try parameters outside these values
upper_bounds = [1 1 Inf Inf];

%initialise stuff
num_params = numel(lower_bounds);   %Number of parameter variables (not combinations of their values)
num_combos = size(parameters,1);

%Loop through parameter combos
% it_acc = 1;
% columnNames = {'Pid', 'Sequence number', 'Suspect', 'Prior', 'Alpha', 'Beta', 'Guilt probability'};
% results_configured = array2table(zeros(0, length(columnNames)), 'VariableNames', columnNames);

for combo = 1:num_combos;

    %I consider responses to the bundle of stimuli of both genders to be a
    %"participant", to which the model will be fitted.
    for stimulus = 1:num_stim;

        if suspects(stimulus) == 0;
            these_params = [parameters(combo,1) parameters(combo,3:end)];
        elseif suspects(stimulus) == 1;
            these_params = [parameters(combo,2) parameters(combo,3:end)];
        end;

        this_sequence = stimuli(stimulus,:);

        %simulate responses to this sequence of agent using this combo of parameters
        configured_probabilities(stimulus,:) = prob_guilt_delta(these_params,this_sequence);

        %         % Assign the row vector to the first six columns
        %         results_configured(it_acc, 1:6) = array2table([combo stimulus suspects(stimulus) these_params ]);
        %
        %         % Assign the column vector to the seventh column
        %         results_configured(i, 7) = configured_probabilities;
        %
        %         it_acc = it_acc+1;
    end;    %loop through stimuli / sequences

    %Now fit the "participant" (combo) data you just created
    
    clear data_to_fit;
    data_to_fit.stimuli = stimuli;
    data_to_fit.suspects = suspects;
    data_to_fit.configured_probabilities = configured_probabilities;
    
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

    fitted_params(combo,:) = new_params;

end;    %loop through parameter combinations


[R_params p_params] = corr(parameters,fitted_params);

% Create a heatmap of the param correlation matrix
f1 = figure('Color',[1 1 1]);

%parameters
% subplot(1,2,1);
hm_param = heatmap(R_params, 'Colormap', cool, 'ColorbarVisible', 'on',  'XLabel', 'Fitted parameters' , 'YLabel', 'Configured parameters');
hm_param.CellLabelFormat = '%2.2f';
ticklabels = {'Prior male', 'Prior female', 'Alpha', 'Beta'};
hm_param.XDisplayLabels = ticklabels;
hm_param.YDisplayLabels = ticklabels;
caxis([-1, 1]);


figure; set(gcf,'Color',[1 1 1])
subplot(2,2,1);
make_param_scatter(parameters(:,1),fitted_params(:,1),ticklabels{1});
subplot(2,2,2);
make_param_scatter(parameters(:,2),fitted_params(:,2),ticklabels{2})
subplot(2,2,3);
make_param_scatter(parameters(:,3),fitted_params(:,3),ticklabels{3})
subplot(2,2,4);
make_param_scatter(parameters(:,4),fitted_params(:,4),ticklabels{4})

disp('audi5000');
%%%%%%%%%%%%%%%%%%end, forensic_beads_parameter_recovery%%%%%%%%%%%%%%%%%%%%%%



function [] = make_param_scatter(x,y, name)

% Create the scatter plot with black points
scatter(x, y, 'k', 'filled');
hold on;

% Label the axes
xlabel(sprintf('Configured: %s',name), 'FontName', 'Arial', 'FontSize', 12);
ylabel(sprintf('Fitted: %s',name), 'FontName', 'Arial', 'FontSize', 12);

% Set the font for the axes tick marks
set(gca, 'FontName', 'Arial', 'FontSize', 12);

% Get the automatic x and y limits
xlim_auto = xlim; % Get current x-axis limits
ylim_auto = ylim; % Get current y-axis limits

% Create five equally spaced tick marks on each axis, rounded to the nearest integer
xticks(round(linspace(xlim_auto(1), xlim_auto(2), 5), 2, 'significant'));
yticks(round(linspace(ylim_auto(1), ylim_auto(2), 5), 2, 'significant'));

%chance diagonal
plot([xlim_auto(1) xlim_auto(2)], [ylim_auto(1) ylim_auto(2)], 'Color', [.5 .5 .5], 'LineWidth',1)

% % Optionally, you can set the grid
% grid on;





% %%%%%%%%%%%%%BEGIN get_behaviour_for_this_ps_sequences%%%%%%%%%%%%%%%%
% function sequence_probabilities = get_behaviour_for_this_ps_sequences(this_ps_data, params);
% 
% %Operates on one sequence and 
% 
%     %What suspect parameter do I need to use for this sequence?
%     this_suspect_code = this_ps_data.suspect(seq_start_indices(seq));
%     
%     %modify param vector to pick out the prior for this sequences suspect
%     this_params = [params(this_suspect_code+1) params(3:end)];
%     
%     %Get data for this one sequence
%     this_seq_data = this_ps_data( seq_start_indices(seq):seq_start_indices(seq)+10, :);
%     
%     %Hand just the one sequence to get_model_behaviour, together with
%     %suspect-specific parameter vector and then accumulate it with the
%     %other sequences for this participant to be returned by function
%     sequence_probabilities = [ ...
%         sequence_probabilities; ...
%         get_model_behaviour(this_params,this_seq_data)*100 ...
%         ];
%     
% end;    %seq: loop through this participant's sequences
% 
% %%%%%%%%%%%%%END get_behaviour_for_this_ps_sequences%%%%%%%%%%%%%%%%













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
function loss = get_model_loss(params, data_to_fit);

%loop through sequences to get model probability predictions for whatever
%test parameters are sent here by the optimiser.
for stimulus = 1:size(data_to_fit.stimuli,1);

    %which suspect (modify parameters accordingly)
    if data_to_fit.suspects(stimulus) == 0;
        these_params = [params(1) params(3:end)];
    elseif data_to_fit.suspects(stimulus) == 1;
        these_params = [params(2) params(3:end)];
    end;

    this_sequence = data_to_fit.stimuli(stimulus,:);

    %simulate responses to this sequence of agent using this combo of parameters
    fitted_probabilities(stimulus,:) = prob_guilt_delta(these_params,this_sequence);
        
end;      %loop through sequences

%sum squared loss function
loss = sum(sum((fitted_probabilities - data_to_fit.configured_probabilities).^2));
%%%%%%%%%%%%%%%%%%end, get_model_loss%%%%%%%%%%%%%%%%%%%%%%




% %%%%%%%%%%%%%%%%%%start, get_model_behaviour%%%%%%%%%%%%%%%%%%%%%%
% function model_probabilities = get_model_behaviour(params, this_ps_suspect_data)
% 
% prior = params(1);
% split = params(2);
% response_bias = params(3);
% response_noise = params(4);
% 
% %on which indices is the display screen 0 (prior rating prompt so first rating)
% seq_start_indices = find(this_ps_suspect_data.seq_pos==0);
% 
% %initialise output
% model_probabilities = nan(size(this_ps_suspect_data,1),1);
% 
% %For each start index, loop through sequence and get model predictions
% for seq = 1:numel(seq_start_indices);
%     
%     %Loop through this sequence
%     for claim = 1:11;
%         
%         %what's the current index into this_ps_suspect_data?
%         index = seq_start_indices(seq)+claim-1;
%         
%         q=split;
%         
%         %get number of guilts (i.e., the number of 1s)
%         ng = sum( this_ps_suspect_data.claims(seq_start_indices(seq)+1:index) ) ;
%         
%         %get number of draws so far
%         nd = claim-1;
%         
%         %condition probability
%         noiseless_p = (1/(1 + ((1-prior)/prior)*(q/(1-q))^(nd-2*ng)));
%         
%         %add noise and response bias
%         noise_p =        response_bias + response_noise*noiseless_p;
%         
%         %add some Gaussian noise, using std of the residuals of model fitting
%         std_resid = 73.77/100;
%         noise_p = noise_p + randn(1,1)*std_resid;
%         
%         if noise_p <= 0;
%             noise_p = 0;
%         elseif noise_p >= 1;
%             noise_p = 1;
%         end;
%         
%         model_probabilities(index,1) = noise_p;
%         
%     end;    %loop through this sequence (claim)
%     
% end;    %loop through sequences
% %%%%%%%%%%%%%%%%%%end, get_model_behaviour%%%%%%%%%%%%%%%%%%%%%%
% 
% 
% 
% 
% 
% 








% %%%%%%%%%start, get_sub_data%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function data = get_sub_data(study_num_to_analyse);
% 
% data = xlsread('C:\matlab_files\fiance\forensic_beads_pub_repo\Forensic-beads-paper-1\Study2\data_trunc.xlsx');
% 
% %In raw, sequence positions 0 have NaNs in place of condition labels
% %for contexts (col 9) and sometimes suspects (col 6). Put
% %them back in or you'll have troubles later
% nan_indices = find(data(:,5) == 0);    %find NaNs
% data(nan_indices,9) = data(nan_indices+1,9);  %assign the missing values at pos 0 with the values at pos 1
% data(nan_indices,6) = data(nan_indices+1,6);  %assign the missing values at pos 0 with the values at pos 1
% %%%%%%%%%end, get_sub_data%%%%%%%%%%%%%%%%%%%%%%%%%%%%



%%%%%%%%%%%%%%%%%%start, prob_guilt_delta%%%%%%%%%%%%%%%%%%%%%%
function model_probabilities = prob_guilt_delta(params, this_ps_suspect_data)


prior = params(1);
alpha = params(2);
beta = params(3);

%initialise output
model_probabilities = nan(size(this_ps_suspect_data,2),1);

%For each start index, loop through sequence and get model predictions
for seq = 1:numel(this_ps_suspect_data);

    draws = this_ps_suspect_data;

    %Loop through this sequence
    for claim = 1:11;

        if claim == 1; 
            
             q_hat = prior;  %first update, without info, based just on prior
        
        else;

            q_hat = q_hat + alpha * (draws(claim) - q_hat);
        
        end;

         model_probabilities(claim,1) = exp(beta*q_hat)/(exp(beta*q_hat) + exp(beta*(1-q_hat)));    


%         prob_prelim = exp(beta*q_hat)/(exp(beta*q_hat) + exp(beta*(1-q_hat)));
% 
%         %Add a little noise to simulate participants' noisiness
%          std_resid = 12/100; %in forensic_beads_study2_2024_v3.m I took all the fitted probabilities, subtracted from them the human probabilities, took the abs of these differences, then found the std over all trials in all participants. Then I divide by 100, as that program was in percentages and here we have proportions.
%         noise_p = prob_prelim + randn(1,1)*std_resid;
% 
%          %Make sure probability stays between 0 and 1
%         if noise_p <= 0;
%             noise_p = 0;
%         elseif noise_p >= 1;
%             noise_p = 1;
%         end;
% 
%         model_probabilities(claim,1) = noise_p;
% 
    end;    %loop through this sequence (claim)
end;    %loop through sequences

%%%%%%%%%%%%%%%%%%end, guilt_prob_delta%%%%%%%%%%%%%%%%%%%%%%



