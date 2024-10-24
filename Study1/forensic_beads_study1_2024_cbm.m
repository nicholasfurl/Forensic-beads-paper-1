function [sus_means_p_ss sus_key_p_ss] = forensic_beads_study1_2024_cbm;

%forensic_beads_study1_2024_cbm is a rewrite of forensic_beads_study1_2024
%but intended to test if I can get good model fits / param recovery using
%cbm toolbox instead of fminsearchbnd.m. 

%I'm now adapting forensic_beads_study2_2024_v3.m to be used for study 1 by
%integrating in elements from forensic_beads_splitTerm_study1_simplify,

%forensic_beads_study2_2024_v3: So it appears that primacy and recency /
%delta can lead to oddball bias. I am adding model fitting and comparison
%functionality for these models.

%forensic_beads_study2_2024_v2: OK I radically chanmged things. forensic_beads_study2_2024 (i.e., v1) is
%now obsolete and its functionally has been moved here, which is now enturely self-contained.
%All comments below refer only to forensic_beads_study2_prior_model or its
%predecessors, designed just to run the prior model. Now it's code has been
%expanded here to make study 2 figures in response to JEP:LMC reviews.

%sus_key_p_ss is also returned, which is (rows: pids*suspect*context, so rows
% match those in sus_means_p_ss) by 3 (codes for pid, then suspect, then
% context). This is used for grouping the data when running ttests and
% averaging. It should be identitical and parallel to sus_key_ss for human participants
% in forensic_beads_study2_2024.

%I've modified this code to now only return sus_means_p_ss, which takes a
%(rows: pid*suspect*context) * (cols: seq pos) format, similar to one used
%to compute ttests and plot by seq position in forensic_beads_study2_2024
%for humans and models. So here we'll return a similar format into that
%code to be treated similarly to those.

%returns sus_cis_p, which is like sus_means_p, except it's 2*11*2 with the
%third dimension capturing the upper (1) and lower (2) bounds of the ci, instead of just the mean value.

%returns sus_means_p, which is the 2*11 matrix of prior mean model probabilities with sequence
%positions in columns and suspect*context combinations in rows (male suspect = 0,
%innocent context = 0; row1: 0 0; row2: 0 1; row3: 1 0; row4: 1 1).

%Modified from forensic_beads_simplify_study2_split2Suspects.m to be called
%by forensic_beads_study2_2024 and to support its graphing.

%In forensic_beads_simplify_study2_split2Suspects, prior and split vary
%by suspect but the two response parameters do not

%Fits a probability estimation model to human probability estimates with
%prior, split, bias and noise parameters. For each participant, fits split
%noise and bias plus fits a prior to each suspect within participant.

tic

%addpath(genpath('C:\matlab_files\fiance\forensic_beads_pub_repo\Forensic-beads-paper-1\FMINSEARCHBND'))
addpath(genpath('C:\matlab_files\fiance\forensic_beads_pub_repo\Forensic-beads-paper-1\klabhub-bayesFactor-3d1e8a5'))
addpath(genpath('C:\matlab_files\fiance\forensic_beads_pub_repo\Forensic-beads-paper-1\shaded_plots'))
addpath(genpath('C:\matlab_files\fiance\parameter_recovery\beta_fixed_code\Model_fitting_hybrid_study\plotSpread'));
addpath(genpath('C:\matlab_files\cbm-master'));

skip_models = 1;    %will run models but not plot them


%%%%%%%%
%1: event index,
%2:participant private id,
%3:RT,
%4: human probability
%5: sequence position (which witness is it?)
%6: suspect religion (1=Christian),
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
    'Suspect',  %suspect religion (0=Atheist)
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
participant_list = unique(data.Pid,'stable');    %Vitally important participant num order is maintained or it'll get mismatched with results arrays later!!!
participant_list = participant_list(1);
num_participants = numel(participant_list);

%Before I do fitting, I need to construct cell array of participant data to
%deliver to cbm_lap
% P_data = cell(num_participants,1);
p_it_a = 1;
p_it_c = 1;
for participant = 1:num_participants;

    clear this_ps_data temp;

    %get data for this participant
    temp = ...
        data( ...
        data.Pid == participant_list(participant),...
        :);

    if temp.Suspect(1) == 0;    %atheist
    
        P_data_A{p_it_a} = temp;
        p_it_a = p_it_a + 1;

    else;   %Christian

        P_data_C{p_it_c} = temp;
        p_it_c = p_it_c + 1;

    end;

end;    %loop through participant to create data cell array

%So I will keep all model probabilities between 0 and 1 to make easier use
%of sigmoid and keep everything consistant. I will divide participant
%behaviour by 100 when computing loss and multiple model probabilities by
%100 when plotting in comparison with participants to put on same scale.

v     = 6.25;
lap_out_string = 'lap_';

% %recency
% var_names = {'Pid','Suspect','Loss','Prior','Window','Bias','Noise'};
% model_struct(1).name = 'recency';
% model_struct(1).model_fitting_results = array2table(nan(0,numel(var_names)), 'VariableNames',var_names);
% model_struct(1).model_behaviour_results = []; %recency
% model_struct(1).sequences = [];
% model_struct(1).prior = struct('mean',[0 13 -10 10]','variance',v);
% model_struct(1).loss_function = @get_loss_recency;
% model_struct(1).behavior_function = @prob_guilt_recency;
% 
% %primacy
% var_names = {'Pid','Suspect','Loss','Prior','Window','Bias','Noise'};
% model_struct(2).name = 'primacy';
% model_struct(2).model_fitting_results = array2table(nan(0,numel(var_names)), 'VariableNames',var_names);
% model_struct(2).model_behaviour_results = [];  %primacy
% model_struct(2).sequences = [];
% model_struct(2).prior = struct('mean',[0 13 -10 10]','variance',v);
% model_struct(2).loss_function = @get_loss_primacy;
% model_struct(2).behavior_function = @prob_guilt_primacy;

%delta
var_names = {'Pid','Suspect','Loss','Prior','Alpha','Beta','Beta_var'};
model_struct(1).name = 'delta';
model_struct(1).model_fitting_results = array2table(nan(0,numel(var_names)), 'VariableNames',var_names);
model_struct(1).model_behaviour_results = []; %delta
model_struct(1).sequences = [];
model_struct(1).prior = struct('mean',[0 0 0 0]','variance',v); % note dimension of 'mean'
model_struct(1).loss_function = @get_loss_delta;
model_struct(1).behavior_function = @prob_guilt_delta;

% %split
% var_names = {'Pid','Suspect','Loss','Prior','Split','Bias','Noise'};
% model_struct(4).name = 'split';
% model_struct(4).model_fitting_results = array2table(nan(0,numel(var_names)), 'VariableNames',var_names);
% model_struct(4).model_behaviour_results = [];  %split
% model_struct(4).sequences = [];
% model_struct(4).prior = struct('mean',[0 .424 -10 10]','variance',v);
% model_struct(4).loss_function = @get_loss_split;
% model_struct(4).behavior_function = @prob_guilt_split;

num_models = size(model_struct,2);

for model = 1:num_models;

%     cbm_lap(P_data_A, model_struct(model).loss_function, model_struct(model).prior, [lap_out_string model_struct(model).name '_atheist.mat']);
    cbm_lap(P_data_C, model_struct(model).loss_function, model_struct(model).prior, [lap_out_string model_struct(model).name '_christian.mat']);

end;    %loop through models



%atheist model comparison
models = {@get_loss_recency,@get_loss_primacy,@get_loss_delta,@get_loss_split};
fcbm_maps = {'lap_recency_atheist.mat','lap_primacy_atheist.mat','lap_delta_atheist.mat','lap_split_atheist.mat'};
fname_hbi = 'forensic_beads_cbm_modelComparison_A.mat';
cbm_hbi(P_data_A,models,fcbm_maps,fname_hbi);

%christian model comparison
models = {@get_loss_recency,@get_loss_primacy,@get_loss_delta,@get_loss_split};
fcbm_maps = {'lap_recency_christian.mat','lap_primacy_christian.mat','lap_delta_christian.mat','lap_split_christian.mat'};
fname_hbi = 'forensic_beads_cbm_modelComparison_C.mat';
cbm_hbi(P_data_C,models,fcbm_maps,fname_hbi);


cbm_hbi_plot('forensic_beads_cbm_modelComparison.mat',{'recency','primacy','delta'},{'prior','alpha','beta'},{'sigmoid','sigmoid','exp'})
cbm_hbi_plot('forensic_beads_cbm_modelComparison_C.mat',{'recency','primacy','delta'},{'prior','alpha','beta'},{'sigmoid','sigmoid','exp'})

load forensic_beads_cbm_modelComparison_A;

winning_model_struct = model_struct(winning_model_it);

%Add the accumulated model probabilities as column of data table
data.ModelProbability = winning_model_struct.model_behaviour_results;

%Add the accumulated model probabilities as column of data table
data.GroundTruthProbability = ground_truth_behaviour_results;

%add adjustment columns to data table

%This is a distraction but I need it if I'm going to repeatedly use this
%code during debugging
columns_in_question = {'HumanAdjust','ModelAdjust','GroundTruthAdjust'};
columns_exist = ismember(columns_in_question, data.Properties.VariableNames);
if any(columns_exist)
    data = removevars(data, columns_in_question(columns_exist));
end;
data = [ ...
    data, ...
    array2table(get_adjustments(data),'VariableNames',{'HumanAdjust','ModelAdjust','GroundTruthAdjust'}) ...
    ];


suspect_0 = winning_model_struct.model_fitting_results.Prior(winning_model_struct.model_fitting_results.Suspect==0);
suspect_1 = winning_model_struct.model_fitting_results.Prior(winning_model_struct.model_fitting_results.Suspect==1);

%trad t-test
[th tpvals tci tstats] =  ttest2(suspect_0, suspect_1);

%paired Bayesian test
[bf10,bfpvals,bfci,bfstats] = bf.ttest2( suspect_0, suspect_1);

disp(sprintf('Suspect effect on prior: trad pval: %0.2f bf10: %7.4f',tpvals,bf10));

if skip_models == 1; winning_model_name = []; end;

%Plot probabilities by sequence position for different contexts and suspects
plot_sequences(data, winning_model_name);

%Plot adjustments by suspect and preceding context
plot_adjustments(data, winning_model_name);



if ~isempty(winning_model_name);

    %Some analysis of residuals
    figure;

    %normally distributed probabilities (all trials, ignoring participant)?
    subplot(2,2,1);
    hist(winning_model_struct.model_behaviour_results);
    title('Probabilities, all trials')
    xlabel('Model probability');
    ylabel('Number of trials');

    %normally distributed probabilities (participant averages)?
    subplot(2,2,2)
    hist(grpstats(winning_model_struct.model_behaviour_results,data.Pid,{'mean'}));
    title('Probabilities, participants')
    xlabel('Mean(model probability)');
    ylabel('Number of participants');

    %all residuals
    residuals_trials = winning_model_struct.model_behaviour_results - data.HumanProbability;

    subplot(2,2,3)
    hist(residuals_trials);
    title('Raw residuals, trials');
    xlabel('Model probability - human probability');
    ylabel('Number of trials');

    subplot(2,2,4)
    hist(grpstats(residuals_trials,data.Pid,{'mean'}));
    title('Raw residuals, participants');
    xlabel('Mean(model probability - human probability)');
    ylabel('Number of participants');

end;    %skip models?


disp('audi5000');

toc
%%%%%%%%%%%%%%%%%%end, forensic_beads_study2_2024_v2 (main function body)%%%%%%%%%%%%%%%%%%%%%%











%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [] = plot_sequence_panels(plot_details)

figure(plot_details.fig);

line_colors = [ones(2,3); repmat(plot_details.color,2,1)];
line_widths = [.5 .5 1 1];

y_pos = [30 90 25 80];
markersize = 4;
asterisk_y_displacement = 5;
asterisk_size = 6;

sequence_type_string = {'Innocent sequences' 'Guilty sequences'};

%There will always be one plot for innocent and one for guilty
for plot = 1:2; %innocent and guilty sequences

    subplot(plot_details.figRows,plot_details.figCols, plot_details.subplot_num+plot-1); hold on;

    if  plot_details.subplot_num == 1;
        title(sequence_type_string{plot},'Fontname','Arial','Fontsize',plot_details.graph_font,'FontWeight','normal');
    end;

    if plot == 1;
        its = [1 3];
    else
        its = [2 4];
    end;


    for i=its;    %atheist and christian lines

        %To keep lines/points clear of each other I add a tiny jitter in x
        if i == 1;
            jitter = 0.05;
        else
            jitter = -0.05;
        end;

        %Plot errors
        %use straight lines for errors
        %Since you're looping through sequence positions now
        %anyway, plot the significance asterisks whilst you're at it
        for j=1:size(plot_details.cis,2);   %sequences positions

            %error bar
            hl = line( ...
                [(j-1)+jitter (j-1)+jitter] ...
                ,[plot_details.cis(i,j,2) plot_details.cis(i,j,1)] ...
                , 'Color',plot_details.color ...
                , 'Marker','none' ...
                ,'LineStyle','-' ...
                ,'LineWidth',1 ...
                );
            hl.Annotation.LegendInformation.IconDisplayStyle = 'off';

            %asterisk, a little above the atheist (higher) error bar
            if i == 1 | i == 2;
                if ~isempty(plot_details.posthocs);
                    if plot_details.posthocs(plot,j) == 1;
                        h2 = line( ...
                            (j-1)+jitter ...
                            , plot_details.cis(i,j,2) + asterisk_y_displacement ...
                            , 'Color',[1 0 0] ...
                            , 'Marker','*' ...
                            , 'MarkerSize',asterisk_size ...
                            ,'LineStyle','none' ...
                            );
                        h2.Annotation.LegendInformation.IconDisplayStyle = 'off';
                    end;    %If posthoc is sig
                end;    %if there are posthocs to plot
            end;    %If plotting the top line (atheists)

        end; %Loop through sequence positions

        %Plot means and lines between them
        line([0:10]+jitter, plot_details.means(i,:) ...
            , 'Color',plot_details.color ...
            , 'Marker','o' ...
            ,'MarkerSize',markersize ...
            ,'MarkerFaceColor',line_colors(i,:) ...
            ,'MarkerEdgeColor',plot_details.color ...
            ,'LineStyle','-' ...
            ,'LineWidth',line_widths(i) ...
            );

    end;    %loop through atheist and christian lines
    box off

    %plot an equal prior line for visualisation's sake
    line([0 10], [50 50] ...
        , 'Color',[.75 .75 .75] ...
        ,'LineStyle','-' ...
        ,'LineWidth',1 ...
        );

    set(gca ...
        ,'Xtick',[0:10] ...
        ,'Xticklabel',{'0' '1' '2' '3' '4' '5' '6' '7' '8' '9' '10'} ...
        ,'XTickLabelRotation',0 ...
        ,'YTick',[0:20:100] ...
        ,'Fontname','Arial' ...
        ,'Fontsize',plot_details.graph_font ...
        );
    ylim([0 100]);
    xlim([-1 10.5])
    if plot_details.subplot_num == 5;
        xlabel('Sequence position');
    end;
    if plot == 1;
        ylabel("'Guilt' probability");
    end;

end;    %loop through plots

legend_strs = { ...
    sprintf('Atheist, %s', plot_details.label) ...
    sprintf('Christian, %s', plot_details.label) ...
    };
lgd = legend(legend_strs, 'Location', 'southwest');
legend boxoff;
lgd.FontSize = plot_details.graph_font;
lgd.FontName = 'Arial';
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%END, plot_sequence_panels%%%%%%%%%%%%%%%%%










%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [] = plot_sequences(data, name);

[sus_means_p_ss sus_means_ss sus_means_m_ss seq_key] = reformat_sequence_data(data);

%---------prepare sequence plots
graph_font = 12;
plot_details.fig = figure; set(gcf,'Color',[1 1 1]);
if isempty(name);
    plot_details.figRows = 2;
else;
    plot_details.figRows = 3;
end;    %skip models?
plot_details.figCols = 2;
plot_details.graph_font = graph_font;


%--------sequence plot for ground truth performance

%get participant averages (There is no variability from participant to participant because every one
%saw the same sequences so we skip the cis and post hocs)
[sus_means_m sus_grps_m] = grpstats( sus_means_m_ss, seq_key(:,[2 3]),{'mean' 'gname'} );

plot_details.subplot_num = 1;
plot_details.means = sus_means_m;
plot_details.cis = [];
plot_details.posthocs = [];
plot_details.label = 'ground truth';
plot_details.color = [0 0 1];

%make pair of plots
plot_sequence_panels(plot_details);



%------sequence plot for human performance

%get participant averages and run posthocs on them
% [sus_means_ss sus_grps_ss] = grpstats( sequences_h, seq_key(:,[1 2 3]),{'mean' 'gname'} );
% sus_key_ss = str2double(sus_grps_ss);
posthocs = get_sequence_ttests(sus_means_ss, seq_key);

%Now get means and cis over those participants for plots
[sus_means sus_cis sus_grps] = grpstats( sus_means_ss, seq_key(:,[2 3]),{'mean' 'meanci','gname'} );

%now plot them
plot_details.subplot_num = 3;
plot_details.means = sus_means;
plot_details.cis = sus_cis;
plot_details.posthocs = posthocs;
plot_details.label = 'participants';
plot_details.color = [0 0 0];

%make pair of plots
plot_sequence_panels(plot_details);

if ~isempty(name);

    %-------sequence plot for prior model performance

    %get prior model averages and run posthocs on them
    posthocs = get_sequence_ttests(sus_means_p_ss, seq_key);

    %Now get means and cis over those participants for plots
    [sus_means_p sus_cis_p sus_grps_p] = grpstats( sus_means_p_ss, seq_key(:,[2 3]),{'mean' 'meanci','gname'} );

    plot_details.subplot_num = 5;
    plot_details.means = sus_means_p;
    plot_details.cis = sus_cis_p;
    plot_details.posthocs = posthocs;
    plot_details.label = name;
    plot_details.color = [0 1 0];

    %make pair of plots
    plot_sequence_panels(plot_details);

end; %skip models?


%%%%%%%%%%%%%%%%%%END, plot_sequence%%%%%%%%%%%%%%%%%%%%%%%%%










%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function posthocs = get_sequence_ttests(sus_means_ss, sus_key_ss)

%inputs should have 104(pids)*4(suspoects and contexts) rows. sus_means_ss
%should have 11 columns (sequence positions) of mean probabilities and
%sus_key_ss should have 3 cols of pids, suspect and context codes

corrected_p_value = 0.05/22;

%Innocent sequences
M_MI = sus_means_ss(find(sus_key_ss(:,2)==0 & sus_key_ss(:,3)==0),:);
F_MI = sus_means_ss(find(sus_key_ss(:,2)==1 & sus_key_ss(:,3)==0),:);
Diff_MI = M_MI-F_MI;
[h_MI,p,ci,stats] = ttest(Diff_MI,0,'alpha',corrected_p_value);

%Guilty sequences
M_MG = sus_means_ss(find(sus_key_ss(:,2)==0 & sus_key_ss(:,3)==1),:);
F_MG = sus_means_ss(find(sus_key_ss(:,2)==1 & sus_key_ss(:,3)==1),:);
Diff_MG = M_MG-F_MG;
[h_MG,p,ci,stats] = ttest(Diff_MG,0,'alpha',corrected_p_value);

posthocs = [h_MI; h_MG];
%%%%%%%%%%%%%%END, get_sequence_ttests%%%%%%%%%%%%%%%%%%%%%%%










%%%%%%%%%%%%%BEGIN reformat_sequence_data%%%%%%%%%%%%%%%%
function [sus_means_p_ss sus_means_ss sus_means_m_ss key_ss] = reformat_sequence_data(data)

temp = grpstats(data,{'Suspect', 'Context','SequencePosition','Pid'},{'mean'});

unique_pids = unique(data.Pid);
num_pids = length(unique_pids);  % Should be 104

sus_means_p_ss = nan(num_pids * 4, 11);  % 416x11 matrix for prior
sus_means_ss = nan(num_pids * 4, 11);  % 416x11 matrix for participants
sus_means_m_ss = nan(num_pids * 4, 11);  % 416x11 matrix for ground truth (legacy: m is for model)
key_ss = nan(num_pids * 4, 3);  % 416x3 (pid, suspect, context codes) matrix

rowIndex = 1;  % Initialize row index

for i = 1:num_pids
    pid_value = unique_pids(i);

    for suspect = 0:1
        for context = 0:1

            %update key
            key_ss(rowIndex,:) = [pid_value, suspect, context];

            % Extract relevant subset of data
            subset = temp(temp.Pid == pid_value & temp.Suspect == suspect & temp.Context == context, :);

            for seqPos = 0:10

                % Calculate the mean for the current SequencePosition
                sus_means_ss(rowIndex, seqPos + 1) = mean(subset.mean_HumanProbability(subset.SequencePosition == seqPos));  % Replace with the column you want to calculate mean for

                % Calculate the mean for the current SequencePosition
                sus_means_p_ss(rowIndex, seqPos + 1) = mean(subset.mean_ModelProbability(subset.SequencePosition == seqPos));  % Replace with the column you want to calculate mean for

                % Calculate the mean for the current SequencePosition
                sus_means_m_ss(rowIndex, seqPos + 1) = mean(subset.mean_GroundTruthProbability(subset.SequencePosition == seqPos));  % Replace with the column you want to calculate mean for

            end

            rowIndex = rowIndex + 1;  % Move to the next row in the matrix
        end
    end
end
%%%%%%%%%%%%%END reformat_sequence_data%%%%%%%%%%%%%%%%














% %%%%%%%%%%%%%BEGIN get_behaviour_for_this_ps_sequences%%%%%%%%%%%%%%%%
% function sequence_probabilities = get_behaviour_for_this_ps_sequences(this_params,model_struct);
% 
% %This version is simplified compared to its counterpart in
% %foresnic_beads_study2_2024. Here, every participant has one sequence and
% %one suspect so we just proceed only to routing this one sequence into its
% %model appropriate function for estimating probabilities
% 
% %Get data for this one sequence
% this_seq_data = model_struct.this_ps_data;
% 
% %Hand just the one sequence to get_model_behaviour, together with
% %suspect-specific parameter vector and then accumulate it with the
% %other sequences for this participant to be returned by function
% if strcmp(model_struct.model_name,'Split')
% 
%     sequence_probabilities = prob_guilt_prior(this_params,this_seq_data)*100;
% 
% elseif strcmp(model_struct.model_name,'Ground truth');
% 
%     sequence_probabilities = prob_guilt_groundTruth(this_params,this_seq_data)*100;
% 
% elseif strcmp(model_struct.model_name,'Recency');
% 
%     sequence_probabilities = prob_guilt_recency(this_params,this_seq_data)*100;
% 
% elseif strcmp(model_struct.model_name,'Primacy');
% 
%     sequence_probabilities = prob_guilt_primacy(this_params,this_seq_data)*100;
% 
% elseif strcmp(model_struct.model_name,'Delta');
% 
%     sequence_probabilities = prob_guilt_delta(this_params,this_seq_data)*100;
% 
% end    %which model?
% 
% %%%%%%%%%%%%%END get_behaviour_for_this_ps_sequences%%%%%%%%%%%%%%%%
% 







% %%%%%%%%%%%%%%DISCONFIRMATORY ADJUSTMENT PLOT BEGIN%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function plot_adjustments(data, name);

%setup plot
input_struct.fig = figure('Color',[1 1 1]);
input_struct.graph_font = 12;
input_struct.series_colours = [0.75 0.75 0.75; 0.5 0.5 0.5];
input_struct.model_colours = [0 0 1; 0 1 0];
if isempty(name); input_struct.skip_models = 1; end;


%-------suspect*claim plot

input_struct.subplot = [2,1,1];

%For some reason, the claims column still retains the 3's for the priors
%(where there is no adjustment data) and so the next step creates NaN
%entries for those, which clog things up when plotting later. I'll just
%replace the 3s with NaNs within the scope of this function.
data.Claim(data.Claim == 3) = NaN;

%Get means over trials for each participant then get means and cis over participants for errorbars
%input_struct.means now contains table with all the mean and ci data
input_struct.means = grpstats( ...
    grpstats(data,{'Claim', 'Suspect','Pid'},{'mean'}), ...
    {'Claim','Suspect'},{'mean' 'meanci'});


input_struct.legend = {'Atheist suspect', 'Christian suspect', 'Ground truth', name};

b = make_grouped_bar_with_errors(input_struct);


%-------context*claim plot

input_struct.subplot = [2,1,2];

%Get means over trials for each participant then get means and cis over participants for errorbars
%input_struct.means now contains table with all the mean and ci data
% context_type = 'Context';
context_type = 'ContextCat';
input_struct.means = grpstats( ...
    grpstats(data,{'Claim', context_type,'Pid'},{'mean'}), ...
    {'Claim',context_type},{'mean' 'meanci'});

input_struct.legend = {'Preceding innocent context', 'Preceding guilty context', 'Ground truth', name};

b = make_grouped_bar_with_errors(input_struct);
% %%%%%%%%%%%%%%DISCONFIRMATORY ADJUSTMENT PLOT ENDS%%%%%%%%%%%%%%%%%%%%%%%%%%%%































%%%%%%%%%%%%%%%%%%start, get adjustments%%%%%%%%%%%%%%%%%%%%%%
function adjustment_cols = get_adjustments(data);

%on which indices is the display screen 0 (prior rating prompt so first rating)
seq_start_indices = find(data.SequencePosition==0);

%initialise (human and model adjustment cols)
adjustment_cols = nan(size(data,1),3);

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

        %ground truth adjustments
        adjustment_cols(index,3) = ...
            data.GroundTruthProbability(index) -  data.GroundTruthProbability(index-1);

    end;    %claims

end;    %sequences
%%%%%%%%%%%%%%%%%%end, get adjustments%%%%%%%%%%%%%%%%%%%%%%







%%%%%%%%%%%%%%%%%begin, make_grouped_bar_with_errors%%%%%%%%%%%%%%%%%%%%%%
function sp_h = make_grouped_bar_with_errors(input_struct);

figure(input_struct.fig);

sp_h = subplot(input_struct.subplot(1), input_struct.subplot(2), input_struct.subplot(3));

%---------Human participants

%get human participant data
human_means = reshape(input_struct.means.mean_mean_HumanAdjust,2,2)';    %after transpose, claims in rows, suspects in cols
human_cis = human_means - reshape(input_struct.means.meanci_mean_HumanAdjust(:,1),2,2)';    %after transpose, claims in rows, suspects in cols#

%get ground truth data
groundTruth_means = reshape(input_struct.means.mean_mean_GroundTruthAdjust,2,2)';    %after transpose, claims in rows, suspects in cols
groundTruth_cis = groundTruth_means - reshape(input_struct.means.meanci_mean_GroundTruthAdjust(:,1),2,2)';    %after transpose, claims in rows, suspects in cols#

%get prior model data
prior_means = reshape(input_struct.means.mean_mean_ModelAdjust,2,2)';    %after transpose, claims in rows, suspects in cols
prior_cis = prior_means - reshape(input_struct.means.meanci_mean_ModelAdjust(:,1),2,2)';    %after transpose, claims in rows, suspects in cols#

%add human participants in bars
b = bar(human_means,'grouped');
hold on;

% Customize the face colors
b(1).FaceColor = [1 1 1];  % Grey color for the first bar series
b(2).FaceColor = [1 1 1];        % White color for the second bar series

% Optional: Customize edge colors or other properties
b(1).EdgeColor = input_struct.series_colours(1,:);  % Black edge color for the first series
b(2).EdgeColor = input_struct.series_colours(2,:);  % Black edge color for the second series

% Optional: Customize edge colors or other properties
b(1).LineWidth = 2;  % Black edge color for the first series
b(2).LineWidth = 2;  % Black edge color for the second series

% Calculate the number of groups and number of bars in each group
[ngroups,nbars] = size(human_means);

% Get the x coordinate of the bars
% x = nan(nbars, ngroups);
if ngroups == 1;
    x = b.XEndPoints;
else
    for i = 1:nbars
        x(i,:) = b(i).XEndPoints;
    end
end;


model_line_width = 0.075;
x= x';
jitter = .025;
for x_pos = 1:numel(x)


    % Plot the human data error bars
    %     e1 = errorbar(x(pos),human_means(x_pos),human_cis(x_pos),'LineWidth', 2, 'Color', [0.5 0.5 0.5],'linestyle','none');
    if x_pos == 1 | x_pos == 3;
        colors = input_struct.series_colours(1,:);
    else;
        colors = input_struct.series_colours(2,:);
    end;
    line( ...
        [x(x_pos) x(x_pos)] ...
        ,[human_means(x_pos) - human_cis(x_pos) human_means(x_pos) + human_cis(x_pos)] ...
        , 'Color',colors ...
        , 'Marker','none' ...
        ,'LineStyle','-' ...
        ,'LineWidth',2 ...
        );

    %ground truth mean line
    line( ...
        [x(x_pos)-model_line_width-jitter x(x_pos)+model_line_width-jitter] ...
        ,[groundTruth_means(x_pos) groundTruth_means(x_pos)] ...
        , 'Color',input_struct.model_colours(1,:) ...
        , 'Marker','none' ...
        ,'LineStyle','-' ...
        ,'LineWidth',2 ...
        );

    %ground truth errorbar
    line( ...
        [x(x_pos)-jitter x(x_pos)-jitter] ...
        ,[groundTruth_means(x_pos) - groundTruth_cis(x_pos) groundTruth_means(x_pos) + groundTruth_cis(x_pos)] ...
        , 'Color',input_struct.model_colours(1,:) ...
        , 'Marker','none' ...
        ,'LineStyle','-' ...
        ,'LineWidth',2 ...
        );

    if input_struct.skip_models ~= 1;

        %winning model mean line
        line( ...
            [x(x_pos)-model_line_width+jitter x(x_pos)+model_line_width+jitter] ...
            ,[prior_means(x_pos) prior_means(x_pos)] ...
            , 'Color',input_struct.model_colours(2,:) ...
            , 'Marker','none' ...
            ,'LineStyle','-' ...
            ,'LineWidth',2 ...
            );

        %winning model errorbar
        line( ...
            [x(x_pos)+jitter x(x_pos)+jitter] ...
            ,[prior_means(x_pos) - prior_cis(x_pos) prior_means(x_pos) + prior_cis(x_pos)] ...
            , 'Color',input_struct.model_colours(2,:) ...
            , 'Marker','none' ...
            ,'LineStyle','-' ...
            ,'LineWidth',2 ...
            );

    end;    %skip models?

end;    %loop through claim*suspect conditions / pocitions on x axis (x_pos)

box off;
xlim([.5 2.5]);
ylim([-16 16]);
ylabel('Mean adjustment towards guilty');
set(gca ...
    ,'YTick',[-15:3:15] ...
    ,'Xtick',[1 2] ...
    ,'Xticklabel',{sprintf('Innocent claims')  sprintf('Guilty claims')} ...
    ,'Fontname','Arial' ...
    ,'Fontsize',input_struct.graph_font ...
    );
box off;

text(.75,15,input_struct.legend{1},'Color',input_struct.series_colours(1,:),'FontName','Arial','FontSize',input_struct.graph_font);
text(.75,12,input_struct.legend{2},'Color',input_struct.series_colours(2,:),'FontName','Arial','FontSize',input_struct.graph_font);

if input_struct.skip_models ~= 1
    text(.75,9,input_struct.legend{3},'Color',input_struct.model_colours(1,:),'FontName','Arial','FontSize',input_struct.graph_font);
    text(.75,6,input_struct.legend{4},'Color',input_struct.model_colours(2,:),'FontName','Arial','FontSize',input_struct.graph_font);
end;    %skip models?

fprintf('')
%%%%%%%%%%%%%%%%%end, make_grouped_bar_with_errors%%%%%%%%%%%%%%%%%%%%%%






%updated for study 1
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







%updated for study 1
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










%%%%%%%%%%%%%%%%%%START, fit_model%%%%%%%%%%%%%%%%%%%%%%
function model_struct = fit_model(params, model_struct);

% %fit params to data
% options = optimset('MaxFunEvals',10000);
% [new_params, new_loss, flag search] = ...
%     fminsearchbnd( ...
%     @(params) get_model_loss(params, model_struct), ...
%     params, ...
%     model_struct.lower_bounds, ...  %lower parameter bounds
%     model_struct.upper_bounds, ... %upper parameter bounds
%     options ...
%     );

model_struct.sequences = [ ...
    model_struct.sequences; ...
    [NaN; model_struct.this_ps_data.Claim(2:end)]' ...
    ];

%save parameter fitting output (see line before loop for col names)
model_struct.model_fitting_results = ...
    [model_struct.model_fitting_results; ...
    num2cell([model_struct.participant model_struct(1).this_ps_data.Suspect(1) new_loss new_params])
    ];

%get performance for this model
%This needs to be divided by suspect, once for each of the
%suspect-specific parameters,both of which are estimated on this loop

%Accumulate results for plotting. Need to re-loop suspects if Study 2
model_struct.model_behaviour_results = [ ...
    model_struct.model_behaviour_results;
    get_behaviour_for_this_ps_sequences(new_params,model_struct) ...
    ];
%%%%%%%%%%%%%%%%%%END, fit_model%%%%%%%%%%%%%%%%%%%%%%









% %%%%%%%%%%%%%%%%%%START, prob_guilt_groundTruth%%%%%%%%%%%%%%%%%%%%%%
function model_probabilities = prob_guilt_groundTruth(params, this_ps_suspect_data)

%Was previously get_model_behaviour. This prior model (which can't explain
%disconfirmatory bias) is deprecated in forensic_beads_study2_2024_v3.

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

        %conditional probability
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
%%%%%%%%%%%%%%%%%%START, prob_guilt_groundTruth%%%%%%%%%%%%%%%%%%%%%%


















%%%%%%%%%%%%%%%%%%start, prob_guilt_recency%%%%%%%%%%%%%%%%%%%%%%
function model_probabilities = prob_guilt_recency(params, this_ps_suspect_data)

%Was previously get_model_behaviour. This prior model (which can't explain
%disconfirmatory bias) is deprecated in forensic_beads_study2_2024_v3.

% if params(1) < 0; prior = 0; elseif params(1) > 1; prior = 1; else prior = params(1); end;
% if params(2) < 1; window = 1; elseif params(2) > 10; window = 10; else window = round(params(2)); end;
% if params(3) < 1; response_bias = 1; elseif params(3) > 10; response_bias = 10; else response_bias = params(3); end;
% if params(4) < 0; response_noise = 0; elseif params(4) > 1; response_noise = 1; else response_noise = params(4); end;

prior = 1./(1+exp(-params(1)));   %bound between zero and one
window = round(sigmoid_oneToten(params(2)));
response_bias = 1./(1+exp(-params(3)));
response_noise = 1./(1+exp(-params(4)));

%on which indices is the display screen 0 (prior rating prompt so first rating)
seq_start_indices = find(this_ps_suspect_data.SequencePosition==0);

%initialise output
model_probabilities = nan(size(this_ps_suspect_data,1),1);

%For each start index, loop through sequence and get model predictions
for seq = 1:numel(seq_start_indices);

    draws = this_ps_suspect_data.Claim(seq_start_indices(seq):seq_start_indices(seq)+10);

    %Loop through this sequence
    for claim = 1:11;

        %change the actual sequence up to t
        draws_recency = draws(max(1, (claim+1) - window):claim);

        %leading nan doesn't count
        ng = nansum(draws_recency);
        nd = sum(~isnan(draws_recency));

        q = .6;

        %conditional probability
        noiseless_p = (1/(1 + ((1-prior)/prior)*(q/(1-q))^(nd-2*ng)));

        %add noise and response bias
        noise_p =        response_bias + response_noise*noiseless_p;

        if noise_p <= 0;
            noise_p = 0;
        elseif noise_p >= 1;
            noise_p = 1;
        end;

        model_probabilities(claim,1) = noise_p;

    end;    %loop through this sequence (claim)

end;    %loop through sequences
%%%%%%%%%%%%%%%%%%end, guilt_prob_recency%%%%%%%%%%%%%%%%%%%%%%








%%%%%%%%%%%%%%%%%%start, prob_guilt_primacy%%%%%%%%%%%%%%%%%%%%%%
function model_probabilities = prob_guilt_primacy(params, this_ps_suspect_data)

%Was previously get_model_behaviour. This prior model (which can't explain
%disconfirmatory bias) is deprecated in forensic_beads_study2_2024_v3.

% if params(1) < 0; prior = 0; elseif params(1) > 1; prior = 1; else prior = params(1); end;
% if params(2) < 1; window = 1; elseif params(2) > 10; window = 10; else window = round(params(2)); end;
% if params(3) < 1; response_bias = 1; elseif params(3) > 10; response_bias = 10; else response_bias = params(3); end;
% if params(4) < 0; response_noise = 0; elseif params(4) > 1; response_noise = 1; else response_noise = params(4); end;

prior = 1./(1+exp(-params(1)));   %bound between zero and one
window = round(sigmoid_oneToten(params(2)));
response_bias = 1./(1+exp(-params(3)));
response_noise = 1./(1+exp(-params(4)));


%on which indices is the display screen 0 (prior rating prompt so first rating)
seq_start_indices = find(this_ps_suspect_data.SequencePosition==0);

%initialise output
model_probabilities = nan(size(this_ps_suspect_data,1),1);

%For each start index, loop through sequence and get model predictions
for seq = 1:numel(seq_start_indices);

    draws = this_ps_suspect_data.Claim(seq_start_indices(seq):seq_start_indices(seq)+10);

    %Loop through this sequence
    for claim = 1:11;

        %change the actual sequence so that samples between N and t are dropped out
        if claim > (window + 1);   %+1 because we want to ignore the leading nan
            draws_primacy = [draws(1:1+window); draws(claim)];   %The second sequence item is the first sample after the nan prior, then N past that prior nan
        else
            draws_primacy = draws(1:claim); % If inside window, just count everything
        end
        % You can now use draws_primacy for your computations

        %leading nan doesn't count
        ng = nansum(draws_primacy);
        nd = sum(~isnan(draws_primacy));

        q = .6;

        %conditional probability
        noiseless_p = (1/(1 + ((1-prior)/prior)*(q/(1-q))^(nd-2*ng)));

        %add noise and response bias
        noise_p =        response_bias + response_noise*noiseless_p;

        if noise_p <= 0;
            noise_p = 0;
        elseif noise_p >= 1;
            noise_p = 1;
        end;

        model_probabilities(claim,1) = noise_p;

    end;    %loop through this sequence (claim)

end;    %loop through sequences
%%%%%%%%%%%%%%%%%%end, guilt_prob_primacy%%%%%%%%%%%%%%%%%%%%%%





%%%%%%%%%%%%%%%%%%start, prob_guilt_delta%%%%%%%%%%%%%%%%%%%%%%
function model_probabilities = prob_guilt_delta(params, this_ps_suspect_data)


prior = 1./(1+exp(-params(1)));   %bound between zero and one
alpha = 1./(1+exp(-params(2)));
beta = exp(params(3));  %bound to be greater than 1

% if params(1) < 0; prior = 0; elseif params(1) > 1; prior = 1; else prior = params(1); end;
% if params(2) < 0; alpha = 0; elseif params(2) > 1; alpha = 1; else alpha = params(2); end;
% if params(3) < 0; beta = 0; else beta = params(3); end;

%on which indices is the display screen 0 (prior rating prompt so first rating)
seq_start_indices = find(this_ps_suspect_data.SequencePosition==0);

%initialise output
model_probabilities = nan(size(this_ps_suspect_data,1),1);

%For each start index, loop through sequence and get model predictions
for seq = 1:numel(seq_start_indices);

    draws = this_ps_suspect_data.Claim(seq_start_indices(seq):seq_start_indices(seq)+10);

    %Loop through this sequence
    for claim = 1:11;

        if claim == 1;

            q_hat = prior;  %first update, without info, based just on prior

        else;

            q_hat = q_hat + alpha * (draws(claim) - q_hat);

        end;

        model_probabilities(claim,1) = exp(beta*q_hat)/(exp(beta*q_hat) + exp(beta*(1-q_hat)));

    end;    %loop through this sequence (claim)

end;    %loop through sequences

%%%%%%%%%%%%%%%%%%end, guilt_prob_delta%%%%%%%%%%%%%%%%%%%%%%




%%%%%%%%%%%%%%%%%%start, prob_guilt_split%%%%%%%%%%%%%%%%%%%%%%
function model_probabilities = prob_guilt_split(params, this_ps_suspect_data)

%Was previously get_model_behaviour. This prior model (which can't explain
%disconfirmatory bias) is deprecated in forensic_beads_study2_2024_v3.
% 
% if params(1) < 0; prior = 0; elseif params(1) > 1; prior = 1; else prior = params(1); end;
% if params(2) < 0; split = 0; elseif params(2) > 1; split = 1; else split = params(2); end;
% if params(3) < 1; response_bias = 1; elseif params(3) > 10; response_bias = 10; else response_bias = params(3); end;
% if params(4) < 0; response_noise = 0; elseif params(4) > 1; response_noise = 1; else response_noise = params(4); end;

prior = 1./(1+exp(-params(1)));   %bound between zero and one
split = 1./(1+exp(-params(2)));   %bound between zero and one
response_bias = 1./(1+exp(-params(3)));
response_noise = 1./(1+exp(-params(4)));

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

        %conditional probability
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
%%%%%%%%%%%%%%%%%%end, guilt_prob_split%%%%%%%%%%%%%%%%%%%%%%







%%%%%%%%%%%%%%%%%%start, get_loss_recency%%%%%%%%%%%%%%%%%%%%%%
function loss = get_loss_recency(params, this_ps_suspect_data)

fitted_probabilities = prob_guilt_recency(params, this_ps_suspect_data);

loss = get_model_loss(fitted_probabilities,this_ps_suspect_data);

%%%%%%%%%%%%%%%%%%end, get_loss_delta%%%%%%%%%%%%%%%%%%%%%%




%%%%%%%%%%%%%%%%%%start, get_loss_primacy%%%%%%%%%%%%%%%%%%%%%%
function loss = get_loss_primacy(params, this_ps_suspect_data)

fitted_probabilities = prob_guilt_primacy(params, this_ps_suspect_data);

loss = get_model_loss(fitted_probabilities,this_ps_suspect_data);

%%%%%%%%%%%%%%%%%%end, get_loss_primacy%%%%%%%%%%%%%%%%%%%%%%




%%%%%%%%%%%%%%%%%%start, get_loss_delta%%%%%%%%%%%%%%%%%%%%%%
function loss = get_loss_delta(params, this_ps_suspect_data)

beta_var = exp(params(end));

fitted_probabilities = prob_guilt_delta(params, this_ps_suspect_data);

loss = get_model_loss(fitted_probabilities,this_ps_suspect_data, beta_var);

%%%%%%%%%%%%%%%%%%end, get_loss_delta%%%%%%%%%%%%%%%%%%%%%%





%%%%%%%%%%%%%%%%%%start, get_loss_split%%%%%%%%%%%%%%%%%%%%%%
function loss = get_loss_split(params, this_ps_suspect_data)

fitted_probabilities = prob_guilt_split(params, this_ps_suspect_data);

loss = get_model_loss(fitted_probabilities,this_ps_suspect_data);

%%%%%%%%%%%%%%%%%%end, get_loss_delta%%%%%%%%%%%%%%%%%%%%%%






%%%%%%%%%%%%%%%%%%start, get_model_loss%%%%%%%%%%%%%%%%%%%%%%
function loss = get_model_loss(fitted_probabilities,this_ps_data,beta_var);
% 
% fitted_probabilities = ...
%     get_behaviour_for_this_ps_sequences(params,model_struct);

% %sum squared error loss
% loss = 0;
% for trial = 1:size(this_ps_data,1);
% 
%     %squared error loss, accumulating over both suspects
%     y_hat = fitted_probabilities(trial,end);
%     y = this_ps_data.HumanProbability(trial,1)/100;  %human / participant probability is col 4.
%     loss = loss + (y_hat - y)^2;
% 
% end;    %loop through trials, this suspect

phi = 1/beta_var;
%negative log-likelihood loss, based on beta distribution
%var_beta = 1/phi;
ll = 0;
for trial = 1:size(this_ps_data,1);

    nu = (fitted_probabilities(trial)*(1-fitted_probabilities(trial)))/beta_var;
    a = fitted_probabilities(trial)*nu;
    b = (1-fitted_probabilities(trial))*nu;

    log_pdf = log(betapdf(this_ps_data.HumanProbability(trial,1)/100,a,b));

    if isnan(log_pdf); %can happen if beta gets too huge inside of exp
        log_pdf = Inf;
    end;

    ll = ll - log_pdf;

end;    %trials / draws / samples / whatever

% ll_phi = 0;
% for trial = 1:size(this_ps_data,1);
% 
%     a = fitted_probabilities(trial)*phi;
%     b = (1-fitted_probabilities(trial))*phi;
% 
%     log_pdf_phi = log(betapdf(this_ps_data.HumanProbability(trial,1)/100,a,b))
% 
%     ll_phi = ll_phi - log_pdf_phi;
% 
% end;    %trials / draws / samples / whatever
% 
% ll_phi_g = 0;
% for trial = 1:size(this_ps_data,1);
% 
%     a = fitted_probabilities(trial)*phi;
%     b = (1-fitted_probabilities(trial))*phi;
% 
%     log_pdf_phi_g = log(betapdf(this_ps_data.HumanProbability(trial,1)/100,a,b))
% 
%     log_pdf = log(gamma(a + b)) - log(gamma(a)) - log(gamma(b)) ...
%               + (a - 1) * log(this_ps_data.HumanProbability(trial,1)/100) + (b - 1) * log(1 - this_ps_data.HumanProbability(trial,1)/100);
% 
%     ll_phi_g = ll_phi_g - log_pdf_phi_g;
% 
% end;    %trials / draws / samples / whatever

loss = ll;

%%%%%%%%%%%%%%%%%%end, get_model_loss%%%%%%%%%%%%%%%%%%%%%%



%%%%%%%%%%%%%%%%%%%%%%start, sigmoid_oneToten%%%%%%%%%%%%%%%%%
function out = sigmoid_oneToten(in);

out = (9./(1+exp(-in)))+1;
%%%%%%%%%%%%%%%%%%%%%%end, sigmoid_oneToten%%%%%%%%%%%%%%%%%


