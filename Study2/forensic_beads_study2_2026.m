function [sus_means_p_ss sus_key_p_ss] = forensic_beads_study2_2026;

% Reduced black-and-white Study 2 version:
% - Fits no participant-level/free-parameter models.
% - Computes only the objective ground-truth / ideal-observer probabilities.
% - Sequence-position figure averages over suspect gender, so there are no
%   separate male/female lines.
% - Sequence-position figure overlays participants and ground truth on the
%   same plots, with one panel for mostly innocent and one for mostly guilty
%   sequences.
% - Gender figure has two subplots:
%       left  = claim * suspect-gender adjustment plot
%       right = sequence-position-0 / prior probability by suspect gender
% - Claim * context adjustment plot is a separate figure.
% - All plots are black/white/grayscale only.
% - Keeps dummy ModelProbability/ModelAdjust columns as NaNs so legacy table
%   helper structure remains stable.
%
% NOTE:
% The returned output name sus_means_p_ss is retained for compatibility with
% the old function signature. In this reduced version, it is assigned to the
% participant sequence means, not fitted-model means.

tic

% Model-fitting dependencies removed.
% This path is probably unnecessary unless you use functions from it elsewhere.
addpath(genpath('C:\matlab_files\fiance\forensic_beads_pub_repo\Forensic-beads-paper-1\shaded_plots'))

%%%%%%%%
%1: event index,
%2: participant private id,
%3: RT,
%4: human probability
%5: sequence position
%6: suspect gender (1=female, 0=male)
%7: witness gender (1=female)
%8: guilty claim (1=guilty)
%9: context (mostly innocent / mostly guilty)
%10 preceding context degree
%11 preceding context category

var_names = ...
    {'Event', ...
    'Pid', ...
    'RT', ...
    'HumanProbability', ...
    'SequencePosition', ...
    'Suspect', ...
    'WitnessGender', ...
    'Claim', ...
    'Context', ...
    'ContextDegree', ...
    'ContextCat'};

data = ...
    array2table( ...
    get_sub_data, ...
    'VariableNames', var_names ...
    );

%Who are the participants?
participant_list = unique(data.Pid, 'stable');
num_participants = numel(participant_list);

%ground truth only
ground_truth_behaviour_results = [];

%compute ground truth participant by participant
for participant = 1:num_participants

    clear this_ps_data

    %get data for this participant
    this_ps_data = ...
        data( ...
        data.Pid == participant_list(participant), ...
        :);

    disp(sprintf('computing ground truth for participant %d', participant))

    %-----------Ground truth model-----------------
    % Parameter vector:
    % [prior_male prior_female split_male split_female response_bias response_noise]
    %
    % Ground truth:
    % prior = .5 for both suspect genders
    % split = .7 for both suspect genders
    % response bias = 0
    % response noise = 1
    ground_truth_params = [0.5 0.5 0.7 0.7 0 1];

    ground_truth_struct.this_ps_data = this_ps_data;
    ground_truth_struct.model_name = 'Ground truth';

    ground_truth_behaviour_results = [ ...
        ground_truth_behaviour_results; ...
        get_behaviour_for_this_ps_sequences(ground_truth_params, ground_truth_struct) ...
        ];

end

%Add ground-truth probabilities as column of data table
data.GroundTruthProbability = ground_truth_behaviour_results;

%Keep a dummy model column so existing helper functions remain stable.
%No fitted models are computed or plotted.
data.ModelProbability = nan(height(data), 1);

%Add adjustment columns to data table
columns_in_question = {'HumanAdjust', 'ModelAdjust', 'GroundTruthAdjust'};
columns_exist = ismember(columns_in_question, data.Properties.VariableNames);

if any(columns_exist)
    data = removevars(data, columns_in_question(columns_exist));
end

data = [ ...
    data, ...
    array2table( ...
    get_adjustments(data), ...
    'VariableNames', {'HumanAdjust', 'ModelAdjust', 'GroundTruthAdjust'} ...
    ) ...
    ];

%Empty model name tells plotting functions to skip fitted-model rows/lines.
winning_model_name = [];

%Plot probabilities by sequence position, averaged over suspect gender
plot_sequences(data, winning_model_name);

%Plot participant male/female sequence-position difference figure
plot_gender_sequence_difference(data);

%Plot gender figure and context figure
plot_adjustments(data, winning_model_name);

%Return participant sequence means and key for compatibility with old output.
%The first output name is legacy.
[~, sus_means_p_ss, ~, sus_key_p_ss] = reformat_sequence_data(data);

disp('audi5000');

toc

end
%%%%%%%%%%%%%%%%%%end main function%%%%%%%%%%%%%%%%%%%%%%











%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [] = plot_sequence_panels_average(plot_details)

% This plots participants and ground truth on the same axes.
% Data have already been averaged over suspect gender.
% There are two panels:
%   left  = mostly innocent sequences
%   right = mostly guilty sequences

figure(plot_details.fig);

participant_marker_size = 4;
ground_truth_marker_size = 4;

sequence_type_string = {'Mostly innocent sequences' 'Mostly guilty sequences'};

for plot = 1:2

    subplot(1, 2, plot);
    hold on;

    title(sequence_type_string{plot}, ...
        'Fontname', 'Arial', ...
        'Fontsize', plot_details.graph_font, ...
        'FontWeight', 'normal');

    % Row index:
    % 1 = mostly innocent context
    % 2 = mostly guilty context
    row_to_plot = plot;

    % Participant CIs
    if ~isempty(plot_details.human_cis)

        for j = 1:size(plot_details.human_cis, 2)

            hl = line( ...
                [j-1 j-1], ...
                [plot_details.human_cis(row_to_plot,j,2) plot_details.human_cis(row_to_plot,j,1)], ...
                'Color', [0 0 0], ...
                'Marker', 'none', ...
                'LineStyle', '-', ...
                'LineWidth', 1);
            hl.Annotation.LegendInformation.IconDisplayStyle = 'off';

        end

    end

    % Participant means: solid line, open circles
    h_participants = line([0:10], plot_details.human_means(row_to_plot,:), ...
        'Color', [0 0 0], ...
        'Marker', 'o', ...
        'MarkerSize', participant_marker_size, ...
        'MarkerFaceColor', [1 1 1], ...
        'MarkerEdgeColor', [0 0 0], ...
        'LineStyle', '-', ...
        'LineWidth', 1.5);

    % Ground truth means: dashed line, filled squares
    h_ground_truth = line([0:10], plot_details.ground_truth_means(row_to_plot,:), ...
        'Color', [0 0 0], ...
        'Marker', 's', ...
        'MarkerSize', ground_truth_marker_size, ...
        'MarkerFaceColor', [0 0 0], ...
        'MarkerEdgeColor', [0 0 0], ...
        'LineStyle', '--', ...
        'LineWidth', 1.5);

    box off

    % Equal-prior reference line
    h_prior = line([0 10], [50 50], ...
        'Color', [0 0 0], ...
        'LineStyle', ':', ...
        'LineWidth', 1);
    h_prior.Annotation.LegendInformation.IconDisplayStyle = 'off';

    set(gca, ...
        'Xtick', [0:10], ...
        'Xticklabel', {'0' '1' '2' '3' '4' '5' '6' '7' '8' '9' '10'}, ...
        'XTickLabelRotation', 0, ...
        'YTick', [0:20:100], ...
        'Fontname', 'Arial', ...
        'Fontsize', plot_details.graph_font);

    ylim([0 100]);
    xlim([-1 10.5])

    xlabel('Sequence position');

    if plot == 1
        ylabel("'Guilt' probability");
    end

    lgd = legend( ...
        [h_participants h_ground_truth], ...
        {'Participants', 'Ground truth'}, ...
        'Location', 'southwest');
    legend boxoff;
    lgd.FontName = 'Arial';
    lgd.FontSize = plot_details.graph_font;

end

end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%END, plot_sequence_panels_average%%%%%%%%%%%%%%%%%










%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [] = plot_sequences(data, name)

% This version averages over suspect gender before plotting.
% Participants and ground truth are plotted on the same axes.
% The figure has two columns:
%   left  = mostly innocent sequences
%   right = mostly guilty sequences

[human_means_ss, ground_truth_means_ss, seq_key] = reformat_sequence_data_average_over_suspect(data);

%---------prepare sequence plots
graph_font = 12;
plot_details.fig = figure;
set(gcf, 'Color', [1 1 1]);
plot_details.graph_font = graph_font;

%--------ground truth means

% Get ground-truth means over participants, grouped by context only.
% There is no true participant variability here because everyone saw the
% same sequences, so no ground-truth CIs are plotted.
[ground_truth_means, ground_truth_grps] = grpstats( ...
    ground_truth_means_ss, ...
    seq_key(:,2), ...
    {'mean' 'gname'} ...
    );

%--------human participant means and CIs

[human_means, human_cis, human_grps] = grpstats( ...
    human_means_ss, ...
    seq_key(:,2), ...
    {'mean' 'meanci' 'gname'} ...
    );

plot_details.human_means = human_means;
plot_details.human_cis = human_cis;
plot_details.ground_truth_means = ground_truth_means;

plot_sequence_panels_average(plot_details);

end
%%%%%%%%%%%%%%%%%%END, plot_sequences%%%%%%%%%%%%%%%%%%%%%%%%%










%%%%%%%%%%%%%BEGIN reformat_sequence_data_average_over_suspect%%%%%%%%%%%%%%%%
function [human_means_ss, ground_truth_means_ss, key_ss] = reformat_sequence_data_average_over_suspect(data)

% This function averages over suspect gender within each participant,
% context, and sequence position. It is used only for the sequence-position
% figure.
%
% Output rows are pid * context:
% key_ss columns:
%   1 = pid
%   2 = context

temp = grpstats(data, {'Context', 'SequencePosition', 'Pid'}, {'mean'});

unique_pids = unique(data.Pid);
num_pids = length(unique_pids);

human_means_ss = nan(num_pids * 2, 11);
ground_truth_means_ss = nan(num_pids * 2, 11);
key_ss = nan(num_pids * 2, 2);

rowIndex = 1;

for i = 1:num_pids

    pid_value = unique_pids(i);

    for context = 0:1

        key_ss(rowIndex,:) = [pid_value, context];

        subset = temp( ...
            temp.Pid == pid_value & ...
            temp.Context == context, ...
            :);

        for seqPos = 0:10

            human_means_ss(rowIndex, seqPos + 1) = ...
                mean(subset.mean_HumanProbability(subset.SequencePosition == seqPos));

            ground_truth_means_ss(rowIndex, seqPos + 1) = ...
                mean(subset.mean_GroundTruthProbability(subset.SequencePosition == seqPos));

        end

        rowIndex = rowIndex + 1;

    end

end

end
%%%%%%%%%%%%%END reformat_sequence_data_average_over_suspect%%%%%%%%%%%%%%%%










%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function posthocs = get_sequence_ttests(sus_means_ss, sus_key_ss)

% Retained for compatibility with older versions, but not used by the
% averaged-over-suspect sequence plot.

corrected_p_value = 0.05/22;

%Innocent sequences
M_MI = sus_means_ss(find(sus_key_ss(:,2)==0 & sus_key_ss(:,3)==0),:);
F_MI = sus_means_ss(find(sus_key_ss(:,2)==1 & sus_key_ss(:,3)==0),:);
Diff_MI = M_MI - F_MI;
[h_MI, p, ci, stats] = ttest(Diff_MI, 0, 'alpha', corrected_p_value);

%Guilty sequences
M_MG = sus_means_ss(find(sus_key_ss(:,2)==0 & sus_key_ss(:,3)==1),:);
F_MG = sus_means_ss(find(sus_key_ss(:,2)==1 & sus_key_ss(:,3)==1),:);
Diff_MG = M_MG - F_MG;
[h_MG, p, ci, stats] = ttest(Diff_MG, 0, 'alpha', corrected_p_value);

posthocs = [h_MI; h_MG];

end
%%%%%%%%%%%%%%END, get_sequence_ttests%%%%%%%%%%%%%%%%%%%%%%%










%%%%%%%%%%%%%BEGIN reformat_sequence_data%%%%%%%%%%%%%%%%
function [sus_means_p_ss, sus_means_ss, sus_means_m_ss, key_ss] = reformat_sequence_data(data)

% Retained for output compatibility and possible downstream analyses.
% This version still preserves suspect-gender * context rows.

temp = grpstats(data, {'Suspect', 'Context', 'SequencePosition', 'Pid'}, {'mean'});

unique_pids = unique(data.Pid);
num_pids = length(unique_pids);

sus_means_p_ss = nan(num_pids * 4, 11);  %legacy model matrix; now dummy NaNs
sus_means_ss = nan(num_pids * 4, 11);    %participants
sus_means_m_ss = nan(num_pids * 4, 11);  %ground truth
key_ss = nan(num_pids * 4, 3);           %pid, suspect, context

rowIndex = 1;

for i = 1:num_pids

    pid_value = unique_pids(i);

    for suspect = 0:1

        for context = 0:1

            %update key
            key_ss(rowIndex,:) = [pid_value, suspect, context];

            % Extract relevant subset of data
            subset = temp( ...
                temp.Pid == pid_value & ...
                temp.Suspect == suspect & ...
                temp.Context == context, ...
                :);

            for seqPos = 0:10

                % Human participants
                sus_means_ss(rowIndex, seqPos + 1) = ...
                    mean(subset.mean_HumanProbability(subset.SequencePosition == seqPos));

                % Dummy fitted model values; retained for compatibility
                if ismember('mean_ModelProbability', subset.Properties.VariableNames)
                    sus_means_p_ss(rowIndex, seqPos + 1) = ...
                        mean(subset.mean_ModelProbability(subset.SequencePosition == seqPos));
                else
                    sus_means_p_ss(rowIndex, seqPos + 1) = NaN;
                end

                % Ground truth
                sus_means_m_ss(rowIndex, seqPos + 1) = ...
                    mean(subset.mean_GroundTruthProbability(subset.SequencePosition == seqPos));

            end

            rowIndex = rowIndex + 1;

        end

    end

end

end
%%%%%%%%%%%%%END reformat_sequence_data%%%%%%%%%%%%%%%%










%%%%%%%%%%%%%BEGIN get_behaviour_for_this_ps_sequences%%%%%%%%%%%%%%%%
function sequence_probabilities = get_behaviour_for_this_ps_sequences(ground_truth_params, model_struct)

% Computes objective ground-truth behaviour only.
%
% ground_truth_params:
% [prior_male prior_female split_male split_female response_bias response_noise]

this_ps_data = model_struct.this_ps_data;

seq_start_indices = find(this_ps_data.SequencePosition == 0);

sequence_probabilities = [];

for seq = 1:numel(seq_start_indices)

    clear this_suspect_code this_params this_seq_data

    % suspect code: 0 = male, 1 = female
    this_suspect_code = this_ps_data.Suspect(seq_start_indices(seq));

    if strcmp(model_struct.model_name, 'Ground truth')

        %Pick suspect-specific prior and split.
        %For suspect 0: prior param 1, split param 3.
        %For suspect 1: prior param 2, split param 4.
        this_params = [ ...
            ground_truth_params(this_suspect_code + 1), ...
            ground_truth_params(this_suspect_code + 3), ...
            ground_truth_params(5), ...
            ground_truth_params(6) ...
            ];

        this_seq_data = this_ps_data(seq_start_indices(seq):seq_start_indices(seq)+10, :);

        sequence_probabilities = [ ...
            sequence_probabilities; ...
            prob_guilt_groundTruth(this_params, this_seq_data) * 100 ...
            ];

    else

        error('Unknown model_name. This reduced script only supports Ground truth.');

    end

end

end
%%%%%%%%%%%%%END get_behaviour_for_this_ps_sequences%%%%%%%%%%%%%%%%








%%%%%%%%%%%%%%DISCONFIRMATORY ADJUSTMENT PLOT BEGIN%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function plot_adjustments(data, name)

% In this version:
% - Gender figure has two subplots:
%     left  = adjustment by suspect gender
%     right = sequence-position-0 probability by suspect gender
% - Context adjustment plot is made as a separate figure.
% - Plots are black/white/grayscale only.

input_struct.graph_font = 12;
input_struct.series_colours = [0 0 0; 0 0 0];
input_struct.model_colours = [0 0 0; 0 0 0];
input_struct.skip_models = 0;

if isempty(name)
    input_struct.skip_models = 1;
end

%For some reason, the claims column may retain the 3's for the priors
%(where there is no adjustment data). Replace with NaNs within this function.
data.Claim(data.Claim == 3) = NaN;

%-------suspect-gender figure: adjustment plot + sequence-position-0 plot

gender_fig = figure('Color', [1 1 1]);

%Left subplot: suspect-gender * claim adjustment plot
input_struct.fig = gender_fig;
input_struct.subplot = [1, 2, 1];
input_struct.figure_title = 'Adjustment by suspect gender';
input_struct.group_labels = {'Male suspect', 'Female suspect'};
input_struct.ground_truth_label = 'Ground truth';

input_struct.means = grpstats( ...
    grpstats(data, {'Claim', 'Suspect', 'Pid'}, {'mean'}), ...
    {'Claim', 'Suspect'}, ...
    {'mean' 'meanci'} ...
    );

make_grouped_bar_with_errors(input_struct);

%Right subplot: sequence-position-0 probability by suspect gender
prior_struct.fig = gender_fig;
prior_struct.subplot = [1, 2, 2];
prior_struct.figure_title = 'Initial probability by suspect gender';
prior_struct.graph_font = input_struct.graph_font;
prior_struct.group_labels = {'Male suspect', 'Female suspect'};

make_prior_bar_with_errors(data, prior_struct);

%-------context*claim plot in its own figure

context_type = 'ContextCat';

input_struct.fig = figure('Color', [1 1 1]);
input_struct.subplot = [];
input_struct.figure_title = 'Adjustment by preceding context';
input_struct.group_labels = {'Preceding innocent context', 'Preceding guilty context'};
input_struct.ground_truth_label = 'Ground truth';

input_struct.means = grpstats( ...
    grpstats(data, {'Claim', context_type, 'Pid'}, {'mean'}), ...
    {'Claim', context_type}, ...
    {'mean' 'meanci'} ...
    );

make_grouped_bar_with_errors(input_struct);

end
%%%%%%%%%%%%%%DISCONFIRMATORY ADJUSTMENT PLOT ENDS%%%%%%%%%%%%%%%%%%%%%%%%%%%%











%%%%%%%%%%%%%%%%%%start, get adjustments%%%%%%%%%%%%%%%%%%%%%%
function adjustment_cols = get_adjustments(data)

%Indices where display screen/sequence position is 0: prior rating prompt.
seq_start_indices = find(data.SequencePosition == 0);

%Columns:
%1 human adjustments
%2 fitted model adjustments; dummy NaNs in this reduced script
%3 ground truth adjustments
adjustment_cols = nan(size(data,1), 3);

%For each sequence, loop through claims
for seq = 1:numel(seq_start_indices)

    for claim = 2:11

        index = seq_start_indices(seq) + claim - 1;

        %human adjustments
        adjustment_cols(index,1) = ...
            data.HumanProbability(index) - data.HumanProbability(index-1);

        %dummy fitted-model adjustments
        adjustment_cols(index,2) = ...
            data.ModelProbability(index) - data.ModelProbability(index-1);

        %ground truth adjustments
        adjustment_cols(index,3) = ...
            data.GroundTruthProbability(index) - data.GroundTruthProbability(index-1);

    end

end

end
%%%%%%%%%%%%%%%%%%end, get adjustments%%%%%%%%%%%%%%%%%%%%%%







%%%%%%%%%%%%%%%%%begin, make_grouped_bar_with_errors%%%%%%%%%%%%%%%%%%%%%%
function sp_h = make_grouped_bar_with_errors(input_struct)

figure(input_struct.fig);

if isfield(input_struct, 'subplot') && ~isempty(input_struct.subplot)
    sp_h = subplot(input_struct.subplot(1), input_struct.subplot(2), input_struct.subplot(3));
else
    sp_h = axes;
end

hold on;

%---------Human participants

%get human participant data
human_means = reshape(input_struct.means.mean_mean_HumanAdjust, 2, 2)';
human_cis = human_means - reshape(input_struct.means.meanci_mean_HumanAdjust(:,1), 2, 2)';

%get ground truth data
groundTruth_means = reshape(input_struct.means.mean_mean_GroundTruthAdjust, 2, 2)';
groundTruth_cis = groundTruth_means - reshape(input_struct.means.meanci_mean_GroundTruthAdjust(:,1), 2, 2)';

%get fitted-model data only if used
if input_struct.skip_models ~= 1
    prior_means = reshape(input_struct.means.mean_mean_ModelAdjust, 2, 2)';
    prior_cis = prior_means - reshape(input_struct.means.meanci_mean_ModelAdjust(:,1), 2, 2)';
end

%add human participants in bars
b = bar(human_means, 'grouped');
hold on;

%Black/white/grayscale bar styling.
%First grouping-variable level: white bar.
%Second grouping-variable level: light gray bar.
b(1).FaceColor = [1 1 1];
b(2).FaceColor = [.75 .75 .75];

b(1).EdgeColor = [0 0 0];
b(2).EdgeColor = [0 0 0];

b(1).LineWidth = 1.5;
b(2).LineWidth = 1.5;

%Get x coordinates of bars
[ngroups, nbars] = size(human_means);

if ngroups == 1
    x = b.XEndPoints;
else
    x = nan(nbars, ngroups);
    for i = 1:nbars
        x(i,:) = b(i).XEndPoints;
    end
end

model_line_width = 0.075;
x = x';
jitter = .025;

for x_pos = 1:numel(x)

    %Human data error bars: black
    line( ...
        [x(x_pos) x(x_pos)], ...
        [human_means(x_pos) - human_cis(x_pos), human_means(x_pos) + human_cis(x_pos)], ...
        'Color', [0 0 0], ...
        'Marker', 'none', ...
        'LineStyle', '-', ...
        'LineWidth', 1.5);

    %Ground truth mean line: black horizontal line, offset slightly left
    line( ...
        [x(x_pos)-model_line_width-jitter, x(x_pos)+model_line_width-jitter], ...
        [groundTruth_means(x_pos), groundTruth_means(x_pos)], ...
        'Color', [0 0 0], ...
        'Marker', 'none', ...
        'LineStyle', '-', ...
        'LineWidth', 2.5);

    %Ground truth error bar: black
    line( ...
        [x(x_pos)-jitter, x(x_pos)-jitter], ...
        [groundTruth_means(x_pos) - groundTruth_cis(x_pos), groundTruth_means(x_pos) + groundTruth_cis(x_pos)], ...
        'Color', [0 0 0], ...
        'Marker', 'none', ...
        'LineStyle', '-', ...
        'LineWidth', 1.5);

    if input_struct.skip_models ~= 1

        %Winning model mean line
        line( ...
            [x(x_pos)-model_line_width+jitter, x(x_pos)+model_line_width+jitter], ...
            [prior_means(x_pos), prior_means(x_pos)], ...
            'Color', [0 0 0], ...
            'Marker', 'none', ...
            'LineStyle', '--', ...
            'LineWidth', 2);

        %Winning model error bar
        line( ...
            [x(x_pos)+jitter, x(x_pos)+jitter], ...
            [prior_means(x_pos) - prior_cis(x_pos), prior_means(x_pos) + prior_cis(x_pos)], ...
            'Color', [0 0 0], ...
            'Marker', 'none', ...
            'LineStyle', '--', ...
            'LineWidth', 1.5);

    end

end

box off;
xlim([.5 2.5]);
ylim([-16 16]);

title(input_struct.figure_title, ...
    'FontName', 'Arial', ...
    'FontSize', input_struct.graph_font, ...
    'FontWeight', 'normal');

ylabel('Mean adjustment towards guilty');

set(gca, ...
    'YTick', [-15:3:15], ...
    'Xtick', [1 2], ...
    'Xticklabel', {'Innocent claims' 'Guilty claims'}, ...
    'Fontname', 'Arial', ...
    'Fontsize', input_struct.graph_font);

box off;

legend_labels = { ...
    input_struct.group_labels{1}, ...
    input_struct.group_labels{2} ...
    };

lgd = legend(b, legend_labels, 'Location', 'southwest');
legend boxoff;
lgd.FontName = 'Arial';
lgd.FontSize = input_struct.graph_font;

text(.75, 14.5, 'Black horizontal ticks = ground truth', ...
    'Color', [0 0 0], ...
    'FontName', 'Arial', ...
    'FontSize', input_struct.graph_font);

end
%%%%%%%%%%%%%%%%%end, make_grouped_bar_with_errors%%%%%%%%%%%%%%%%%%%%%%







%%%%%%%%%%%%%%%%%begin, make_prior_bar_with_errors%%%%%%%%%%%%%%%%%%%%%%
function sp_h = make_prior_bar_with_errors(data, input_struct)

% This plots sequence-position-0 probability ratings by suspect gender.
% This is not technically an "adjustment", because no preceding claim has
% occurred at sequence position 0. It is the initial/prior probability rating.

figure(input_struct.fig);

sp_h = subplot(input_struct.subplot(1), input_struct.subplot(2), input_struct.subplot(3));
hold on;

%Restrict to sequence position 0
prior_data = data(data.SequencePosition == 0, :);

%Get means over trials/sequences for each participant, then means and CIs over participants
prior_means_table = grpstats( ...
    grpstats(prior_data, {'Suspect', 'Pid'}, {'mean'}), ...
    {'Suspect'}, ...
    {'mean' 'meanci'} ...
    );

%The two rows should correspond to Suspect = 0 then Suspect = 1.
human_means = prior_means_table.mean_mean_HumanProbability;
human_cis = human_means - prior_means_table.meanci_mean_HumanProbability(:,1);

groundTruth_means = prior_means_table.mean_mean_GroundTruthProbability;
groundTruth_cis = groundTruth_means - prior_means_table.meanci_mean_GroundTruthProbability(:,1);

%Participant bars
b = bar(human_means);
hold on;

b.FaceColor = [1 1 1];
b.EdgeColor = [0 0 0];
b.LineWidth = 1.5;

x = b.XEndPoints;

for x_pos = 1:numel(x)

    %Participant error bars
    line( ...
        [x(x_pos) x(x_pos)], ...
        [human_means(x_pos) - human_cis(x_pos), human_means(x_pos) + human_cis(x_pos)], ...
        'Color', [0 0 0], ...
        'Marker', 'none', ...
        'LineStyle', '-', ...
        'LineWidth', 1.5);

    %Ground truth mean tick
    line( ...
        [x(x_pos)-0.075, x(x_pos)+0.075], ...
        [groundTruth_means(x_pos), groundTruth_means(x_pos)], ...
        'Color', [0 0 0], ...
        'Marker', 'none', ...
        'LineStyle', '-', ...
        'LineWidth', 2.5);

    %Ground truth error bar
    line( ...
        [x(x_pos)-0.025, x(x_pos)-0.025], ...
        [groundTruth_means(x_pos) - groundTruth_cis(x_pos), groundTruth_means(x_pos) + groundTruth_cis(x_pos)], ...
        'Color', [0 0 0], ...
        'Marker', 'none', ...
        'LineStyle', '-', ...
        'LineWidth', 1.5);

end

%Equal prior/reference line
h_prior = line([0.5 2.5], [50 50], ...
    'Color', [0 0 0], ...
    'LineStyle', ':', ...
    'LineWidth', 1);
h_prior.Annotation.LegendInformation.IconDisplayStyle = 'off';

box off;
xlim([.5 2.5]);
ylim([0 100]);

title(input_struct.figure_title, ...
    'FontName', 'Arial', ...
    'FontSize', input_struct.graph_font, ...
    'FontWeight', 'normal');

ylabel("'Guilt' probability at sequence position 0");

set(gca, ...
    'YTick', [0:20:100], ...
    'Xtick', [1 2], ...
    'Xticklabel', input_struct.group_labels, ...
    'Fontname', 'Arial', ...
    'Fontsize', input_struct.graph_font);

box off;

text(.65, 92, 'Black horizontal ticks = ground truth', ...
    'Color', [0 0 0], ...
    'FontName', 'Arial', ...
    'FontSize', input_struct.graph_font);

end
%%%%%%%%%%%%%%%%%end, make_prior_bar_with_errors%%%%%%%%%%%%%%%%%%%%%%







%%%%%%%%%%%%%%%%%%start, get_contexts%%%%%%%%%%%%%%%%%%%%%%
function contexts = get_contexts(raw)

%This loop takes each sequence, finds the cumulative sum/proportion of guilt
%claims, and creates:
% contexts(:,1): preceding context degree
% contexts(:,2): preceding context category

degree_half = NaN;

seq_starts = find(raw(:,5) == 0);

contexts = nan(size(raw,1), 2);

for sequence = 1:size(seq_starts,1)

    clear this_sequence_claims* temp_MG

    %extract claims for this sequence, positions 1 through 10
    this_sequence_claims = raw(seq_starts(sequence,1)+1 : seq_starts(sequence,1)+10, 8);

    %cumulative guilt claims and proportions
    this_sequence_claims_cumsum = cumsum(this_sequence_claims);
    this_sequence_claims_cumpro = this_sequence_claims_cumsum ./ [1:10]';

    %preceding context degree
    contexts(seq_starts(sequence,1):seq_starts(sequence,1)+10, 1) = ...
        [degree_half; degree_half; this_sequence_claims_cumpro(1:9,:)];

    %preceding context category
    temp_MG = NaN(11,1);

    for position = 1:9

        if this_sequence_claims_cumpro(position,1) == .5
            temp_MG(position+2) = NaN;
        elseif this_sequence_claims_cumpro(position,1) > .5
            temp_MG(position+2) = 1;
        else
            temp_MG(position+2) = 0;
        end

    end

    contexts(seq_starts(sequence,1):seq_starts(sequence,1)+10, 2) = temp_MG;

end

end
%%%%%%%%%%%%%%%%%%end, get_contexts%%%%%%%%%%%%%%%%%%%%%%







%%%%%%%%%start, get_sub_data%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function data = get_sub_data

% data_trunc.xlsx was formatted from data_exp_11596-v23_task-mwjx (1).csv.
%
%1: event index
%2: participant private id
%3: RT
%4: human probability
%5: sequence position
%6: suspect gender (1=female, 0=male)
%7: witness gender (1=female)
%8: guilty claim (1=guilty)
%9: context (mostly innocent / mostly guilty)

data = xlsread('C:\matlab_files\fiance\forensic_beads_pub_repo\Forensic-beads-paper-1\Study2\data_trunc.xlsx');

%At sequence position 0, fill missing condition labels using position 1.
nan_indices = find(data(:,5) == 0);

data(nan_indices,9) = data(nan_indices+1,9);
data(nan_indices,6) = data(nan_indices+1,6);

%Get context operationalised according to preceding claims only
[data(:,[10 11])] = get_contexts(data);

end
%%%%%%%%%end, get_sub_data%%%%%%%%%%%%%%%%%%%%%%%%%%%%










%%%%%%%%%%%%%%%%%%START, prob_guilt_groundTruth%%%%%%%%%%%%%%%%%%%%%%
function model_probabilities = prob_guilt_groundTruth(params, this_ps_suspect_data)

% params:
% [prior split response_bias response_noise]

prior = params(1);
split = params(2);
response_bias = params(3);
response_noise = params(4);

%indices where sequence position is 0
seq_start_indices = find(this_ps_suspect_data.SequencePosition == 0);

%initialise output
model_probabilities = nan(size(this_ps_suspect_data,1), 1);

%For each sequence, loop through sequence and get ground-truth probabilities
for seq = 1:numel(seq_start_indices)

    for claim = 1:11

        %current index into this_ps_suspect_data
        index = seq_start_indices(seq) + claim - 1;

        q = split;

        %number of guilt claims so far
        ng = sum(this_ps_suspect_data.Claim(seq_start_indices(seq)+1:index));

        %number of draws so far
        nd = claim - 1;

        %Bayesian posterior probability
        noiseless_p = ...
            1 / (1 + ((1-prior)/prior) * (q/(1-q))^(nd - 2*ng));

        %response mapping; for ground truth this is identity
        noise_p = response_bias + response_noise * noiseless_p;

        if noise_p <= 0
            noise_p = 0;
        elseif noise_p >= 1
            noise_p = 1;
        end

        model_probabilities(index,1) = noise_p;

    end

end

end
%%%%%%%%%%%%%%%%%%END, prob_guilt_groundTruth%%%%%%%%%%%%%%%%%%%%%%








%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [] = plot_gender_sequence_difference(data)

% Participant-only sequence-position plot split by suspect gender.
% This is like the old middle row of the original sequence plot:
%   left  = mostly innocent sequences
%   right = mostly guilty sequences
% with male/female suspect lines and significant differences marked.

[sus_means_p_ss, sus_means_ss, sus_means_m_ss, seq_key] = reformat_sequence_data(data);

graph_font = 12;

% Run paired t-tests comparing male vs female suspect at each sequence position
posthocs = get_sequence_ttests(sus_means_ss, seq_key);

% Get participant means and CIs over participants,
% grouped by Suspect and Context.
[sus_means, sus_cis, sus_grps] = grpstats( ...
    sus_means_ss, ...
    seq_key(:,[2 3]), ...
    {'mean' 'meanci' 'gname'} ...
    );

fig_h = figure;
set(gcf, 'Color', [1 1 1]);

sequence_type_string = {'Mostly innocent sequences' 'Mostly guilty sequences'};

markersize = 4;
asterisk_y_displacement = 5;
asterisk_size = 6;

for plot = 1:2

    subplot(1, 2, plot);
    hold on;

    title(sequence_type_string{plot}, ...
        'Fontname', 'Arial', ...
        'Fontsize', graph_font, ...
        'FontWeight', 'normal');

    % Rows in sus_means/sus_cis should be:
    % 1 = male suspect, mostly innocent
    % 2 = male suspect, mostly guilty
    % 3 = female suspect, mostly innocent
    % 4 = female suspect, mostly guilty
    %
    % For mostly innocent, use rows 1 and 3.
    % For mostly guilty, use rows 2 and 4.
    if plot == 1
        rows_to_plot = [1 3];
    else
        rows_to_plot = [2 4];
    end

    for r = 1:numel(rows_to_plot)

        row_index = rows_to_plot(r);

        % Jitter keeps male/female error bars and markers visually separate.
        if r == 1
            jitter = 0.05;
            line_style = '-';
            marker_style = 'o';
            marker_face = [1 1 1];     % open circle
            line_width = 1.5;
        else
            jitter = -0.05;
            line_style = '--';
            marker_style = 's';
            marker_face = [.75 .75 .75]; % filled grey square
            line_width = 1.5;
        end

        % Error bars
        for j = 1:size(sus_cis, 2)

            hl = line( ...
                [(j-1)+jitter (j-1)+jitter], ...
                [sus_cis(row_index,j,2) sus_cis(row_index,j,1)], ...
                'Color', [0 0 0], ...
                'Marker', 'none', ...
                'LineStyle', '-', ...
                'LineWidth', 1);
            hl.Annotation.LegendInformation.IconDisplayStyle = 'off';

        end

        % Means and connecting lines
        h_lines(r) = line([0:10] + jitter, sus_means(row_index,:), ...
            'Color', [0 0 0], ...
            'Marker', marker_style, ...
            'MarkerSize', markersize, ...
            'MarkerFaceColor', marker_face, ...
            'MarkerEdgeColor', [0 0 0], ...
            'LineStyle', line_style, ...
            'LineWidth', line_width);

    end

    % Significance asterisks for male vs female differences.
    % posthocs row 1 = mostly innocent; row 2 = mostly guilty.
    for j = 1:size(posthocs, 2)

        if posthocs(plot,j) == 1

            % Put asterisk slightly above the higher of the two CIs.
            male_row = rows_to_plot(1);
            female_row = rows_to_plot(2);

            y_ast = max([ ...
                sus_cis(male_row,j,2), ...
                sus_cis(female_row,j,2) ...
                ]) + asterisk_y_displacement;

            h_ast = line( ...
                j-1, ...
                y_ast, ...
                'Color', [0 0 0], ...
                'Marker', '*', ...
                'MarkerSize', asterisk_size, ...
                'LineStyle', 'none');
            h_ast.Annotation.LegendInformation.IconDisplayStyle = 'off';

        end

    end

    % Equal-prior reference line
    h_prior = line([0 10], [50 50], ...
        'Color', [0 0 0], ...
        'LineStyle', ':', ...
        'LineWidth', 1);
    h_prior.Annotation.LegendInformation.IconDisplayStyle = 'off';

    box off;

    set(gca, ...
        'Xtick', [0:10], ...
        'Xticklabel', {'0' '1' '2' '3' '4' '5' '6' '7' '8' '9' '10'}, ...
        'XTickLabelRotation', 0, ...
        'YTick', [0:20:100], ...
        'Fontname', 'Arial', ...
        'Fontsize', graph_font);

    ylim([0 100]);
    xlim([-1 10.5]);

    xlabel('Sequence position');

    if plot == 1
        ylabel("'Guilt' probability");
    end

    lgd = legend( ...
        h_lines, ...
        {'Male suspect', 'Female suspect'}, ...
        'Location', 'southwest');
    legend boxoff;
    lgd.FontSize = graph_font;
    lgd.FontName = 'Arial';

end

end
%%%%%%%%%%%%%%%%%%END, plot_gender_sequence_difference%%%%%%%%%%%%%%%%%%%%%%%%%