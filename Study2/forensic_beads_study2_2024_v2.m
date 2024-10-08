function [sus_means_p_ss sus_key_p_ss] = forensic_beads_study2_prior_model;

%forensic_beads_study2_prior_model.m: OK I radically chanmged things. forensic_beads_study2_2024 (i.e., v1) is
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

addpath(genpath('C:\matlab_files\fiance\forensic_beads_pub_repo\Forensic-beads-paper-1\FMINSEARCHBND'))
addpath(genpath('C:\matlab_files\fiance\forensic_beads_pub_repo\Forensic-beads-paper-1\klabhub-bayesFactor-3d1e8a5'))
addpath(genpath('C:\matlab_files\fiance\forensic_beads_pub_repo\Forensic-beads-paper-1\shaded_plots'))

%%%%%%
%Set-able stuff!

ground_truth_params = [0.5 0.5 0.7 0.7 0 1];

%initial value of free params
params(1) = .5; %prior, initialised to optimal value (ground truth of paradigm)
params(2) = .5; %prior, initialised to optimal value (ground truth of paradigm)
params(3) = .7;  %split term, initialised to optimal value
params(4) = .7;  %split term, initialised to optimal value
params(5) = 0;  %bias term, intialised to optimal value
% params(6) = 0;  %bias term, intialised to optimal value
params(6) = 1;  %noise term, initialised to optimal value
% params(8) = 1;  %noise term, intialised to optimal value
    
%fitted parameters constrained to be between these values
%technically all params will be free, but you can fix some by forcing their
%range to be one value only.
lower_bounds = [0 0 0.7 .7 0 0];
upper_bounds = [1 1  0.7  .7 1 1];
%     lower_bounds = [0 0 .5 .5 0 0 1 1];
%     upper_bounds = [1 1  1  1 0 0 1 1];
    


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

% data = sortrows(data,{'Pid','Suspect','Context'});
% data = sortrows(data,{'Pid'});


%Who are the participants?
participant_list = unique(data.Pid,'stable');    %Vitally important participant num order is maintained or it'll get mismatched with results arrays later!!!
num_participants = numel(participant_list);

%initialise data matrix to hold modelling results
% var_names = {'Pid','Loss','PriorMale','PriorFemale','SplitMale','SplitFemale','BiasMale','BiasFemale','NoiseMale','NoiseFemale'};
var_names = {'Pid','Loss','PriorMale','PriorFemale','SplitMale','SplitFemale','Bias','Noise'};

%initialise table to hold modelling results
model_fitting_results = array2table(nan(0,numel(var_names)), 'VariableNames',var_names);

%initialise vector to hold models' probabilities
model_behaviour_results = [];

%initialise vector to hold ground truth probabilities
ground_truth_behaviour_results = [];

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

    %That was prior model, now do the same for ground truth model
    ground_truth_behaviour_results = [ ...
        ground_truth_behaviour_results;
        get_behaviour_for_this_ps_sequences(this_ps_data,ground_truth_params) ...
        ];

end;    %loop through participants

%Add the accumulated model probabilities as column of data table
data.ModelProbability = model_behaviour_results;

%Add the accumulated model probabilities as column of data table
data.GroundTruthProbability = ground_truth_behaviour_results;

%add adjustment columns to data table
data = [ ...
    data, ...
    array2table(get_adjustments(data),'VariableNames',{'HumanAdjust','ModelAdjust','GroundTruthAdjust'}) ...
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


% %Bar plots of parameter values and loss
% make_parameter_plot(model_fitting_results);

%Kernel distribution plot of the prior parameter values
% make_kernel_dists(model_fitting_results);

plot_sequences(data);

% %Plot adjustments by suspect and context
plot_adjustments(data);

fprintf('');

disp('audi5000');
%%%%%%%%%%%%%%%%%%end, forensic_beads_study2_2024_v2 (main function body)%%%%%%%%%%%%%%%%%%%%%%











%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [] = plot_sequence_panels(plot_details)

figure(plot_details.fig);

line_colors = [1 1 1; 1 1 1; 0 0 0; 0 0 0];

y_pos = [30 90 25 80];
markersize = 4;
asterisk_y_displacement = 3;
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


    for i=its;    %male and female lines

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
                , 'Color',[0 0 0] ...
                , 'Marker','none' ...
                ,'LineStyle','-' ...
                ,'LineWidth',1 ...
                );
            hl.Annotation.LegendInformation.IconDisplayStyle = 'off';

            %asterisk, a little above the male (higher) error bar
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
            end;    %If plotting the top line (males)

        end; %Loop through sequence positions

        %Plot human means and lines between them
        line([0:10]+jitter, plot_details.means(i,:) ...
            , 'Color',[0 0 0] ...
            , 'Marker','o' ...
            ,'MarkerSize',markersize ...
            ,'MarkerFaceColor',line_colors(i,:) ...
            ,'MarkerEdgeColor',[0 0 0] ...
            ,'LineStyle','-' ...
            ,'LineWidth',1 ...
            );

    end;    %loop through male and female lines
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
    sprintf('Male, %s', plot_details.label) ...
    sprintf('Female, %s', plot_details.label) ...
    };
lgd = legend(legend_strs, 'Location', 'southeast');
legend boxoff;
lgd.FontSize = plot_details.graph_font;
lgd.FontName = 'Arial';
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%END, plot_sequence_panels%%%%%%%%%%%%%%%%%










%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [] = plot_sequences(data);

[sus_means_p_ss sus_means_ss sus_means_m_ss seq_key] = reformat_sequence_data(data);

%---------prepare sequence plots
graph_font = 12;
plot_details.fig = figure; set(gcf,'Color',[1 1 1]);
plot_details.figRows = 3;
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

%make pair of plots
plot_sequence_panels(plot_details);



%-------sequence plot for prior model performance

%get prior model averages and run posthocs on them
posthocs = get_sequence_ttests(sus_means_p_ss, seq_key);

%Now get means and cis over those participants for plots
[sus_means_p sus_cis_p sus_grps_p] = grpstats( sus_means_p_ss, seq_key(:,[2 3]),{'mean' 'meanci','gname'} );

plot_details.subplot_num = 5;
plot_details.means = sus_means_p;
plot_details.cis = sus_cis_p;
plot_details.posthocs = posthocs;
plot_details.label = 'prior';

%make pair of plots
plot_sequence_panels(plot_details);
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
%     this_params = [new_params(this_suspect_code+1) new_params(this_suspect_code + 3) new_params(this_suspect_code + 5) new_params(this_suspect_code + 7)];
    this_params = [new_params(this_suspect_code+1) new_params(this_suspect_code + 3) new_params(5:end)];
    
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








% %%%%%%%%%%%%%%DISCONFIRMATORY ADJUSTMENT PLOT BEGIN%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function plot_adjustments(data);

%setup plot
input_struct.fig = figure('Color',[1 1 1]);
input_struct.graph_font = 12;
input_struct.series_colours = [0.75 0.75 0.75; 0.5 0.5 0.5];
input_struct.model_colours = [0 0 1; 0 1 0];


%-------suspect*claim plot

input_struct.subplot = [2,1,1];

%Get means over trials for each participant then get means and cis over participants for errorbars
%input_struct.means now contains table with all the mean and ci data
input_struct.means = grpstats( ...
    grpstats(data,{'Claim', 'Suspect','Pid'},{'mean'}), ...
    {'Claim','Suspect'},{'mean' 'meanci'});

input_struct.legend = {'Man suspect', 'Woman suspect', 'Ground truth', 'Prior model'};

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

input_struct.legend = {'Preceding innocent context', 'Preceding guilty context', 'Ground truth', 'Prior model'}

b = make_grouped_bar_with_errors(input_struct);
% %%%%%%%%%%%%%%DISCONFIRMATORY ADJUSTMENT PLOT ENDS%%%%%%%%%%%%%%%%%%%%%%%%%%%%











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
    param_means.mean_SplitMale',
    param_means.mean_SplitFemale,
    param_means.mean_Bias',
    param_means.mean_Noise'];
%     param_means.mean_BiasMale',
%     param_means.mean_NoiseMale',


input_struct.input_ci = ...
    [param_means.meanci_PriorMale(:,2)',
    param_means.meanci_PriorFemale(:,2)',
    param_means.meanci_SplitMale(:,2)',
    param_means.meanci_SplitFemale(:,2)',
    param_means.meanci_Bias(:,2)',
    param_means.meanci_Noise(:,2)'] - ...
    input_struct.input_means;
%     param_means.meanci_BiasFemale(:,2)',
%     param_means.meanci_NoiseMale(:,2)',

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


%While I'm at it, get the grand standard deviation of the residuals for use
%in parameter recovery

%average all rows (sequences / one sequence per participant)
model_fitting_results.residual = sqrt( model_fitting_results.Loss );
model_fitting_results.Ones = ones(size(model_fitting_results,1),1);
stds = grpstats(model_fitting_results, 'Ones',{'std'});

disp(sprintf('std of all residuals is: %4.4f',stds.std_residual));
%%%%%%%%%%%%%%MAKE_PARAMETER_PLOT ENDS%%%%%%%%%%%%%%%%%%%%%%%%%%%%%





% 
% 
% 
% %%%%%%%%%%%%%%MAKE KERNEL DISTS BEGIN%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function make_kernel_dists(model_fitting_results);
% 
% cmap = [0 0 0;
%     .25 .25 .25
%     ];
% 
% 
% alphas = [0 .75];
% 
% h1 = figure('Color',[1 1 1]);
% 
% % suspects = unique(model_fitting_results.Suspect);
% % suspects_num = numel(suspects);
% 
% for suspect = [1, 2];  %suspect parameters loop
%     
%     %Kernel density
%     [kernel(:,suspect) xi(:,suspect)] = ...
%         ksdensity( ...
%         table2array(model_fitting_results(:,suspect+2)) ...
%         );
%     
%     plot_shaded(xi(:,suspect),kernel(:,suspect),'Alpha',alphas(suspect),'Color',cmap(suspect,:),'LineWidth',2);
%     hold on;
%     
% end;    %suspect parameters loop
% 
% box off
% xlabel('Estimated prior parameter value');
% ylabel('Probability (kernel) density');
% %%%%%%%%%%%%%%MAKE KERNEL DISTS ENDS%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 










%%%%%%%%%%%%%%%%%%start, get_model_loss%%%%%%%%%%%%%%%%%%%%%%
function loss = get_model_loss(params, this_ps_data);

fitted_probabilities = ...
    get_behaviour_for_this_ps_sequences(this_ps_data,params);
    
    %Loop through the trials (claims) for this suspect
    loss = 0;
    for trial = 1:size(this_ps_data,1);
        
        %squared error loss, accumulating over both suspects
        y_hat = fitted_probabilities(trial,end);
        y = this_ps_data.HumanProbability(trial,1);  %human / participant probability is col 4.
        loss = loss + (y_hat - y)^2;
        
    end;    %loop through trials, this suspect
%%%%%%%%%%%%%%%%%%end, get_model_loss%%%%%%%%%%%%%%%%%%%%%%











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
            [x(x_pos)-model_line_width x(x_pos)+model_line_width] ...
            ,[groundTruth_means(x_pos) groundTruth_means(x_pos)] ...
            , 'Color',input_struct.model_colours(1,:) ...
            , 'Marker','none' ...
            ,'LineStyle','-' ...
            ,'LineWidth',2 ...
            );
        
    %ground truth errorbar
    line( ...
        [x(x_pos) x(x_pos)] ...
        ,[groundTruth_means(x_pos) - groundTruth_cis(x_pos) groundTruth_means(x_pos) + groundTruth_cis(x_pos)] ...
        , 'Color',input_struct.model_colours(1,:) ...
        , 'Marker','none' ...
        ,'LineStyle','-' ...
        ,'LineWidth',2 ...
        );

    %prior model mean line
    line( ...
            [x(x_pos)-model_line_width x(x_pos)+model_line_width] ...
            ,[prior_means(x_pos) prior_means(x_pos)] ...
            , 'Color',input_struct.model_colours(2,:) ...
            , 'Marker','none' ...
            ,'LineStyle','-' ...
            ,'LineWidth',2 ...
            );
        
    %prior model errorbar
    line( ...
        [x(x_pos) x(x_pos)] ...
        ,[prior_means(x_pos) - prior_cis(x_pos) prior_means(x_pos) + prior_cis(x_pos)] ...
        , 'Color',input_struct.model_colours(2,:) ...
        , 'Marker','none' ...
        ,'LineStyle','-' ...
        ,'LineWidth',2 ...
        );

end;    %loop through claim*suspect conditions / pocitions on x axis (x_pos)

box off;
xlim([.5 2.5]);
ylim([-15 15]);
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
text(.75,9,input_struct.legend{3},'Color',input_struct.model_colours(1,:),'FontName','Arial','FontSize',input_struct.graph_font);
text(.75,6,input_struct.legend{4},'Color',input_struct.model_colours(2,:),'FontName','Arial','FontSize',input_struct.graph_font);

fprintf('')
%%%%%%%%%%%%%%%%%end, make_grouped_bar_with_errors%%%%%%%%%%%%%%%%%%%%%%







%%%%%%%%%%%%%%%%%%start, get_context%%%%%%%%%%%%%%%%%%%%%%
function contexts = get_contexts(raw);

%This next loop takes out each sequence, finds the cumulative sum of guilt
%claims for each sequence position and finds
%cumulative proportion of guilt claims

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
%%%%%%%%%%%%%%%%%%end, get_model_behaviour%%%%%%%%%%%%%%%%%%%%%%

