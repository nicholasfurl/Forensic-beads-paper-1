function [] = forensic_beads_study2_2024;

%I've combined forensic_beads_naina_analysis_v10_2023,
%analysis_naina_use_for_context and
%forensic_beads_simplify_study2_split2Suspects to make new results figures
%to respond to JEP:LM&C reviews

%raw:
%1: event index
%2: participant private id
%3: RT
%4: human guilt probability
%5: sequence position (including prior 0-10)
%6: suspect gender (0=Male)
%7: witness gender (0=Male)
%8: witness claim (1=guilt)
%9: context (1=mosty guilty)
%Then in the code below we will create:
%10: ground truth guilt probability
%11: human adjustment
%12: ground truth adjustment
%13: adjustement "bias" (human - ground truth adjustment)
%14: preceding context (degree)
%15: preceding context (category)
raw = xlsread('C:\matlab_files\fiance\forensic_beads_pub_repo\Forensic-beads-paper-1\Study2\data_trunc.xlsx');

num_subs = numel(unique(raw(:,2)));

%Assign context condition to rows where Gorilla didn't fill it in (i.e., seq pos 0)
nan_indices = find( raw(:,5) == 0);    %find NaNs
raw(nan_indices,9) = raw(nan_indices+1,9);  %assign the missing values at pos 0 with the values at pos 1


%now get sub averages of probabilities for the four conditions at every sequence position
%subject, then suspect then context then sequence position
[seq_analysis_ss seq_analysis_ss_grps] = ...
    grpstats( [raw(:,2) raw(:,6) raw(:,9) raw(:,5) raw(:,4)], {  raw(:,2) raw(:,5) raw(:,6) raw(:,9) },{'mean', 'gname'});

sub_names = unique(seq_analysis_ss(:,1));
seq_cols = [];
for i=1:size(sub_names,1);
    this_sub_data = ...
        [ sub_names(i,:) seq_analysis_ss(find(seq_analysis_ss(:,1)==sub_names(i,:)),5)' ];
    seq_cols = [seq_cols; this_sub_data];
end;



%Time for some figures

%Place here code from *_v4 that computes human and model adjustments

%on which indices is the display screen 0 (prior rating prompt so first
%rating
seq_start_indices = find(raw(:,5)==0);

%For each start index, loop through sequence and get model predictions and
%adjustments for model and human
for seq = 1:numel(seq_start_indices);

    %initialize starting values (prompt screen) for model estimates (col 10) and
    %human (col 11) and model (col 12) adjustments
    raw(seq_start_indices,10) = 50;     %model assumes 50/50 chance
    raw(seq_start_indices,11) = NaN;    %the first (prior) rating cannot be adjusted by anything previous to it
    raw(seq_start_indices,12) = NaN;    %the first (prior) rating cannot be adjusted by anything previous to it

    %let's make a new dataset that includes the sequences for each
    %condition that I can plot later
    %dom1: sequence position, dim2:witness(M/F); dim3:suspect(M/F); dim4:claim(I/G); dim5:context(I/G)
    %....and assign the zeroth (prior) claim to it

    %raw probabilities, organized by sequence (assignment of prior)
    sequences_h(seq, 1) = raw(seq_start_indices(seq),4);   %assign data for prior rating
    sequences_m(seq, 1) = raw(seq_start_indices(seq),10);   %assign data for prior rating
    %adjustments, organized by sequence (assignment of prior)
    sequences_adj_h(seq, 1) = NaN;   %There is no prior adjustment
    sequences_adj_m(seq, 1) = NaN;   %There is no prior adjustment

    %data that are general to each sequence, just suspect and context
    %note: gorilla doesn't assign all the condition data to the prior
    %screen so it's extracted from the first sequence screen instead (+1)
    seq_key(seq,:) = [...
        raw(seq_start_indices(seq)+1,2) ...    %participant private id
        raw(seq_start_indices(seq)+1,6) ...    %suspectMF
        raw(seq_start_indices(seq)+1,9) ...    %contextIG
        ];

    %Loop through this sequence
    for claim = 2:11;

        pg=[];
        q=0.7;

        ng = sum( raw(seq_start_indices(seq)+1:seq_start_indices(seq)+claim-1,8) );
        nd = claim-1;
        raw(seq_start_indices(seq)+claim-1,10) = (1/(1 + (q/(1-q))^(nd-2*ng)))*100;

        %human adjustments (original human probs in 4, human adjustments in 11)
        raw(seq_start_indices(seq)+claim-1,11) = raw(seq_start_indices(seq)+claim-1,4) - raw(seq_start_indices(seq)+claim-2,4);
        %model adjustments (original model probs in 10, model adjustments in 12)
        raw(seq_start_indices(seq)+claim-1,12) = raw(seq_start_indices(seq)+claim-1,10) - raw(seq_start_indices(seq)+claim-2,10);
        %(human - model) adjustments (in 13)
        raw(seq_start_indices(seq)+claim-1,13) = raw(seq_start_indices(seq)+claim-1,11) - raw(seq_start_indices(seq)+claim-1,12);

        %Another array, this time data are organised into sequences so
        %sequence positions are easily indexable. Used mainly for timescale plots

        %assign sequence position (this is raw probability)
        sequences_h(seq, claim) =  raw(seq_start_indices(seq)+claim-1,4);   %assign data
        %model raw probabilites
        sequences_m(seq, claim) =  raw(seq_start_indices(seq)+claim-1,10);   %assign data
        %assign sequence position (human adjustments)
        sequences_adj_h(seq, claim) = raw(seq_start_indices(seq)+claim-1,11);   %assign data
        %assign sequence position (model adjustments)
        sequences_adj_m(seq, claim) = raw(seq_start_indices(seq)+claim-1,12);   %assign data

    end;    %loop through this sequence (claim)
end;    %loop through seq starts

%Alright so looks like the above code gives us four sequences* matrices
%that give us data by sequence position for human and model probabilities
%and adjustments. They're 832*11 so 832 = 2 suspect genders * 2 contexts *
%2 sequences * 104 participants and 11 = prior + 10 sequence positions. How
%do we know which of the 832 are what? We have a 832*3 sequence key that gives the
%participant, suspect and context of each row of sequences*.

%Add two more columns to raw to operationalise context according to
%preceding claims (degree and category)
[raw(:,[14 15])] = ...
    get_contexts(raw);








% figure; set(gcf,'Color',[1 1 1]);
graph_font = 12;
% dot_colors = [0.5 0.5 0.5; 0 0 0];
dot_colors = [0.75 0.75 0.75; 0.5 0.5 0.5];
% model_color = [0 0 1];
model_color = [0 0 0];

x_locs = [1 4; 2 5];

model_jitter = 0.1;
model_line_width = 0.4;

figure; set(gcf,'Color',[1 1 1]);
% subplot(1,4,2); hold on;

%get averages of claims and witness genders and subjects
[ss ss_grps] = grpstats(raw(:,11),[raw(:,2) raw(:,6) raw(:,8)],{'mean' 'gname'});
ss_grps = [str2num(cell2mat(ss_grps(:,2))) str2num(cell2mat(ss_grps(:,3)))  ]; %witnessMF then claimIG
[h_mean h_ci grps] = grpstats(ss,[ss_grps(:,1) ss_grps(:,2)],{'mean' 'meanci','gname'});
grps = [str2num(cell2mat(grps(:,1))) str2num(cell2mat(grps(:,2)))  ]; %witnessMF then claimIG
h_ci = h_ci(:,2) - h_mean;

%Do t-tests between suspects for each witness claim to paint on graph

%Do t-tests between suspects for each witness claim to paint on graph
[h_claimI,p_claim(1),ci,stats_claimI] = ...
    ttest2(...
    ss(find(ss_grps(:,1)==0 & ss_grps(:,2)==0),1), ...   %get data for mostly innocent, innocent claims
    ss(find(ss_grps(:,1)==1 & ss_grps(:,2)==0),1) ...    %subtract sata for mostly guilty, innocent claims
    ,'alpha',0.05/2);

[h_claimG,p_claim(2),ci,stats_claimG] = ...
    ttest2(...
    ss(find(ss_grps(:,1)==0 & ss_grps(:,2)==1),1), ...   %get data for mostly innocent, innocent claims
    ss(find(ss_grps(:,1)==1 & ss_grps(:,2)==1),1) ...    %subtract sata for mostly guilty, innocent claims
    ,'alpha',0.05/2);

%And for model
%get averages of claims and witness genders and subjects
[ms ms_grps] = grpstats(raw(:,12),[raw(:,2) raw(:,6) raw(:,8)],{'mean' 'gname'});
ms_grps = [str2num(cell2mat(ms_grps(:,2))) str2num(cell2mat(ms_grps(:,3)))  ]; %witnessMF then claimIG
[m_mean m_ci m_grps] = grpstats(ms,[ms_grps(:,1) ms_grps(:,2)],{'mean' 'meanci','gname'});
m_grps = [str2num(cell2mat(m_grps(:,1))) str2num(cell2mat(m_grps(:,2)))  ]; %witnessMF then claimIG
% ci = ci(:,2) - mean;

for suspect = 0:1
    for claim = 0:1;

        bar(x_locs(suspect+1,claim+1),h_mean(find(grps(:,1)==suspect & grps(:,2)==claim),1),...
            'FaceColor',[1 1 1],'EdgeColor',dot_colors(suspect+1,:),'LineWidth',2);

        hold on;

        errorbar(x_locs(suspect+1,claim+1),...
            h_mean(find(grps(:,1)==suspect & grps(:,2)==claim),1),...
            h_ci(find(grps(:,1)==suspect & grps(:,2)==claim),1), ...
            'Color',dot_colors(suspect+1,:),'LineStyle','none','LineWidth',2);

        %model points (Make it a horizontal line, since the plot function won't)
        line( ...
            [(x_locs(suspect+1,claim+1)+model_jitter)-model_line_width (x_locs(suspect+1,claim+1)+model_jitter)+model_line_width] ...
            ,[m_mean(find(m_grps(:,1)==suspect & m_grps(:,2)==claim),1) m_mean(find(m_grps(:,1)==suspect & m_grps(:,2)==claim),1)] ...
            , 'Color',model_color ...
            , 'Marker','none' ...
            ,'LineStyle','-' ...
            ,'LineWidth',1 ...
            );

        %model errors
        line( ...
            [x_locs(suspect+1,claim+1)+model_jitter x_locs(suspect+1,claim+1)+model_jitter] ...
            ,[m_ci(find(m_grps(:,1)==suspect & m_grps(:,2)==claim),1) m_ci(find(m_grps(:,1)==suspect & m_grps(:,2)==claim),2)] ...
            , 'Color',model_color ...
            , 'Marker','none' ...
            ,'LineStyle','-' ...
            ,'LineWidth',2 ...
            );

        %         xlim([0 6])
        %         ylim([-12 12]);
        ylabel('Mean adjustment towards guilty');
        set(gca ...
            ,'Xtick',[1.5 4.5] ...
            ,'Xticklabel',{sprintf('Innocent claims')  sprintf('Guilty claims')} ...
            ,'YTick',[-100:2:100] ...
            ,'Fontname','Arial' ...
            ,'Fontsize',graph_font ...
            );
        box off;

        text(.25,10,'Male suspect','Color',dot_colors(1,:),'FontName','Arial','FontSize',graph_font);
        text(.25,8.5,'Female suspect','Color',dot_colors(2,:),'FontName','Arial','FontSize',graph_font);
        text(.25,7,'Optimal','Color',model_color,'FontName','Arial','FontSize',graph_font);

    end;    %claim
end;    %suspect











%---------prepare sequence plots
plot_details.fig = figure; set(gcf,'Color',[1 1 1]);
plot_details.figRows = 3;
plot_details.figCols = 2;
plot_details.graph_font = graph_font;


%--------sequence plot for ground truth performance

%get participant averages (There is no variability from participant to participant because every one
%saw the same sequences so we skip the cis and post hocs)
[sus_means_m sus_grps_m] = grpstats( sequences_m, seq_key(:,[2 3]),{'mean' 'gname'} );

plot_details.subplot_num = 1;
plot_details.means = sus_means_m;
plot_details.cis = [];
plot_details.posthocs = [];
plot_details.label = 'ground truth';

%make pair of plots
plot_sequences(plot_details);



%------sequence plot for human performance

%get participant averages and run posthocs on them
[sus_means_ss sus_grps_ss] = grpstats( sequences_h, seq_key(:,[1 2 3]),{'mean' 'gname'} );
sus_key_ss = str2double(sus_grps_ss);
posthocs = get_sequence_ttests(sus_means_ss, sus_key_ss);

%Now get means and cis over those participants for plots
[sus_means sus_cis sus_grps] = grpstats( sus_means_ss, cell2mat(sus_grps_ss(:,[2 3])),{'mean' 'meanci','gname'} );

%now plot them
plot_details.subplot_num = 3;
plot_details.means = sus_means;
plot_details.cis = sus_cis;
plot_details.posthocs = posthocs;
plot_details.label = 'participants';

%make pair of plots
plot_sequences(plot_details);



%-------sequence plot for prior model performance

%get prior model averages and run posthocs on them
[sus_means_p_ss, sus_key_p_ss] = forensic_beads_study2_prior_model;
posthocs = get_sequence_ttests(sus_means_p_ss, sus_key_p_ss);

%Now get means and cis over those participants for plots
[sus_means_p sus_cis_p sus_grps_p] = grpstats( sus_means_p_ss, sus_key_p_ss(:,[2 3]),{'mean' 'meanci','gname'} );

plot_details.subplot_num = 5;
plot_details.means = sus_means_p;
plot_details.cis = sus_cis_p;
plot_details.posthocs = posthocs;
plot_details.label = 'prior';

%make pair of plots
plot_sequences(plot_details);

disp('audi5000');



%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%start, get_context%%%%%%%%%%%%%%%%%%%%%%
function contexts = get_contexts(raw);

%This next loop takes out each sequence, finds the cumulative sum of guilt
%claims for each sequence position

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












%% %%%%%%%%%%%%%%%%%%%%%%%%%%
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







%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [] = plot_sequences(plot_details)

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


