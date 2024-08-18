clear all;

%22-Aug 2023 This looks like a version of the code used to create
%the JEP:LMC paper figure. Creates just the adjustment disconfirmatory
%figure. Look for *_10_2023 for the other plots. 
%Let's build off this to revise submission
%figures and analysis. Am putting this into newly-created repo. *_2022 Was
%originally in C:\matlab_files\fiance\forensics_beads_MvF_naina.

%v10: builds on v7. Completes output file. Realised that I was regressing
%trials and not participants when do degree analysis. v10 starts to rectify
%that. I'll begin by outputting individual participant/model slopes for
%analysis later.

%v7: reapproaching v6 and hopefully revamping it with fresh eyes after a
%couple of months away from it. Makes bar graph of oddball effect and
%scatterplot of guilt degree effect for figures. Use v9 for sequence
%position and athiest bar graph figures.

%v6 continues down the road of using context preceding current claim. Now
%want to treat as continuous measure and also I want to use a separate
%anaysis of sequence with the standard definition of MI & MG, a separate
%pre-registered analysis of suspect without any context variable, and a
%continuous analysis of context, without those other variables. This
%program does continuous analysis of context.

%v5 redefines mostly innocent and mostly guilty according to the context
%preceding each rating

%analysis_rujiba_v3: This is Naina's analysis code for her study's figure
%forensics_beads_naina_analysis_v5.m but the input is adapted to accept
%Rujiba's slightly differently formatted data files.

%data_trunc.xlsx, I formated from data_exp_11596-v23_task-mwjx (1).csv,
%which Naina acquired from the Gorilla box for this part of the study flowchart.
%1: event index, 2:participant private id, 3:RT, 4: response 5: display (which witness is it?)
%6: suspect gender (1=female), 7:witness gender (1=female),
%8: guilty claim (1=guilty)
%we add later 9: context degree 10: context category 11: model probability
%12: adj humans 13: adj model 14: probability bias (human - model) 15: human context degree slope 16: model context degree slope
raw = xlsread('C:\matlab_files\fiance\forensics_beads_MvF_naina\data_trunc.xlsx');


%This next loop takes out each sequence, finds the cumulative sum of guilt
%claims for each sequence position and puts it in 9th col of raw, and finds
%cumulative proportion of guilt claims and puts it in 10th column of raw
seq_starts = find(raw(:,5) == 0);   %first element of each sequence
for sequence=1:size(seq_starts,1); %for every start of a sequence
    
    clear this_sequence_claims* temp_MG;
    
    %extract claims for this sequence (positions 1 through 10)
    this_sequence_claims = raw(seq_starts(sequence,1)+1:seq_starts(sequence,1)+10,8);
    %get cumsum of guilt claims (=1 each)
    this_sequence_claims_cumsum = cumsum( this_sequence_claims);
    this_sequence_claims_cumpro = this_sequence_claims_cumsum./[1:10]';
    %save proportions to array called raw
%     raw(seq_starts(sequence,1):seq_starts(sequence,1)+10,9) = [NaN; NaN; this_sequence_claims_cumpro(1:9,:)]; %NEW RAW COL 9, CONTEXT DEGREE
    
    
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
    
    raw(seq_starts(sequence,1):seq_starts(sequence,1)+10,10) = temp_MG; %NEW RAW COL 10, CONTEXT CATEGORY
    
end;    %loop seq starts


num_subs = numel(unique(raw(:,2)));

%on which indices is the display screen 0 (prior rating prompt so first
%rating
seq_start_indices = find(raw(:,5)==0);

%For each start index, loop through sequence and get model predictions and
%adjustments for model and human
for seq = 1:numel(seq_start_indices);
    
    %initialise starting values for sequences
    %I've reoganised, now 9: context degree 10: context category 11: model probability 12: adj humans 13: adj model 14: bias (human - model)
    
    raw(seq_start_indices,11) = 50;     %model assumes 50/50 chance
    raw(seq_start_indices,12) = NaN;    %the first (prior) human's rating cannot be adjusted by anything previous to it
    raw(seq_start_indices,13) = NaN;    %the first (prior) model's rating cannot be adjusted by anything previous to it
    
    %Loop through this sequence
    for claim = 2:11;
        
        %get model prediction for every seq position
        pg=[];
        q=0.7;
        
        %get number of guilts (i.e., the number of 1s)
        ng = sum( raw(seq_start_indices(seq)+1:seq_start_indices(seq)+claim-1,8) );
        %get number of draws so far
        nd = claim-1;
        
        %assign model probability to column 11
        raw(seq_start_indices(seq)+claim-1,11) = (1/(1 + (q/(1-q))^(nd-2*ng)))*100;        
        
        %4: human probs in 4 9: context degree 10: context category 11: model probability; 
        %NEW: 12: adj humans 13: adj model 14: bias (human - model)
        raw(seq_start_indices(seq)+claim-1,12) = raw(seq_start_indices(seq)+claim-1,4) - raw(seq_start_indices(seq)+claim-2,4);
        %model adjustments (original model probs in 10, model adjustments in 13)
        raw(seq_start_indices(seq)+claim-1,13) = raw(seq_start_indices(seq)+claim-1,11) - raw(seq_start_indices(seq)+claim-2,11);
        %human - model adjustments (in 14)
        raw(seq_start_indices(seq)+claim-1,14) = raw(seq_start_indices(seq)+claim-1,4) - raw(seq_start_indices(seq)+claim-1,11);
        
    end;    %loop through this sequence (claim)
end;    %loop through seq starts


%%%%WHAT WE'LL USE:
%I've reoganised, now 8: claim type 9: context degree 10: context category
%11: model probability 12: adj humans 13: adj model 14: bias (human -
%model) 15: human context degree slope 16: model context degree slope 


% %%%%Now a reorganised and colored bar graph with majority preceding context.
graph_font = 12;
figure; set(gcf,'Color',[1 1 1]);

% Now the means
dot_colors = [0.5 0.5 0.5; 0 0 0];
dot_colors = [0.75 0.75 0.75; 0.5 0.5 0.5];
model_color = [0 0 1];
model_color = [0 0 0];


% ADJUSTMENTS: CONTEXT
% this plot has separate bars for claims
% subplot(2,2,4); hold on;
x_locs = [1 4; 2 5];
p_y_pos = [1; -1];

%get averages of subject (2), context category (10) and claims (8)
which_context = 9;
[ss ss_grps] = grpstats(raw(:,12),[raw(:,2) raw(:,which_context) raw(:,8)],{'mean' 'gname'});
ss_grps = [str2num(cell2mat(ss_grps(:,2))) str2num(cell2mat(ss_grps(:,3)))  ]; %witnessMF then claimIG
[h_mean h_ci grps] = grpstats(ss,[ss_grps(:,1) ss_grps(:,2)],{'mean' 'meanci','gname'});
grps = [str2num(cell2mat(grps(:,1))) str2num(cell2mat(grps(:,2)))  ]; %witnessMF then claimIG
h_ci = h_ci(:,2) - h_mean;

% Do t-tests between suspects for each witness claim to paint on graph
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


% And for model
% get averages of claims and witness genders and subjects
[ms ms_grps] = grpstats(raw(:,13),[raw(:,2) raw(:,which_context) raw(:,8)],{'mean' 'gname'});
ms_grps = [str2num(cell2mat(ms_grps(:,2))) str2num(cell2mat(ms_grps(:,3)))  ]; %witnessMF then claimIG
[m_mean m_ci m_grps] = grpstats(ms,[ms_grps(:,1) ms_grps(:,2)],{'mean' 'meanci','gname'});
m_grps = [str2num(cell2mat(m_grps(:,1))) str2num(cell2mat(m_grps(:,2)))  ]; %witnessMF then claimIG
% ci = ci(:,2) - mean;

model_jitter = 0.1;
model_line_width = 0.4;

for context = 0:1;
    for claim = 0:1;
        hold on;
        
        bar(x_locs(context+1,claim+1),h_mean(find(grps(:,1)==context & grps(:,2)==claim),1),...
            'FaceColor',[1 1 1],'EdgeColor',dot_colors(context+1,:),'LineWidth',2);
        
        errorbar(x_locs(context+1,claim+1),...
            h_mean(find(grps(:,1)==context & grps(:,2)==claim),1),...
            h_ci(find(grps(:,1)==context & grps(:,2)==claim),1), ...
            'Color',dot_colors(context+1,:),'LineStyle','none','LineWidth',2);
        
%         model points (Make it a horizontal line, since the plot function won't)
        line( ...
            [(x_locs(context+1,claim+1)+model_jitter)-model_line_width (x_locs(context+1,claim+1)+model_jitter)+model_line_width] ...
            ,[m_mean(find(m_grps(:,1)==context & m_grps(:,2)==claim),1) m_mean(find(m_grps(:,1)==context & m_grps(:,2)==claim),1)] ...
            , 'Color',model_color ...
            , 'Marker','none' ...
            ,'LineStyle','-' ...
            ,'LineWidth',1 ...
            );
        
%         model errorbars
        line( ...
            [x_locs(context+1,claim+1)+model_jitter x_locs(context+1,claim+1)+model_jitter] ...
            ,[m_ci(find(m_grps(:,1)==context & m_grps(:,2)==claim),1) m_ci(find(m_grps(:,1)==context & m_grps(:,2)==claim),2)] ...
            , 'Color',model_color ...
            , 'Marker','none' ...
            ,'LineStyle','-' ...
            ,'LineWidth',2 ...
            );
        
        xlim([0 6])
%                 ylim([-12 12]);
%                 ylim([-100 100]);
        ylabel('Mean adjustment towards guilty');
        set(gca ...
            ,'Xtick',[1.5 4.5] ...
            ,'Xticklabel',{sprintf('Innocent claims')  sprintf('Guilty claims')} ...
            ,'Fontname','Arial' ...
            ,'Fontsize',graph_font ...
            );
%                             ,'YTick',[-100:2:100] ...
        box off;
        
        text(.25,15,'Preceding innocent context','Color',dot_colors(1,:),'FontName','Arial','FontSize',graph_font);
        text(.25,13,'Preceding guilty context','Color',dot_colors(2,:),'FontName','Arial','FontSize',graph_font);
        text(.25,11,'Optimal','Color',model_color,'FontName','Arial','FontSize',graph_font);
        
    end;    %claim
end;    %context

disp('audi5000');