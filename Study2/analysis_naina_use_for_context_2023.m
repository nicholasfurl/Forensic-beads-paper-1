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
    raw(seq_starts(sequence,1):seq_starts(sequence,1)+10,9) = [NaN; NaN; this_sequence_claims_cumpro(1:9,:)]; %NEW RAW COL 9, CONTEXT DEGREE
    
    
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
%         
%         if switch_majority == 1;
%             
%             %There's no column that codes for the "true" majority in this
%             %updated code (it's always based on previous draws) so compute
%             %it manually now
%             if sum( raw(seq_start_indices(seq)+1:seq_start_indices(seq)+10,8) ) < 5;    %If there are fewer than five guilts
%                 
%                 ng = nd - ng;
%                 
%             end;
%             
%         end;    %swicth majority
        
        %assign model probability to column 11
        raw(seq_start_indices(seq)+claim-1,11) = (1/(1 + (q/(1-q))^(nd-2*ng)))*100;
        
        
        %         else
        %
        %             %get model prediction for every seq position
        %             %         if raw(seq_start_indices(seq)+1,9)==1; %Figure out if it's majority G or I. Careful of NaNs on the display screen
        %             ng = sum( raw(seq_start_indices(seq)+1:seq_start_indices(seq)+claim-1,8) );
        %             nd = claim-1;
        % %             nd = (claim-1)-ng;  %how many innocents so far?
        %             raw(seq_start_indices(seq)+claim-1,11) = (1/(1 + (q/(1-q))^(nd-2*ng)))*100;
        % %         else
        % %              nd = sum( raw(seq_start_indices(seq)+1:seq_start_indices(seq)+claim-1,8) );
        % %             ng = claim-1;
        % % %             ng = (claim-1)-nd;
        % %             raw(seq_start_indices(seq)+claim-1,10) = (1-(1/(1 + (q/(1-q))^(nd-2*ng))))*100;;
        % %         end;
        
        
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


%%%%%%%%%%%%%%%%
%%%%regression analysis of context degree. I've altered it to doit
%%%%participant by participant for now
%%%%%%%%%%%%%%

sub_nums = unique(raw(:,2));

for sub = 1:numel(sub_nums);
    
    %guilt claims
    clear mdl*
    this_sub_data_gc = raw( find( raw(:,2) == sub_nums(sub) & raw(:,8) == 1 ), [9 12 13] );  % contextdegree human adjust, model adjust
    mdl_h = fitlm(this_sub_data_gc(:,1),this_sub_data_gc(:,2));    %human adjust
    mdl_m = fitlm(this_sub_data_gc(:,1),this_sub_data_gc(:,3));    %model adjust
    raw( find( raw(:,2) == sub_nums(sub) & raw(:,8) == 1 ), 15) = table2array(mdl_h.Coefficients(2,1));
    raw( find( raw(:,2) == sub_nums(sub) & raw(:,8) == 1 ), 16) = table2array(mdl_m.Coefficients(2,1));
  
    %innocent claims
    clear mdl*
    this_sub_data_gc = raw( find( raw(:,2) == sub_nums(sub) & raw(:,8) == 0 ), [9 12 13] );  % contextdegree human adjust, model adjust
    mdl_h = fitlm(this_sub_data_gc(:,1),this_sub_data_gc(:,2));    %human adjust
    mdl_m = fitlm(this_sub_data_gc(:,1),this_sub_data_gc(:,3));    %model adjust
    raw( find( raw(:,2) == sub_nums(sub) & raw(:,8) == 0 ), 15) = table2array(mdl_h.Coefficients(2,1));
    raw( find( raw(:,2) == sub_nums(sub) & raw(:,8) == 0 ), 16) = table2array(mdl_m.Coefficients(2,1));
     
end;    %loop through subs


% %data for guilt claims
% adjs_Gclaims = raw( find( raw(:,8) == 1 ), [4 11 12 13] );   %human prob, model prob, human adjust, model adjust
% contexts_Gclaims = raw( find( raw(:,8) == 1 ), 9 );     %context degree for guilt claims
% analysis_labels = {'humProb' 'modProb' 'humAdj' 'modAdj'};
% 
% %run regressions and output
% disp('Guilt Claims:');
% for i=1:size(adjs_Gclaims,2);
%     
%     mdl = fitlm(contexts_Gclaims,adjs_Gclaims(:,i));
%     
%     disp(sprintf( ...
%         '%s b=%2.2f, r=%0.2f t=%2.2f p=%0.3f' ...
%         , analysis_labels{i} ...
%         , table2array(mdl.Coefficients(2,1)) ... %beta
%         , sqrt(mdl.Rsquared.Ordinary)*(table2array(mdl.Coefficients(2,3))/abs(table2array(mdl.Coefficients(2,3)))) ...
%         , table2array(mdl.Coefficients(2,3)) ... %t-value
%         , table2array(mdl.Coefficients(2,4)) ... %p-value
%         ));
%     yhat_Gclaims(:,i) = table2array(mdl.Coefficients(1,1)) + table2array(mdl.Coefficients(2,1))*[0 1];
%     
%     
% end;    %Loop through dependent variables (4: prob, adjust for hum, model)
% 
% 
% %data for innocent claims
% adjs_Iclaims = raw( find( raw(:,8) == 0 ), [4 11 12 13] );   %human prob, model prob, human adjust, model adjust
% contexts_Iclaims = raw( find( raw(:,8) == 0 ), 9 );     %context degree for Innocent claims
% 
% 
% %run regressions and output
% disp('Innocent Claims:');
% for i=1:size(adjs_Iclaims,2);
%     
%     mdl = fitlm(contexts_Iclaims,adjs_Iclaims(:,i));
%     
%     disp(sprintf( ...
%         '%s b=%2.2f, r=%0.2f t=%2.2f p=%0.3f' ...
%         , analysis_labels{i} ...
%         , table2array(mdl.Coefficients(2,1)) ... %beta
%         , sqrt(mdl.Rsquared.Ordinary)*(table2array(mdl.Coefficients(2,3))/abs(table2array(mdl.Coefficients(2,3)))) ...
%         , table2array(mdl.Coefficients(2,3)) ... %t-value
%         , table2array(mdl.Coefficients(2,4)) ... %p-value
%         ));
%     yhat_Iclaims(:,i) = table2array(mdl.Coefficients(1,1)) + table2array(mdl.Coefficients(2,1))*[0 1];
%     
%     
% end;    %Loop through dependent variables (4: prob, adjust for hum, model)
% 
% %scatterplots
% % figure; set(gcf,'Color',[1 1 1]);
% % 
graph_font = 12;
% % fig_colors = [.4 .4 .4; .8 .8 .8];
% % fig_colors = [1 0 0; 0 0 1];
% fig_colors = [1 0 1; 0 1 1];
% line_colors = [1 0 0; 0 0 1];
% % line_colors = fig_colors;
% % line_colors = [1 0 1; 0 1 1];
% % fig_colors = cool(2);
% graph_alpha = .1;
% marker_size = 12;
% marker_pale = .75;   %multiplies by fig_colors for the background dots to make them more pale than the regression line
% jitter_size = .005; %just to separate points along x over so slightly (variance of Gaussian noise)
% reg_line_width = 1;
% 
% 
% 
% %%%%%SAME SCATTERPLOTS BUT NOW PARTICIPANTS AND MODELS COMPARED ON SAME
% %%%%%GRAPHS AND GRAPHS SEPARATED BY CLAIM
% scatterplots
% figure; set(gcf,'Color',[1 1 1]);
% 
% Alright let's try a colored plot with human adjustments
% subplot(1,2,1);
% hold on;
% x_data = contexts_Gclaims;
% x_data = x_data+jitter_size.*randn(size(x_data,1),1);
% y_data = adjs_Gclaims(:,3);
% s = scatter(x_data,y_data,marker_size ...
%     ,'Marker','o' ...
%     ,'MarkerEdgeColor','none' ...
%     ,'MarkerFaceColor', marker_pale*fig_colors(1,:) ...
%     ,'MarkerFaceAlpha', graph_alpha ...
%     ,'MarkerEdgeAlpha',graph_alpha ...
%     );
%     ,'MarkerFaceColor', 'none' ...
%     ,'MarkerEdgeColor',fig_colors(2,:) ...
% 
% x_data = contexts_Gclaims;
% x_data = x_data+jitter_size.*randn(size(x_data,1),1);
% y_data = adjs_Gclaims(:,4);
% s = scatter(x_data,y_data,marker_size+5 ...
%     ,'Marker','o' ...
%     ,'MarkerEdgeColor','none' ...
%     ,'MarkerFaceColor', marker_pale*fig_colors(2,:) ...
%     ,'MarkerFaceAlpha', graph_alpha ...
%     ,'MarkerEdgeAlpha',graph_alpha ...
%     );
%     ,'MarkerFaceColor', 'none' ...
%     ,'MarkerEdgeColor',fig_colors(2,:) ...
% 
% plotnregression lines
% line([0 1],yhat_Gclaims([1 end],3), ...
%     'Color', line_colors(1,:), ...
%     'LineWidth',reg_line_width ...
%     );
% line([0 1],yhat_Gclaims([1 end],4), ...
%     'Color', line_colors(2,:), ...
%     'LineWidth',reg_line_width ...
%     );
% 
% xlim([0 1])
% ylim([-100 100]);
% ylabel('Adjustments, guilt claims');
% xlabel('Proportion guilt context');
% set(gca ...
%     ,'Xtick',[0:.1:1] ...
%     ,'YTick',[-100:25:100] ...
%     ,'Fontname','Arial' ...
%     ,'Fontsize',graph_font ...
%     );
% title('Guilt claims' ...
% ,'Fontname','Arial' ...
%     ,'Fontsize',graph_font ...
%     );
% 
% And innocent claims ....
% subplot(1,2,2);
% figure; set(gcf,'Color',[1 1 1]);
% hold on;
% x_data = contexts_Iclaims;
% x_data = x_data+jitter_size.*randn(size(x_data,1),1);
% y_data = adjs_Iclaims(:,3);
% s = scatter(x_data,y_data,marker_size ...
%     ,'Marker','o' ...
%     ,'MarkerEdgeColor','none' ...
%     ,'MarkerFaceColor', marker_pale*fig_colors(1,:) ...
%     ,'MarkerFaceAlpha', graph_alpha ...
%     ,'MarkerEdgeAlpha',graph_alpha ...
%     );
%     ,'MarkerFaceColor', 'none' ...
%     ,'MarkerEdgeColor',fig_colors(2,:) ...
% 
% x_data = contexts_Iclaims;
% x_data = x_data+jitter_size.*randn(size(x_data,1),1);
% y_data = adjs_Iclaims(:,4);
% s = scatter(x_data,y_data,marker_size ...
%     ,'Marker','o' ...
%     ,'MarkerEdgeColor','none' ...
%     ,'MarkerFaceColor', marker_pale*fig_colors(2,:) ...
%     ,'MarkerFaceAlpha', graph_alpha ...
%     ,'MarkerEdgeAlpha',graph_alpha ...
%     );
%     ,'MarkerFaceColor', 'none' ...
%     ,'MarkerEdgeColor',fig_colors(2,:) ...
% 
% plot regression lines
% line([0 1],yhat_Iclaims([1 end],3), ...
%     'Color', line_colors(1,:), ...
%     'LineWidth',reg_line_width ...
%     );
% line([0 1],yhat_Iclaims([1 end],4), ...
%     'Color', line_colors(2,:), ...
%     'LineWidth',reg_line_width ...
%     );
% 
% xlim([0 1])
% ylim([-100 100]);
% ylabel('Adjustment, innocent claims');
% xlabel('Proportion guilt context');
% set(gca ...
%     ,'Xtick',[0:.1:1] ...
%     ,'YTick',[-100:25:100] ...
%     ,'Fontname','Arial' ...
%     ,'Fontsize',graph_font ...
%     );
% 
% legend({'Participants' 'Optimal' 'Participants' 'Optimal'});
% legend box off;
% 
% 
% %%%%Now a reorganised and colored bar graph with majority preceding context.
figure; set(gcf,'Color',[1 1 1]);
% 
% Now the means
dot_colors = [0.5 0.5 0.5; 0 0 0];
dot_colors = [0.75 0.75 0.75; 0.5 0.5 0.5];
model_color = [0 0 1];
model_color = [0 0 0];
% 
% ADJUSTMENTS: CONTEXT
% this plot has separate bars for claims
% subplot(2,2,4); hold on;
x_locs = [1 4; 2 5];
p_y_pos = [1; -1];

%get averages of subject (2), context category (10) and claims (8)
[ss ss_grps] = grpstats(raw(:,12),[raw(:,2) raw(:,10) raw(:,8)],{'mean' 'gname'});
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
[ms ms_grps] = grpstats(raw(:,13),[raw(:,2) raw(:,10) raw(:,8)],{'mean' 'gname'});
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

%make output datafile with probability, adjustment, bias averages for further analysis, need cols: 
% 1: pid
% 2: human probability, witness claim i context i
% 3: human probability, witness claim i context g
% 4: human probability, witness claim g context i
% 5: human probability, witness claim g context g
% 6: human adjustment, witness claim i context i
% 7: human adjustment, witness claim i context g
% 8: human adjustment, witness claim g context i
% 9: human adjustment, witness claim g context g
% 10: model probability, witness claim i context i
% 11: model probability, witness claim i context g
% 12: model probability, witness claim g context i
% 13: model probability, witness claim g context g
% 14: model adjustment, witness claim i context i
% 15: model adjustment, witness claim i context g
% 16: model adjustment, witness claim g context i
% 17: model adjustment, witness claim g context g
% 18: bias, witness claim i context i
% 19: bias, witness claim i context g
% 20: bias, witness claim g context i
% 21: bias, witness claim g context g
% 22: human degree slope context guilt
% 23: human degree slope context innocent
% 24: human degree slope context guilt
% 25: human degree slope context innocent
% 26: suspect
% 27: Affiliation (I will need to add to exported datafile later)

%I know I probably did some of this above but let's just start from scratch
%to be organised and not get confused and make sure all things are
%processed in the same way in the same place

% %I can't just use grpstats because some participants might not have any
% %valid context categories for some witness claims (e.g., 768105 in rujiba's data) so have to do
% %it the complicated way (of course)
% out_ss = [];
% sub_nums = unique(raw(:,2));
% for sub=1:numel(sub_nums);
%     
%     this_in_data = raw(find(raw(:,2) == sub_nums(sub)),:);
%     this_out_data = [];
%     this_out_data(1,1) = sub_nums(sub);
%     
%     
%     
%     for claim = 0:1;
%         
%         for context = 0:1;
%             
%             %average the performance for this claim type and this context
%             %human prob, human adj, model prob, model adj, bias
%             assigned_data = this_in_data(find(this_in_data(:,8)==claim & this_in_data(:,10)==context),[4 12 11 13 14]);
%             if ~isempty(assigned_data);
%                 this_out_data =  ...
%                     [this_out_data nanmean(assigned_data,1)];
%             else
%                 this_out_data = [this_out_data nan(1,5)];
%             end;
%             
%         end;    %find contexts
%         
%         %Get slopes for the current claim. Above, when I computed slope, I
%         %wrote the same slope to every row with a claim to which that slope
%         %applies. So if we average all slopes for g or i claims for a
%         %participant, we'll be averageing the same number and will get the
%         %right number back again
%         this_out_data = ...
%             [this_out_data nanmean(this_in_data(find(this_in_data(:,8)==claim),[15 16]))];
%         
%     end;    %find claims
%     
%     %add suspect
%     out_ss(sub,:) = [this_out_data this_in_data(1,6)];
%     
% end;    %Loop through subs

% %reorder to match more sensible column order mentioned above
% out_ss = out_ss(:,[1 ... pid
%     2 7 14 19  ... participant probability
%     3 8 15 20  ... participant adjustment
%     4 9 16 21 ...  model probability
%     5 10 17 22 ... model adjustment
%     6 11 18 23 ... bias
%     12 24  ... slopes human
%     24 25  ... slopes model
%     26 ... suspect religion
%     ]);

%%%%The above creates an output file where claims are split between
%%%%preceding context columns and one of these columns often has a NaN
%%%%because of lack of data. This causes JASP to throw out these4 subjects
%%%%even in analyses excluding the preceding context variable and I can't
%%%%exclude preceding context variable in SPSS at all. So write another
%%%%file that outputs without any context so I can do the base
%%%%claim*Suspect analysis without data exclusion. Probably recapitulates
%%%%other previous outputs

%average whole raw dataset by witness claim (8) and participant (2) and
%suspect (6) and context (10)
context_ss = grpstats(raw,[ raw(:,2) raw(:,8) raw(:,6) raw(:,10)],{'mean'});
context_temp_mgg = context_ss( find(context_ss(:,6) == 0 & context_ss(:,10) == 1 & context_ss(:,8) == 1),[2 4 12 11 13 14 15 16 6]);    %male suspect, guilt context guilt claims
context_temp_mgi = context_ss( find(context_ss(:,6) == 0 & context_ss(:,10) == 1 & context_ss(:,8) == 0),[2 4 12 11 13 14 15 16 6]);    %male suspect, guilt context innocent claims
context_temp_mig = context_ss( find(context_ss(:,6) == 0 & context_ss(:,10) == 0 & context_ss(:,8) == 1),[2 4 12 11 13 14 15 16 6]);    %male suspect, innocent context guilt claims
context_temp_mii = context_ss( find(context_ss(:,6) == 0 & context_ss(:,10) == 0 & context_ss(:,8) == 0),[2 4 12 11 13 14 15 16 6]);    %male suspect, innocent context innocent claims
context_temp_fgg = context_ss( find(context_ss(:,6) == 1 & context_ss(:,10) == 1 & context_ss(:,8) == 1),[2 4 12 11 13 14 15 16 6]);    %female suspect, guilt context guilt claims
context_temp_fgi = context_ss( find(context_ss(:,6) == 1 & context_ss(:,10) == 1 & context_ss(:,8) == 0),[2 4 12 11 13 14 15 16 6]);    %female suspect, guilt context innocent claims
context_temp_fig = context_ss( find(context_ss(:,6) == 1 & context_ss(:,10) == 0 & context_ss(:,8) == 1),[2 4 12 11 13 14 15 16 6]);    %female suspect, innocent context guilt claims
context_temp_fii = context_ss( find(context_ss(:,6) == 1 & context_ss(:,10) == 0 & context_ss(:,8) == 0),[2 4 12 11 13 14 15 16 6]);    %female suspect, innocent context innocent claims
%in nocontext_temp: (1,10,19,28,37,46,55,64) participant,
%(2,11,20,29,38,47,56,65) prob, (3,12,21,30,39,48,57,66) adj,
%(4,13,22,31,40,49,58,67) model prob, (5,14,23,32,41,50,59,68) model adj,
%(6,15,24,33,42,51,60,69) bias, (-) slope, (-) model
%slope, (-) suspect (-)
context_temp = [context_temp_mgg context_temp_mgi context_temp_mig context_temp_mii context_temp_fgg context_temp_fgi context_temp_fig context_temp_fii];
out_context_ss = context_temp(:,[1 ... participant
    2 11 20 29 38 47 56 65 ...  prob (four conditions mgg mgi fig fii, suspect then context then claim)
    3 12 21 30 39 48 57 66 ...  adj
    4 13 22 31 40 49 58 67 ... model prob
    5 14 23 32 41 50 59 68 ... model adj
    6 15 24 33 42 51 60 69 ]); %bias

%stack claim columns for jasp/spss - need probability, adjustment and bias
nocontext_ss = grpstats(raw,[ raw(:,2) raw(:,8) raw(:,6)],{'mean'});
nocontext_temp_mg = nocontext_ss( find(nocontext_ss(:,6) == 0 & nocontext_ss(:,8) == 1),[2 4 12 11 13 14 15 16 6]);    %male suspect, guilt claims
nocontext_temp_mi = nocontext_ss( find(nocontext_ss(:,6) == 0 & nocontext_ss(:,8) == 0),[2 4 12 11 13 14 15 16 6]);    %male suspect, innocent claims
nocontext_temp_fg = nocontext_ss( find(nocontext_ss(:,6) == 1 & nocontext_ss(:,8) == 1),[2 4 12 11 13 14 15 16 6]);    %female suspect, guilt claims
nocontext_temp_fi = nocontext_ss( find(nocontext_ss(:,6) == 1 & nocontext_ss(:,8) == 0),[2 4 12 11 13 14 15 16 6]);    %female suspect, innocent claims
nocontext_temp = [nocontext_temp_mg nocontext_temp_mi nocontext_temp_fg nocontext_temp_fi];
%Now reorder cols (pid, slopes and suspect are redundant in the two halves)
%in nocontext_temp: (1,10,19,28) participant, (2,11,20,29) prob, (3,12,21,30) adj, (4,13,22,31) model prob, (5,14,23,32) model adj, (6,15,24,33) bias, (7,16,25,34) slope, (8,17,26,35) model slope, (9,18,27,36) suspect
out_nocontext_ss = nocontext_temp(:,[1 ... participant
    2 11 20 29 ...  prob (four conditions mg mi fg fi)
    3 12 21 30 ...  adj
    4 13 22 31 ... model prob
    5 14 23 32 ... model adj
    6 15 24 33 ]); %bias
    
%I ended up not using slope analysis for Rujiba. And for rujiba, it was
%more straightforward because there was one sequence per participant. Now
%we have a number of them, so we'd either have to loop through sequences
%per participant and average their slopes and then group-analyse averages
%or switch to linear mixed model approach. I'll leave it for now.
%     7 16 ...        male slopes (there's no g and i for slopes but there is male and female suspects)
%     8 17 ]);        % female slopes

%%%%%%%
%%WRITE OUTPUT FILES
%%%%%%%%


% %What do I have to do to get rid of the scientific notation in matlab?
% %dlmwrite writes the notation into the text.
% dlmwrite('C:\matlab_files\fiance\forensic_beads\anova_2022.txt',out_ss,' ');
% 
%Correctly writes out pids instead of sci nottation
xlswrite('C:\matlab_files\fiance\forensics_beads_MvF_naina\anova_2022.xls',out_context_ss);

%Might be used for LMM, especially a LMM with context degree as covariate
xlswrite('C:\matlab_files\fiance\forensics_beads_MvF_naina\LMM_degree_2022.xls',raw);

% %Can be used by base models with no context (so missing context conditions
% %won't get participants thrown out of analysis by JASP)
% xlswrite('C:\matlab_files\fiance\forensics_beads_MvF_naina\anova_nocontext_2022.xls',out_nocontext_ss);



disp('audi5000');