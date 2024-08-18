clear all;


%22-Aug 2023 This looks like a version of the code used to create
%the JEP:LMC paper figure. Creates the sequence position probability plots and 
%the suspect adjustment plot. Look for analysis_naina_use_for_context_2023.m
%for the disconfirmatory adjustment plot. Let's build off this to revise submission
%figures and analysis. Am putting this into newly-created repo. *_2022 Was
%originally in C:\matlab_files\fiance\forensics_beads_MvF_naina.

%I've modified analysis_rujiba_v9 to work with Naina's data. This
%should create sequence position graphs in separate figures and a
%gender*claim graph.

%v9: Modified version of v3 for next rewrite of paper after Ryan's comments. 
%I need the sequence figures and the religion bar graph but the
%claim*context graph now computes context incorrectly.
%I need to update how conditional probabilities are
%computed, and I need to relabel some plots with new terminology. 
%Hopefully can put new scatterplots in space below these four panels.

%analysis_rujiba_v3: This is Naina's analysis code for her study's figure
%forensics_beads_naina_analysis_v5.m but the input is adapted to accept
%Rujiba's slightly differently formatted data files.


% cd('C:\matlab_files\fiance\forensic_beads');

% %reorganise so like naina's code
% fnames = {...
%     'C:\matlab_files\fiance\forensic_beads\01_christi_mostly_innoce.xlsx'...
%     'C:\matlab_files\fiance\forensic_beads\02_atheist_mostly_innoce.xlsx'...
%     'C:\matlab_files\fiance\forensic_beads\03_christi_mostly_guilty.xlsx'...
%     'C:\matlab_files\fiance\forensic_beads\04_atheist_mostly_guilty.xlsx'...
%     };
% sus_rel_codes = [1 0 1 0];
% context_codes = [0 0 1 1];
% raw = [];
% 
% for file = 1:4;
%     
%     clear temp;
%     temp = xlsread(fnames{file});
%     
%     raw = [...
%         raw;
%         temp(:,6) ...                                   %1: event index
%         temp(:,1) ...                                   %2: participant private id
%         temp(:,2) ...                                   %3: RT
%         temp(:,3) ...                                   %4: probability estimate
%         repmat([0:10]',size(temp,1)/11 ,1) ...          %5: sequence position (including prior 0-10)
%         sus_rel_codes(file)*ones(size(temp,1),1) ...    %6: suspect religion (0=Atheist)
%         temp(:,5) ...                                   %7: witness gender (1=female)
%         temp(:,4) ...                                   %8: witness claim (1=guilt)
%         context_codes(file)*ones(size(temp,1),1) ...    %9: context (1=mosty guilty)
%         ];
%     
% end;    %Loop through datafiles (file)

raw = xlsread('C:\matlab_files\fiance\forensics_beads_MvF_naina\data_trunc.xlsx');


% 
% %recode columns 9 (context) so that it is 1 (mostly guilty) if more ones
% %than zeros or 0 if more zeros than ones or NaN is they are equal.
% %maybe not elegent but I'm adding this bit of code after everything else is
% %written so I'm going to leave that code as is and make this change right
% %from the start.
% seq_starts = find(raw(:,5) == 0);   %first element of each sequence
% for sequence=1:size(seq_starts,1); %for every start of a sequence
%     
%     clear num_guilts;
%     
%     %extract claims for this sequence (positions 1 through 10)
%     this_sequence_claims = raw(seq_starts(sequence,1)+1:seq_starts(sequence,1)+10,8);
%     %get cumsum of guilt claims (=1 each)
%     this_sequence_claims_cumsum = cumsum( this_sequence_claims);
%     this_sequence_claims_cumpro = this_sequence_claims_cumsum./[1:10]';
%     %recode as MI (0), MG (1) or neiher (NaN)
%     temp_MG = NaN(11,1);
%     for position = 1:9; %9, because the last position doesn't give a context to any folling claim / rating so should be left off
%         
%         if this_sequence_claims_cumpro(position,1) == .5;
%             
%             temp_MG(position+2) = NaN;
%             
%         elseif this_sequence_claims_cumpro(position,1) > .5;
%             
%             temp_MG(position+2) = 1;
%             
%         else;
%             
%             temp_MG(position+2) = 0;
%             
%         end;    %If/then to assign MG/MI values
%         
%     end;    %loop through this seq positions
%     
% end;    %loop seq starts




num_subs = numel(unique(raw(:,2)));
%Presumably there will be 11 data points for all eight sequences for each of the 104 participants in correct
%sequence order when raw is sorted by event index.
seq_pos = repmat([1:11]',num_subs*8,1);

%LMM

%In raw, sequence positions 0 have NaNs in place of condition labels for contexts (col 9). Put
%them back in or grpstats will exclude the prior from the average
nan_indices = find( raw(:,5) == 0);    %find NaNs
raw(nan_indices,9) = raw(nan_indices+1,9);  %assign the missing values at pos 0 with the values at pos 1

%Keeps all the trials
%subject, then suspect then context then sequence position
seq_out = [  raw(:,2) raw(:,6) raw(:,9) raw(:,5) raw(:,4) ];

%I commented out so future runs don't overwrite the header I manually put in the text file.
dlmwrite('C:\matlab_files\fiance\forensics_beads_MvF_naina\forensic_beads_rujiba_seqdata_forLMM_2021.txt',seq_out,' ');

%OK now try again but average over the two trials in each condition and
%array them in columns for repeated measures analysis in spss/jasp

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

dlmwrite('C:\matlab_files\fiance\forensics_beads_MvF_naina\forensic_beads_rujiba_seqdata_forANOVA_2021.txt',seq_cols,' ');


%Time for some figures

%Place here code from *_v4 that computes human and model adjustments

%on which indices is trhe display screen 0 (prior rating prompt so first
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
        %human - model adjustments (in 13)
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


% figure; set(gcf,'Color',[1 1 1]);
graph_font = 12;
% dot_colors = [0.5 0.5 0.5; 0 0 0];
dot_colors = [0.75 0.75 0.75; 0.5 0.5 0.5];
% model_color = [0 0 1];
model_color = [0 0 0];

x_locs = [1 4; 2 5];

% %ADJUSTMENTS: CONTEXT
% %this plot has separate bars for claims
% % subplot(1,4,1); hold on;

% p_y_pos = [1; -1];
% 
% % %get averages of subject (2), context (9) and claims (8)
% [ss ss_grps] = grpstats(raw(:,11),[raw(:,2) raw(:,9) raw(:,8)],{'mean' 'gname'});
% ss_grps = [str2num(cell2mat(ss_grps(:,2))) str2num(cell2mat(ss_grps(:,3)))  ]; %witnessMF then claimIG
% [h_mean h_ci grps] = grpstats(ss,[ss_grps(:,1) ss_grps(:,2)],{'mean' 'meanci','gname'});
% grps = [str2num(cell2mat(grps(:,1))) str2num(cell2mat(grps(:,2)))  ]; %witnessMF then claimIG
% h_ci = h_ci(:,2) - h_mean;
% 
% %Do t-tests between suspects for each witness claim to paint on graph
% [h_claimI,p_claim(1),ci,stats_claimI] = ...
%     ttest2(...
%     ss(find(ss_grps(:,1)==0 & ss_grps(:,2)==0),1), ...   %get data for mostly innocent, innocent claims
%     ss(find(ss_grps(:,1)==1 & ss_grps(:,2)==0),1) ...    %subtract sata for mostly guilty, innocent claims
%     ,'alpha',0.05/2);
% 
% [h_claimG,p_claim(2),ci,stats_claimG] = ...
%     ttest2(...
%     ss(find(ss_grps(:,1)==0 & ss_grps(:,2)==1),1), ...   %get data for mostly innocent, innocent claims
%     ss(find(ss_grps(:,1)==1 & ss_grps(:,2)==1),1) ...    %subtract sata for mostly guilty, innocent claims
%     ,'alpha',0.05/2);
% 
% 
% %And for model
% %get averages of claims and witness genders and subjects
% [ms ms_grps] = grpstats(raw(:,12),[raw(:,2) raw(:,9) raw(:,8)],{'mean' 'gname'});
% ms_grps = [str2num(cell2mat(ms_grps(:,2))) str2num(cell2mat(ms_grps(:,3)))  ]; %witnessMF then claimIG
% [m_mean m_ci m_grps] = grpstats(ms,[ms_grps(:,1) ms_grps(:,2)],{'mean' 'meanci','gname'});
% m_grps = [str2num(cell2mat(m_grps(:,1))) str2num(cell2mat(m_grps(:,2)))  ]; %witnessMF then claimIG
% % ci = ci(:,2) - mean;
% 
model_jitter = 0.1;
model_line_width = 0.4;
% 
% for context = 0:1;
%     for claim = 0:1;
%         
%         %         %stupid APA style p-value reporting!
%         %         if p_claim(claim+1) < .001;
%         %             tstr = '<.001';
%         %         elseif p_claim(claim+1) > .001 & p_claim(claim+1) < .01;
%         %             tstr = strrep(sprintf('%3.3f',p_claim(claim+1)),'0.','.');
%         %         else
%         %             tstr = strrep(sprintf('%3.2f',p_claim(claim+1)),'0.','.');
%         %         end;
%         %
%         %             text( ...
%         %                 mean( [x_locs(1,claim+1)  x_locs(2,claim+1)] ) ...
%         %                 , p_y_pos(claim+1) ...
%         %                 , tstr ...
%         %                 , 'Fontname','Arial' ...
%         %                 ,'Fontsize', graph_font ...
%         %                 ,'HorizontalAlignment','center' ...
%         %                 ,'VerticalAlignment','middle' ...
%         %                 );
%         
%         bar(x_locs(context+1,claim+1),h_mean(find(grps(:,1)==context & grps(:,2)==claim),1),...
%             'FaceColor',[1 1 1],'EdgeColor',dot_colors(context+1,:),'LineWidth',2);
%         
%         errorbar(x_locs(context+1,claim+1),...
%             h_mean(find(grps(:,1)==context & grps(:,2)==claim),1),...
%             h_ci(find(grps(:,1)==context & grps(:,2)==claim),1), ...
%             'Color',dot_colors(context+1,:),'LineStyle','none','LineWidth',2);
%         
%         %         handles = plotSpread(ss(find(ss_grps(:,1)==context & ss_grps(:,2)==claim),1), ...
%         %             'xValues',x_locs(context+1,claim+1),'distributionColors',dot_colors(context+1,:),'distributionMarkers','.');
%         %
%         %         plot(x_locs(context+1,claim+1),h_mean(find(grps(:,1)==context & grps(:,2)==claim),1),...
%         %            'LineStyle','none','Marker','o','markersize',6,'MarkerFaceColor',[1 1 1],'MarkerEdgeColor',[1 0 0]);
%         %
%         
%         %model points (Make it a horizontal line, since the plot function won't)
%         line( ...
%             [(x_locs(context+1,claim+1)+model_jitter)-model_line_width (x_locs(context+1,claim+1)+model_jitter)+model_line_width] ...
%             ,[m_mean(find(m_grps(:,1)==context & m_grps(:,2)==claim),1) m_mean(find(m_grps(:,1)==context & m_grps(:,2)==claim),1)] ...
%             , 'Color',model_color ...
%             , 'Marker','none' ...
%             ,'LineStyle','-' ...
%             ,'LineWidth',1 ...
%             );
%         %
%         %         plot( ...
%         %             x_locs(context+1,claim+1)+model_jitter ...
%         %             ,m_mean(find(m_grps(:,1)==context & m_grps(:,2)==claim),1),...
%         %             'LineStyle','none','Marker','-','markersize',12,'MarkerFaceColor',[1 1 1],'MarkerEdgeColor',model_color);
%         %
%         %model errors
%         line( ...
%             [x_locs(context+1,claim+1)+model_jitter x_locs(context+1,claim+1)+model_jitter] ...
%             ,[m_ci(find(m_grps(:,1)==context & m_grps(:,2)==claim),1) m_ci(find(m_grps(:,1)==context & m_grps(:,2)==claim),2)] ...
%             , 'Color',model_color ...
%             , 'Marker','none' ...
%             ,'LineStyle','-' ...
%             ,'LineWidth',2 ...
%             );
%         %
%         xlim([0 6])
%         ylim([-12 12]);
%         ylabel('Mean adjustment towards guilty');
%         set(gca ...
%             ,'Xtick',[1.5 4.5] ...
%             ,'Xticklabel',{sprintf('Innocent claims')  sprintf('Guilty claims')} ...
%             ,'YTick',[-100:2:100] ...
%             ,'Fontname','Arial' ...
%             ,'Fontsize',graph_font ...
%             );
%         box off;
%         
%         text(.25,11.5,'Mostly innocent','Color',dot_colors(1,:),'FontName','Arial','FontSize',graph_font);
%         text(.25,9.75,'Mostly guilty','Color',dot_colors(2,:),'FontName','Arial','FontSize',graph_font);
%         text(.25,8,'Optimal','Color',model_color,'FontName','Arial','FontSize',graph_font);
%         
%     end;    %claim
% end;    %context

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
        
        %                 %stupid APA style p-value reporting!
        %         if p_claim(claim+1) < .001;
        %             tstr = '{\it p} <.001';
        %         elseif p_claim(claim+1) > .001 & p_claim(claim+1) < .01;
        %             tstr = ['{\it p} = ' strrep(sprintf('%3.3f',p_claim(claim+1)),'0.','.')];
        %         else
        %             tstr = strrep(sprintf('%3.2f',p_claim(claim+1)),'0.','.');
        %         end;
        
        %             text( ...
        %                 mean( [x_locs(1,claim+1)  x_locs(2,claim+1)] ) ...
        %                 , p_y_pos(claim+1) ...
        %                 , tstr ...
        %                 , 'Fontname','Arial' ...
        %                 ,'Fontsize', graph_font ...
        %                 ,'HorizontalAlignment','center' ...
        %                 ,'VerticalAlignment','middle' ...
        %                 );
        
        bar(x_locs(suspect+1,claim+1),h_mean(find(grps(:,1)==suspect & grps(:,2)==claim),1),...
            'FaceColor',[1 1 1],'EdgeColor',dot_colors(suspect+1,:),'LineWidth',2);
        
        hold on;
        
        errorbar(x_locs(suspect+1,claim+1),...
            h_mean(find(grps(:,1)==suspect & grps(:,2)==claim),1),...
            h_ci(find(grps(:,1)==suspect & grps(:,2)==claim),1), ...
            'Color',dot_colors(suspect+1,:),'LineStyle','none','LineWidth',2);
        
        %         handles = plotSpread(ss(find(ss_grps(:,1)==suspect & ss_grps(:,2)==claim),1), ...
        %             'xValues',x_locs(suspect+1,claim+1),'distributionColors',dot_colors(suspect+1,:),'distributionMarkers','.');
        %
        %         plot(x_locs(suspect+1,claim+1),h_mean(find(grps(:,1)==suspect & grps(:,2)==claim),1),...
        %            'LineStyle','none','Marker','o','markersize',6,'MarkerFaceColor',[1 1 1],'MarkerEdgeColor',[1 0 0]);
        
        %         %model points
        %         plot(...
        %             x_locs(suspect+1,claim+1)+model_jitter ...
        %             ,m_mean(find(m_grps(:,1)==suspect & m_grps(:,2)==claim),1),...
        %             'LineStyle','none','Marker','.','markersize',20,'MarkerFaceColor',model_color,'MarkerEdgeColor',model_color);
        %
        %model points (Make it a horizontal line, since the plot function won't)
        line( ...
            [(x_locs(suspect+1,claim+1)+model_jitter)-model_line_width (x_locs(suspect+1,claim+1)+model_jitter)+model_line_width] ...
            ,[m_mean(find(m_grps(:,1)==suspect & m_grps(:,2)==claim),1) m_mean(find(m_grps(:,1)==suspect & m_grps(:,2)==claim),1)] ...
            , 'Color',model_color ...
            , 'Marker','none' ...
            ,'LineStyle','-' ...
            ,'LineWidth',1 ...
            );
        %
        %         plot( ...
        %             x_locs(context+1,claim+1)+model_jitter ...
        %             ,m_mean(find(m_grps(:,1)==context & m_grps(:,2)==claim),1),...
        %             'LineStyle','none','Marker','-','markersize',12,'MarkerFaceColor',[1 1 1],'MarkerEdgeColor',model_color);
        %
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


%Mostly innocent/guilt sequence timecourses
f_h1 = figure; set(gcf,'Color',[1 1 1]);
f_h2 = figure; set(gcf,'Color',[1 1 1]);

% line_colors = [0.25 0.25 0.25; 0.25 0.25 0.25; 0.5 0.5 0.5; 0.5 0.5 0.5];
line_colors = [1 1 1; 1 1 1; 0 0 0; 0 0 0];

%humans: get the suspects (M & F)
%I think seq_key 2 is suspect and 3 is context
[sus_means sus_grps] = grpstats( sequences_h, seq_key(:,[1 2 3]),{'mean' 'gname'} );

%While you have the subject data for the four conditions, run t-tests comparing males and females separately for MI and MG
[key_temp key_grps] = grpstats( seq_key, seq_key(:,[1 2 3]),{'mean' 'gname'} ); %This is a hack that gives me the group labels as numbers

M_MI = sus_means(find(key_temp(:,2)==0 & key_temp(:,3)==0),:);
F_MI = sus_means(find(key_temp(:,2)==1 & key_temp(:,3)==0),:);
Diff_MI = M_MI-F_MI;
[h_MI,p,ci,stats] = ttest(Diff_MI,0,'alpha',0.05/22);
% [h_MI,p,ci,stats] = ...
%     ttest(M_MI,F_MI,'alpha',0.05/22);

M_MG = sus_means(find(key_temp(:,2)==0 & key_temp(:,3)==1),:);
F_MG = sus_means(find(key_temp(:,2)==1 & key_temp(:,3)==1),:);
Diff_MG = M_MG-F_MG;
[h_MG,p,ci,stats] = ttest(Diff_MG,0,'alpha',0.05/22);
% [h_MG,p,ci,stats] = ttest2(...
%     M_MG,F_MG,'alpha',0.05/22);



[sus_means sus_cis sus_grps] = grpstats( sus_means, cell2mat(sus_grps(:,[2 3])),{'mean' 'meanci','gname'} );
% sus_cis = sus_cis(:,:,2) - mean(sus_cis,3);
[dummy grp_i] = sortrows(cell2mat(sus_grps));
sus_means = sus_means(grp_i,:);
sus_cis = sus_cis(grp_i,:,:);

%model: get the suspects (M & F)
[sus_means_m sus_grps_m] = grpstats( sequences_m, seq_key(:,[1 2 3]),{'mean' 'gname'} );
[sus_means_m sus_cis_m sus_grps_m] = grpstats( sus_means_m, cell2mat(sus_grps_m(:,[2 3])),{'mean' 'meanci','gname'} );
% sus_cis_m = sus_cis_m(:,:,2) - mean(sus_cis_m,3);
%Make sure rows are sorted by condition codes and not some random subject
%order as entered into grpstats
[dummy grp_i] = sortrows(cell2mat(sus_grps_m));
sus_means_m = sus_means_m(grp_i,:);
sus_cis_m = sus_cis_m(grp_i,:,:);

y_pos = [30 90 25 80];
markersize = 4;
asterisk_y_displacement = 3;
asterisk_size = 6;

for plot = 1:2;
    
    if plot == 1;
        its = [1 3];
%         subplot(1,2,2);
        figure(f_h1); hold on;
        title('Innocent sequences','Fontname','Arial','Fontsize',graph_font,'FontWeight','normal');
        limsy = [0 65];
        hold on;
        posthocs = h_MI;
        
    else;
        its = [2 4];
%         subplot(1,2,1);
        figure(f_h2); hold on;
        title('Guilty sequences','Fontname','Arial','Fontsize',graph_font,'FontWeight','normal');
        limsy = [40 100];
        hold on;
        posthocs = h_MG;
    end;
    
    for i=its;
        
        %To keep lines/points clear of each other I add a tiny jitter in x for
        %the human performance and in y dfor the model
        if i < 3;
            jitter = 0.05;
        else
            jitter = -0.05;
        end;
        
        
        
        %Plot human errors
        %use straight lines for errors
        %Since you're looping through sequence positions now
        %anyway, plot the significance asterisks whilst you're at it
        for j=1:size(sus_cis,2);
            
            %model error bars
            h2 = line( ...
                [(j-1)+jitter (j-1)+jitter] ...
                ,[sus_cis_m(i,j,2) sus_cis_m(i,j,1)] ...
                , 'Color',[0.65 0.65 0.65] ...
                , 'Marker','none' ...
                ,'LineStyle','-' ...
                ,'LineWidth',1 ...
                );
            h2.Annotation.LegendInformation.IconDisplayStyle = 'off';
            
            %error bar
            hl = line( ...
                [(j-1)+jitter (j-1)+jitter] ...
                ,[sus_cis(i,j,2) sus_cis(i,j,1)] ...
                , 'Color',[0 0 0] ...
                , 'Marker','none' ...
                ,'LineStyle','-' ...
                ,'LineWidth',1 ...
                );
            hl.Annotation.LegendInformation.IconDisplayStyle = 'off';
            
            %asterisk, a little above the male (higher) error bar
            if i == 1 | i == 2;
                if posthocs(j) == 1;
                    h2 = line( ...
                        (j-1)+jitter ...
                        , sus_cis(i,j,2) + asterisk_y_displacement ...
                        , 'Color',[1 0 0] ...
                        , 'Marker','*' ...
                        , 'MarkerSize',asterisk_size ...
                        ,'LineStyle','none' ...
                        );
                    h2.Annotation.LegendInformation.IconDisplayStyle = 'off';
                    
                end;    %If posthoc is sig
            end;    %If a male suspect
            
        end; %Loop thropugh sequence positions
        
        %Plot model lines
        line(0:10, sus_means_m(i,:)+jitter ...
            , 'Color',[0.65 0.65 0.65] ...
            , 'Marker','o' ...
            ,'MarkerSize',markersize-1 ...
            ,'MarkerFaceColor',line_colors(i,:) ...
            ,'MarkerEdgeColor',[0.65 0.65 0.65] ...
            ,'LineStyle',':' ...
            ,'LineWidth',1 ...
            ,'DisplayName','test' ...
            );
        
        
        %Plot human means and lines between them
        %         errorbar([0:10]+x_jitter, sus_means(i,:),sus_cis(1,:) ...
        line([0:10]+jitter, sus_means(i,:) ...
            , 'Color',[0 0 0] ...
            , 'Marker','o' ...
            ,'MarkerSize',markersize ...
            ,'MarkerFaceColor',line_colors(i,:) ...
            ,'MarkerEdgeColor',[0 0 0] ...
            ,'LineStyle','-' ...
            ,'LineWidth',1 ...
            );
        
    end;
    box off
    
    set(gca ...
        ,'Xtick',[0:10] ...
        ,'Xticklabel',{'0' '1' '2' '3' '4' '5' '6' '7' '8' '9' '10'} ...
        ,'YTick',[-100:10:100] ...
        ,'Fontname','Arial' ...
        ,'Fontsize',graph_font ...
        );
    %     ylim(limsy);
    xlim([-1 10.5])
    xlabel('Sequence position');
    ylabel('Mean probability of guilt');
    
end;

legend_strs = { ...
    'Male, optimal' ...
    'Male, participants' ...
    'Female, optimal' ...
    'Female, participants' ...
    };
lgd = legend(legend_strs);
legend boxoff;
lgd.FontSize = graph_font;
lgd.FontName = 'Arial';
% legend('test' ...
%     , '',  '', '', '', '', '', '', '', '', '', '','test' ...
%    , 'test' ...
%     , '',  '', '', '', '', '', '', '', '', '', '','test' ...
%     );

disp('audi5000');