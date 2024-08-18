
%%
function [] = forensic_beads_simulation_delta;

%First make some stimuli for the model simulation

% Sequences = 1;  %if 1, creates random sequences to test, if 2 uses same 8 sequences as used in Study 2 (gender study)
sequence_compute = 1;  %if 1, creates random sequences to test, if 2 uses same 8 sequences as used in Study 2 (gender study)

if sequence_compute == 1;

    % Define parameters
    num_majority = 7;
    sequence_length = 10; % Length of each vector
    q = num_majority / sequence_length;    %probably of guilt
    num_participants = 600; % Total number of participants

    % Set the random seed
    seed = 42;
    rng(seed);

    % Number of vectors with 6 ones and 4 zeros and vice versa
    half_participants = num_participants / 2;

    % Preallocate the array to hold the vectors
    vectors = zeros(num_participants, sequence_length);

    for i = 1:half_participants

        % Create a vector with 6 ones and 4 zeros
        vec_ones = [ones(1, num_majority), zeros(1, sequence_length-num_majority)];
        % Shuffle the vector
        sequences(i, :) = vec_ones(randperm(sequence_length));

        % Create a vector with 6 zeros and 4 ones
        vec_zeros = [zeros(1, num_majority), ones(1, sequence_length-num_majority)];
        % Shuffle the vector
        sequences(i + half_participants, :) = vec_zeros(randperm(sequence_length));

    end;    %ends loop through simulated participants

    %Add a Nan to the front so prior can be processed same way as in
    %analysis of participants
    sequences = [nan(size(sequences,1),1) sequences];

elseif sequence_compute == 2;   %If you want to use same sequences as in study 2

    num_majority = 7;
    sequence_length = 10; % Length of each vector
    q = num_majority / sequence_length;    %probably of guilt

    sequences = [ ...
        NaN     1     1     0     1     1     1     0     1     1     0
        NaN     0     1     1     0     1     1     0     1     1     1
        NaN     1     0     1     0     1     1     0     1     1     1
        NaN     0     0     1     0     1     0     0     0     1     0
        NaN     0     1     0     0     0     0     1     0     1     0
        NaN     0     1     1     1     1     0     0     1     1     1
        NaN     0     0     1     0     1     0     0     1     0     0
        NaN     0     0     0     0     1     0     1     1     0     0
        ];

    num_participants = size(sequences,1); % Total number of participants

end;

%Get model probabilities for stimuli just created

%model parameters if using delta model
alpha = .9;
beta = 1;

%model parameters if using window model
prior = 0.5;
window = 3;

% Preallocate arrays to hold the average probabilities
avg_adj_Gclaim_Gseq = [];
avg_adj_Iclaim_Gseq = [];
avg_adj_Gclaim_Iseq = [];
avg_adj_Iclaim_Iseq = [];

for participant = 1:num_participants;

    this_sequence = sequences(participant,:);

%     [this_prob this_adj] = prob_guilt_nowindow(this_sequence,q, prior);
%         [this_prob this_adj] = prob_guilt_primacy(this_sequence,q,prior,window);
%         [this_prob this_adj] = prob_guilt_recency(this_sequence,q,prior,window);
        [this_prob this_adj] = prob_guilt_delta(this_sequence,alpha,beta,prior);

    %To make probabilities match ratings
    this_prob = this_prob*100;
    this_adj = this_adj*100;


    %set up some averages for plots later
    if nansum(this_sequence) < sequence_length / 2

        % Mostly innocent
        avg_adj_Gclaim_Iseq(end+1) = nanmean(this_adj(this_sequence == 1));
        avg_adj_Iclaim_Iseq(end+1) = nanmean(this_adj(this_sequence == 0));

    else

        % Mostly guilty
        avg_adj_Gclaim_Gseq(end+1) = nanmean(this_adj(this_sequence == 1));
        avg_adj_Iclaim_Gseq(end+1) = nanmean(this_adj(this_sequence == 0));

    end

    prob_guilt(participant,:) = this_prob;
    adj_guilt(participant,:) = this_adj;

end;    %loop through participants

% Prepare data for bar plot
x = categorical({'Guilty', 'Innocent'});
x = reordercats(x, {'Innocent', 'Guilty'});
y = [
    nanmean(avg_adj_Gclaim_Iseq), nanmean(avg_adj_Iclaim_Iseq);
    nanmean(avg_adj_Gclaim_Gseq), nanmean(avg_adj_Iclaim_Gseq)
]';

% Create the bar plot
figure;
bar(x, y, 'grouped');
legend({'Mostly Innocent','Mostly Guilty'});
ylabel('Average Adjustment');

fprintf('');








%%
%to be used separately for each sequence
function [prob_guilt,adj_guilt] = prob_guilt_delta(draws,alpha,beta,prior)

num_draws = numel(draws);

q_hat = prior;
q_hat_list = nan(1,num_draws);
q_hat_list(1,1) = q_hat;

prob_guilt = nan(1,num_draws);

for t = 1:num_draws;

    if ~isnan(draws(t));

        q_hat = q_hat + alpha * (draws(t) - q_hat);
        q_hat_list(1,t) = q_hat;

    end;    %keep qhat as prior if nan

    prob_guilt(1,t) = exp(beta*q_hat)/(exp(beta*q_hat) + exp(beta*(1-q_hat)));

end;    %loop through draws

adj_guilt = [NaN diff(prob_guilt)];



%%
function [prob_guilt,adj_guilt] = prob_guilt_recency(draws,q,prior,window);

num_draws = numel(draws);

for t = 1:num_draws;

    %change the actual sequence up to t
    draws_recency = draws(max(1, (t+1) - window):t);

%     %leading nan doesn't count
%     ng = sum(draws_recency(2:end));
%     num_samples = numel(draws_recency(2:end));

        %leading nan doesn't count
        ng = nansum(draws_recency);
        num_samples = sum(~isnan(draws_recency));

    prob_guilt(1,t) = prior / (prior + (1-prior)*(q/(1-q))^(num_samples-2*ng));

end;    %loop through draws

adj_guilt = [NaN diff(prob_guilt)];









%%
function [prob_guilt,adj_guilt] = prob_guilt_primacy(draws,q,prior,window);

N = window; %number of items at beginning of sequence, starting with the first sample (ignoring the prior nan)
num_draws = numel(draws);

for t = 1:num_draws;

    clear draws_primacy

    %change the actual sequence so that samples between N and t are dropped out
    if t > N + 1;   %+1 because we want to ignore the leading nan
        draws_primacy = [draws(1:1+N), draws(t)];   %The second sequence item is the first sample after the nan prior, then N past that prior nan
    else
        draws_primacy = draws(1:t); % If inside window, just count everything 
    end
    % You can now use draws_primacy for your computations

    %leading nan doesn't count
    ng = sum(draws_primacy(2:end));
    num_samples = numel(draws_primacy(2:end));

    prob_guilt(1,t) = prior / (prior + (1-prior)*(q/(1-q))^(num_samples-2*ng));

end

adj_guilt = [NaN diff(prob_guilt)];







%%
function [prob_guilt,adj_guilt] = prob_guilt_nowindow(draws,q, prior);

num_draws = numel(draws);

for t = 1:num_draws;

    %ignore the nan
    ng = sum(draws(2:t));
    nd = t-1;

    %conditional probability
    prob_guilt(t,1) = (1/(1 + ((1-prior)/prior)*(q/(1-q))^(nd-2*ng)));

end;    %loop through draws

adj_guilt = [NaN diff(prob_guilt)'];





