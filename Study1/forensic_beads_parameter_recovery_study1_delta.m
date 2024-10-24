%%%%%%%%%%%%%%%%%%start, forensic_beads_parameter_recovery%%%%%%%%%%%%%%%%%%%%%%
function forensic_beads_parameter_recovery_study2_delta;


%forensic_beads_parameter_recovery_study1_delta: I modified
%forensic_beads_parameter_recovery_study2_delta with input from
%forensic_beads_parameter_recovery_study1_simple

%forensic_beads_parameter_recovery_study2_delta: I modified
% forensic_beads_parameter_recovery_study2_simple to do parameter recovery
% with the delta model.

%A simpler, rewritten, debugged and double checked version of parameter recovery for
%study 2, fitting different priors to male and female sequences
%within-participant but one split, bias and noise parameter per participant.

addpath(genpath('C:\matlab_files\fiance\forensic_beads_pub_repo\Forensic-beads-paper-1\FMINSEARCHBND'));
addpath(genpath('C:\matlab_files\fiance\parameter_recovery\beta_fixed_code\Model_fitting_hybrid_study\plotSpread'));

%Settable things

%I'm disabling this and going all participanmt-specific for now
% %1 if parameters derived from empirical study
% %2 if small test grid search
% %3 if large grid search
% config_list = 1;

% num_reps = 10;  %Will get this number of random stimuli and create a "sample" of them for each combo of config_list == 2 or 3

% %random seed
% rng(42);

%Same as used in Study 1
%In study 1, there was one sequence for each of the 599 participants
stimuli = [ ...
   NaN     1     0     1     0     1     0     0     0     0     1
   NaN     0     1     0     1     0     0     0     1     1     0
   NaN     1     0     0     1     1     0     0     0     0     1
   NaN     0     1     0     0     0     1     0     1     1     0
   NaN     1     1     0     0     0     0     0     1     0     1
   NaN     1     1     0     0     0     0     0     1     0     1
   NaN     0     1     0     0     1     1     0     1     0     0
   NaN     0     1     0     1     1     0     0     0     0     1
   NaN     0     1     1     0     0     0     0     1     1     0
   NaN     0     1     0     0     0     1     1     0     0     1
   NaN     0     1     0     1     0     1     0     1     0     0
   NaN     1     0     0     1     1     0     0     0     1     0
   NaN     0     0     0     0     1     0     1     1     0     1
   NaN     0     0     1     0     0     1     1     0     0     1
   NaN     1     1     0     1     0     0     0     0     0     1
   NaN     0     0     1     0     0     0     1     1     1     0
   NaN     1     0     0     0     1     0     1     1     0     0
   NaN     1     0     0     0     1     1     0     1     0     0
   NaN     1     0     0     1     1     0     1     0     0     0
   NaN     0     0     0     0     1     0     1     1     1     0
   NaN     0     1     0     0     1     0     1     0     1     0
   NaN     0     0     1     0     1     0     1     0     0     1
   NaN     0     1     0     0     1     0     1     0     0     1
   NaN     0     0     1     1     0     0     1     0     0     1
   NaN     0     0     1     0     1     1     0     0     0     1
   NaN     1     0     1     1     0     0     1     0     0     0
   NaN     0     0     0     1     0     0     1     1     1     0
   NaN     0     0     1     0     1     1     0     1     0     0
   NaN     1     1     0     0     0     0     0     1     0     1
   NaN     1     0     1     1     0     0     0     0     0     1
   NaN     1     0     0     1     0     0     1     1     0     0
   NaN     0     0     0     0     1     0     0     1     1     1
   NaN     0     0     1     0     1     1     0     0     1     0
   NaN     1     0     0     0     0     1     0     1     1     0
   NaN     0     0     1     0     0     0     0     1     1     1
   NaN     0     0     1     1     1     0     1     0     0     0
   NaN     0     0     0     1     1     0     1     1     0     0
   NaN     1     0     0     0     1     0     0     1     0     1
   NaN     1     0     0     0     0     1     0     0     1     1
   NaN     0     0     0     0     0     1     0     1     1     1
   NaN     0     0     0     1     1     1     0     0     1     0
   NaN     1     1     0     1     0     0     0     1     0     0
   NaN     0     0     0     1     1     0     1     1     0     0
   NaN     1     0     0     0     0     1     0     1     1     0
   NaN     0     0     0     1     1     0     0     1     1     0
   NaN     1     0     0     1     0     0     1     0     1     0
   NaN     0     1     0     0     0     0     1     1     1     0
   NaN     0     1     0     0     0     0     0     1     1     1
   NaN     0     0     0     0     0     1     0     1     1     1
   NaN     0     0     1     1     0     0     0     0     1     1
   NaN     1     0     0     1     0     0     0     1     0     1
   NaN     0     1     0     1     1     0     0     1     0     0
   NaN     0     0     1     0     0     0     0     1     1     1
   NaN     1     0     0     0     0     0     1     1     0     1
   NaN     0     0     1     1     0     1     1     0     0     0
   NaN     0     1     0     1     0     0     0     1     1     0
   NaN     0     0     0     0     1     1     1     1     0     0
   NaN     0     0     0     0     1     1     0     1     1     0
   NaN     0     1     0     1     1     0     0     0     0     1
   NaN     0     1     1     0     1     0     0     0     1     0
   NaN     1     0     1     0     1     0     1     0     0     0
   NaN     1     0     0     1     1     1     0     0     0     0
   NaN     0     1     0     0     1     0     1     0     1     0
   NaN     1     0     0     0     1     0     0     1     0     1
   NaN     0     1     1     1     0     0     1     0     0     0
   NaN     1     0     0     1     0     1     1     0     0     0
   NaN     0     0     0     1     0     0     0     1     1     1
   NaN     1     1     1     0     0     0     0     0     1     0
   NaN     0     0     1     0     1     1     1     0     0     0
   NaN     0     0     0     0     0     0     1     1     1     1
   NaN     1     1     0     1     0     0     0     0     0     1
   NaN     0     1     0     1     0     0     0     0     1     1
   NaN     0     0     1     0     0     0     1     1     1     0
   NaN     0     0     1     0     1     0     1     1     0     0
   NaN     0     0     1     1     0     0     1     0     1     0
   NaN     0     0     0     1     1     0     0     1     1     0
   NaN     0     0     0     1     1     1     0     0     0     1
   NaN     0     1     1     0     1     0     1     0     0     0
   NaN     0     1     0     0     0     0     1     1     1     0
   NaN     0     1     0     1     0     0     1     0     0     1
   NaN     0     0     1     1     1     0     0     0     0     1
   NaN     1     0     0     0     0     1     1     0     0     1
   NaN     0     1     0     1     0     1     0     0     1     0
   NaN     0     0     1     0     0     1     1     0     0     1
   NaN     1     0     0     0     0     1     1     1     0     0
   NaN     0     0     1     0     0     1     1     0     0     1
   NaN     1     0     0     0     0     0     1     1     0     1
   NaN     0     0     1     1     0     0     0     1     0     1
   NaN     0     0     1     1     0     0     1     1     0     0
   NaN     0     0     1     1     0     1     1     0     0     0
   NaN     0     0     0     0     0     0     1     1     1     1
   NaN     0     1     1     0     0     0     0     1     1     0
   NaN     0     1     1     0     1     0     0     1     0     0
   NaN     1     0     0     1     1     0     0     1     0     0
   NaN     1     1     0     0     1     0     1     0     0     0
   NaN     0     1     1     0     0     1     0     1     0     0
   NaN     1     0     1     0     0     1     0     0     1     0
   NaN     0     0     1     1     0     1     0     1     0     0
   NaN     1     0     0     0     0     1     0     0     1     1
   NaN     0     1     1     0     0     0     1     0     1     0
   NaN     1     1     0     0     0     0     1     0     1     0
   NaN     1     0     0     1     0     1     0     1     0     0
   NaN     0     1     0     0     0     0     1     1     1     0
   NaN     0     1     0     1     1     0     0     1     0     0
   NaN     1     1     0     0     0     0     0     1     0     1
   NaN     1     0     0     0     1     1     0     0     0     1
   NaN     1     0     0     0     0     1     0     1     1     0
   NaN     1     1     0     0     1     0     0     0     0     1
   NaN     0     0     1     0     1     0     0     1     0     1
   NaN     1     1     1     0     0     0     0     0     0     1
   NaN     0     1     0     1     1     0     0     0     0     1
   NaN     0     1     0     0     1     0     1     0     1     0
   NaN     0     1     0     0     0     1     1     1     0     0
   NaN     0     1     1     0     0     1     0     0     1     0
   NaN     1     1     1     0     1     0     0     0     0     0
   NaN     0     0     1     0     1     0     1     0     0     1
   NaN     0     1     0     0     1     1     0     0     1     0
   NaN     1     1     1     0     0     0     1     0     0     0
   NaN     1     0     0     1     0     1     1     0     0     0
   NaN     0     0     1     0     1     1     0     0     1     0
   NaN     0     1     0     1     0     0     0     0     1     1
   NaN     0     1     1     0     0     0     0     0     1     1
   NaN     1     1     0     1     1     0     0     0     0     0
   NaN     1     0     0     0     1     1     1     0     0     0
   NaN     1     1     0     1     1     0     0     0     0     0
   NaN     1     0     0     0     1     0     1     1     0     0
   NaN     0     1     0     1     1     0     0     0     0     1
   NaN     0     1     1     0     0     1     0     0     0     1
   NaN     0     1     0     0     0     1     0     0     1     1
   NaN     0     0     0     1     1     0     0     1     0     1
   NaN     1     1     0     0     1     0     0     1     0     0
   NaN     1     0     0     1     0     0     1     1     0     0
   NaN     0     0     0     1     0     0     1     1     0     1
   NaN     1     0     0     1     1     0     1     0     0     0
   NaN     0     1     0     0     0     0     1     1     1     0
   NaN     0     1     0     0     0     1     1     0     0     1
   NaN     0     1     1     0     1     1     0     0     0     0
   NaN     1     0     0     1     0     1     0     0     0     1
   NaN     0     1     0     0     1     0     0     1     1     0
   NaN     0     0     0     0     1     1     0     0     1     1
   NaN     0     0     1     0     0     1     0     1     0     1
   NaN     1     1     0     0     0     1     1     0     0     0
   NaN     0     0     0     0     0     1     0     1     1     1
   NaN     0     1     1     1     0     1     0     0     0     0
   NaN     0     0     1     0     1     0     1     0     1     0
   NaN     0     0     1     1     1     0     0     1     0     0
   NaN     1     0     0     1     0     0     0     0     1     1
   NaN     0     0     0     1     0     1     0     1     0     1
   NaN     0     0     0     0     0     1     1     1     1     0
   NaN     1     0     1     0     0     0     1     1     0     0
   NaN     1     1     1     0     0     0     1     0     0     0
   NaN     1     1     0     0     0     0     1     1     0     0
   NaN     1     1     0     1     1     0     0     0     0     0
   NaN     0     0     1     1     0     0     0     1     1     0
   NaN     1     1     0     0     0     1     0     0     0     1
   NaN     0     0     0     1     1     0     1     0     0     1
   NaN     0     1     0     0     0     0     1     1     1     0
   NaN     1     0     0     0     0     1     1     1     0     0
   NaN     1     1     1     0     0     0     0     0     0     1
   NaN     0     1     0     0     1     0     0     0     1     1
   NaN     0     1     0     0     0     1     1     1     0     0
   NaN     1     1     0     1     0     0     0     1     0     0
   NaN     0     1     1     0     0     1     0     1     0     0
   NaN     0     1     0     1     1     0     0     0     0     1
   NaN     0     1     0     0     1     0     0     1     0     1
   NaN     0     1     1     0     0     0     1     0     0     1
   NaN     0     0     1     1     0     1     1     0     0     0
   NaN     1     0     1     1     0     0     0     1     0     0
   NaN     0     1     0     0     0     0     1     1     1     0
   NaN     1     0     1     0     1     1     0     0     0     0
   NaN     1     1     0     0     0     0     1     1     0     0
   NaN     0     1     0     1     0     0     1     0     0     1
   NaN     1     1     1     0     1     0     0     0     0     0
   NaN     0     0     1     0     1     0     0     0     1     1
   NaN     1     0     1     0     0     0     1     0     0     1
   NaN     1     0     0     0     0     1     1     0     0     1
   NaN     1     0     0     1     1     1     0     0     0     0
   NaN     1     1     0     1     0     1     0     0     0     0
   NaN     0     1     0     0     1     0     1     0     0     1
   NaN     0     0     1     1     1     0     0     1     0     0
   NaN     0     1     0     0     1     0     1     0     1     0
   NaN     0     1     0     0     1     1     0     1     0     0
   NaN     0     0     0     0     1     1     1     1     0     0
   NaN     0     0     0     0     1     1     1     0     0     1
   NaN     1     0     0     1     0     0     0     1     1     0
   NaN     1     0     0     0     0     1     0     1     0     1
   NaN     0     1     1     1     0     0     0     0     0     1
   NaN     1     1     0     1     0     0     0     0     1     0
   NaN     0     0     0     0     0     1     1     1     0     1
   NaN     1     1     0     0     0     0     0     1     0     1
   NaN     0     0     0     1     1     1     0     1     0     0
   NaN     0     1     0     0     1     1     1     0     0     0
   NaN     0     0     0     1     0     0     1     1     0     1
   NaN     0     0     1     0     0     1     0     1     0     1
   NaN     0     1     0     1     0     0     1     0     0     1
   NaN     0     0     0     0     1     1     0     1     0     1
   NaN     0     0     0     1     1     0     1     1     0     0
   NaN     0     1     0     0     1     0     0     1     1     0
   NaN     1     0     1     0     0     0     1     1     0     0
   NaN     0     0     1     1     1     1     0     0     0     0
   NaN     0     0     1     0     1     0     0     1     1     0
   NaN     0     0     0     1     0     1     1     0     1     0
   NaN     0     0     1     1     0     1     0     0     1     0
   NaN     1     1     0     0     0     0     1     1     0     0
   NaN     1     0     0     1     0     0     0     1     0     1
   NaN     1     0     0     1     1     0     0     0     0     1
   NaN     0     0     0     1     0     1     1     0     1     0
   NaN     0     0     0     1     1     1     0     1     0     0
   NaN     1     1     0     0     1     1     0     0     0     0
   NaN     0     1     0     1     1     0     0     0     0     1
   NaN     0     1     0     1     0     0     0     1     0     1
   NaN     0     1     0     1     0     1     0     0     1     0
   NaN     1     0     0     0     1     1     1     0     0     0
   NaN     1     0     1     1     0     0     0     0     1     0
   NaN     0     1     1     0     0     0     1     0     0     1
   NaN     1     0     1     1     1     0     0     0     0     0
   NaN     1     1     1     0     0     0     0     0     1     0
   NaN     1     0     0     0     0     0     1     1     1     0
   NaN     1     1     0     0     1     0     1     0     0     0
   NaN     0     0     0     1     0     1     1     0     1     0
   NaN     1     0     0     0     0     0     1     1     0     1
   NaN     0     0     0     0     0     0     1     1     1     1
   NaN     1     0     1     0     1     0     1     0     0     0
   NaN     1     1     0     0     0     0     1     1     0     0
   NaN     0     0     1     0     0     1     0     1     1     0
   NaN     0     1     1     0     0     0     1     0     1     0
   NaN     0     1     0     0     1     0     1     0     1     0
   NaN     0     1     0     1     0     0     0     1     0     1
   NaN     0     0     1     1     0     1     0     0     0     1
   NaN     0     1     0     1     1     1     0     0     0     0
   NaN     0     1     0     1     1     0     0     0     1     0
   NaN     1     0     0     0     0     1     1     0     1     0
   NaN     1     1     0     0     0     0     1     0     1     0
   NaN     1     0     0     0     0     0     0     1     1     1
   NaN     1     1     0     0     1     0     0     0     0     1
   NaN     1     1     0     0     0     1     1     0     0     0
   NaN     0     1     1     0     1     1     0     0     0     0
   NaN     0     0     0     0     1     1     0     1     0     1
   NaN     0     0     1     0     0     0     1     0     1     1
   NaN     1     0     0     0     0     1     1     0     0     1
   NaN     1     1     0     0     1     0     0     1     0     0
   NaN     1     0     1     0     0     1     0     0     0     1
   NaN     1     0     0     1     0     0     1     0     1     0
   NaN     1     0     1     0     0     0     0     1     0     1
   NaN     0     1     1     1     1     0     0     0     0     0
   NaN     0     0     0     0     1     1     1     0     1     0
   NaN     0     0     1     0     1     0     0     0     1     1
   NaN     0     1     0     0     1     1     0     1     0     0
   NaN     0     1     0     0     0     0     1     1     1     0
   NaN     1     0     1     0     0     1     0     1     0     0
   NaN     0     1     1     1     0     1     0     0     0     0
   NaN     0     1     0     0     1     1     1     0     0     0
   NaN     0     0     0     0     0     0     1     1     1     1
   NaN     0     0     0     1     1     0     1     1     0     0
   NaN     1     0     0     0     1     0     1     1     0     0
   NaN     0     0     1     0     1     0     0     1     1     0
   NaN     0     0     1     1     0     1     0     0     1     0
   NaN     1     0     1     0     0     1     0     1     0     0
   NaN     0     1     1     0     0     1     1     0     0     0
   NaN     0     0     1     0     1     1     1     0     0     0
   NaN     0     0     0     1     1     1     1     0     0     0
   NaN     0     1     0     1     1     0     1     0     0     0
   NaN     0     0     0     1     0     0     0     1     1     1
   NaN     0     0     1     1     1     0     0     0     1     0
   NaN     0     1     1     0     1     0     0     0     0     1
   NaN     0     1     0     0     0     1     1     1     0     0
   NaN     0     1     0     0     0     1     1     0     1     0
   NaN     1     0     0     0     0     1     1     0     1     0
   NaN     0     1     0     0     0     0     1     1     0     1
   NaN     1     0     0     1     1     0     0     0     1     0
   NaN     0     0     0     1     0     1     1     0     0     1
   NaN     0     1     0     0     0     0     1     0     1     1
   NaN     1     0     1     0     0     1     0     0     0     1
   NaN     0     0     0     0     0     1     1     1     0     1
   NaN     0     1     0     1     0     1     0     0     1     0
   NaN     0     1     0     1     0     1     0     0     0     1
   NaN     1     0     1     0     0     0     0     0     1     1
   NaN     1     0     0     1     0     1     1     0     0     0
   NaN     0     1     0     0     1     0     0     0     1     1
   NaN     1     0     0     0     0     0     0     1     1     1
   NaN     0     1     0     0     0     1     1     1     0     0
   NaN     1     0     0     0     1     0     1     0     0     1
   NaN     1     0     1     0     0     0     0     1     0     1
   NaN     0     1     0     0     1     0     1     0     0     1
   NaN     1     0     0     1     1     1     0     0     0     0
   NaN     1     1     0     0     0     0     1     0     1     0
   NaN     1     0     1     0     0     1     0     1     0     0
   NaN     1     1     1     1     0     0     0     0     0     0
   NaN     1     0     0     1     0     0     0     1     1     0
   NaN     1     1     1     0     0     0     0     0     0     1
   NaN     0     1     1     0     0     1     1     0     0     0
   NaN     1     0     0     1     0     0     0     1     1     0
   NaN     0     1     1     1     0     0     0     1     0     0
   NaN     1     0     1     0     1     0     1     0     0     0
   NaN     1     0     0     0     1     0     0     0     1     1
   NaN     1     0     1     0     1     0     0     1     0     0
   NaN     0     1     0     1     0     0     0     1     0     1
   NaN     0     1     1     0     1     0     0     0     0     1
   NaN     0     1     0     0     0     0     1     1     1     0
   NaN     1     1     1     0     1     1     0     0     0     1
   NaN     1     0     1     0     0     1     1     0     1     1
   NaN     0     1     1     1     0     1     0     0     1     1
   NaN     0     0     0     1     1     1     0     1     1     1
   NaN     1     1     0     1     1     0     0     1     0     1
   NaN     1     1     1     0     1     0     1     1     0     0
   NaN     0     0     1     0     1     1     1     1     1     0
   NaN     0     0     1     0     1     1     1     1     1     0
   NaN     1     0     0     0     0     1     1     1     1     1
   NaN     0     1     0     0     0     1     1     1     1     1
   NaN     1     1     0     1     0     1     0     1     0     1
   NaN     0     1     1     1     0     1     1     0     0     1
   NaN     0     0     0     1     1     1     1     1     1     0
   NaN     1     0     0     0     1     1     1     1     1     0
   NaN     0     1     0     0     1     0     1     1     1     1
   NaN     0     1     1     0     1     1     1     0     1     0
   NaN     1     1     0     1     0     1     1     1     0     0
   NaN     1     1     1     0     0     1     0     1     0     1
   NaN     1     0     0     1     1     0     1     0     1     1
   NaN     1     1     1     0     1     0     0     1     0     1
   NaN     0     1     0     1     0     1     1     0     1     1
   NaN     1     0     0     1     0     1     1     0     1     1
   NaN     0     0     1     0     1     0     1     1     1     1
   NaN     0     0     1     0     1     1     1     1     1     0
   NaN     0     1     0     1     1     0     0     1     1     1
   NaN     1     0     0     0     1     1     0     1     1     1
   NaN     1     0     1     1     0     1     0     0     1     1
   NaN     1     1     1     1     1     0     0     0     1     0
   NaN     1     1     1     0     0     1     1     0     1     0
   NaN     1     1     0     1     0     1     0     0     1     1
   NaN     1     1     1     0     0     1     1     1     0     0
   NaN     1     1     0     0     0     1     0     1     1     1
   NaN     1     0     1     1     0     0     1     1     1     0
   NaN     1     0     1     0     1     1     1     1     0     0
   NaN     0     1     1     0     1     0     1     0     1     1
   NaN     1     0     0     1     1     1     0     1     0     1
   NaN     0     1     0     1     0     1     1     1     1     0
   NaN     0     1     0     0     1     1     1     1     1     0
   NaN     1     1     1     1     0     1     1     0     0     0
   NaN     1     0     0     1     1     1     1     0     0     1
   NaN     1     1     1     1     0     1     0     0     1     0
   NaN     0     1     1     0     1     0     1     0     1     1
   NaN     0     1     1     0     0     1     0     1     1     1
   NaN     1     1     1     0     1     1     0     1     0     0
   NaN     0     0     1     1     1     1     1     0     1     0
   NaN     1     1     0     0     0     0     1     1     1     1
   NaN     1     1     0     0     1     1     1     0     0     1
   NaN     1     0     1     0     0     1     1     0     1     1
   NaN     0     1     1     0     1     0     0     1     1     1
   NaN     0     1     0     1     0     1     0     1     1     1
   NaN     0     1     0     1     0     1     1     0     1     1
   NaN     0     1     1     1     0     1     1     1     0     0
   NaN     0     1     0     1     0     1     1     1     0     1
   NaN     0     0     1     1     0     1     0     1     1     1
   NaN     1     1     1     0     0     0     1     0     1     1
   NaN     0     0     1     1     1     1     0     0     1     1
   NaN     0     0     1     1     0     0     1     1     1     1
   NaN     0     1     1     1     1     1     0     0     1     0
   NaN     1     1     0     1     0     0     1     1     1     0
   NaN     1     1     1     1     0     0     1     0     1     0
   NaN     0     1     1     1     1     0     1     0     1     0
   NaN     0     1     1     1     1     0     1     1     0     0
   NaN     1     0     1     1     1     1     0     1     0     0
   NaN     0     0     1     1     1     1     1     1     0     0
   NaN     0     1     0     0     1     0     1     1     1     1
   NaN     1     1     1     1     1     1     0     0     0     0
   NaN     1     1     0     1     0     1     0     0     1     1
   NaN     1     1     1     0     0     1     0     0     1     1
   NaN     1     0     1     1     1     0     0     0     1     1
   NaN     1     0     0     1     1     1     0     0     1     1
   NaN     0     1     1     0     1     0     1     0     1     1
   NaN     0     1     1     0     1     1     0     1     0     1
   NaN     0     0     1     0     1     1     1     0     1     1
   NaN     1     0     0     0     1     1     1     1     0     1
   NaN     0     0     1     0     1     0     1     1     1     1
   NaN     1     0     0     1     0     1     1     1     0     1
   NaN     0     1     1     1     1     0     0     0     1     1
   NaN     1     0     0     1     1     1     0     0     1     1
   NaN     1     1     1     1     0     0     1     0     1     0
   NaN     1     1     0     1     1     0     1     0     0     1
   NaN     1     0     1     1     0     0     0     1     1     1
   NaN     0     1     1     1     1     0     0     1     0     1
   NaN     1     1     1     0     0     0     1     1     0     1
   NaN     0     0     1     1     0     1     0     1     1     1
   NaN     0     1     0     0     1     1     1     0     1     1
   NaN     1     0     1     0     1     1     1     0     1     0
   NaN     1     1     0     1     1     0     1     1     0     0
   NaN     1     0     1     1     1     0     0     1     1     0
   NaN     1     0     1     0     1     0     1     1     0     1
   NaN     1     0     1     0     0     0     1     1     1     1
   NaN     1     1     0     0     1     0     1     1     1     0
   NaN     1     1     1     0     0     1     1     0     0     1
   NaN     1     1     0     1     0     0     1     1     1     0
   NaN     1     0     1     0     1     0     1     1     1     0
   NaN     0     1     1     1     1     1     0     1     0     0
   NaN     0     0     1     0     1     0     1     1     1     1
   NaN     1     0     0     1     1     0     0     1     1     1
   NaN     1     1     0     1     1     0     0     1     0     1
   NaN     1     0     0     1     0     0     1     1     1     1
   NaN     1     0     1     0     0     1     0     1     1     1
   NaN     1     0     1     0     1     0     1     0     1     1
   NaN     1     0     0     1     1     1     0     1     1     0
   NaN     1     0     1     1     0     1     0     1     1     0
   NaN     1     1     1     0     1     0     1     0     1     0
   NaN     1     0     0     0     0     1     1     1     1     1
   NaN     0     1     1     1     1     0     0     1     1     0
   NaN     0     0     1     1     1     1     1     0     1     0
   NaN     1     1     1     1     0     1     0     0     0     1
   NaN     0     0     1     1     1     1     1     0     1     0
   NaN     0     1     1     1     0     1     1     0     0     1
   NaN     0     1     1     0     1     1     1     1     0     0
   NaN     1     1     0     1     1     0     1     0     0     1
   NaN     0     1     1     1     0     0     1     0     1     1
   NaN     1     0     1     1     0     0     1     1     0     1
   NaN     0     1     1     0     1     0     1     1     1     0
   NaN     1     1     0     0     1     0     1     1     0     1
   NaN     1     1     0     1     1     0     0     0     1     1
   NaN     1     0     1     0     1     1     1     1     0     0
   NaN     0     1     0     1     0     1     1     1     1     0
   NaN     0     0     1     1     1     0     1     0     1     1
   NaN     1     1     1     1     1     0     1     0     0     0
   NaN     0     0     1     1     1     1     1     0     0     1
   NaN     1     0     1     1     0     0     1     1     1     0
   NaN     0     1     0     0     1     1     1     0     1     1
   NaN     0     1     0     0     1     1     1     1     0     1
   NaN     1     1     0     0     0     1     1     1     1     0
   NaN     0     1     0     1     1     0     1     1     0     1
   NaN     0     0     1     1     1     0     1     0     1     1
   NaN     0     1     0     1     1     0     1     1     1     0
   NaN     1     0     0     1     0     1     1     0     1     1
   NaN     1     1     1     0     0     0     1     0     1     1
   NaN     1     0     1     0     1     0     1     0     1     1
   NaN     0     1     1     0     1     1     1     1     0     0
   NaN     1     0     1     1     0     0     1     1     0     1
   NaN     1     0     1     1     0     1     1     1     0     0
   NaN     1     0     0     1     1     0     1     1     1     0
   NaN     1     1     0     0     0     0     1     1     1     1
   NaN     0     0     0     1     1     1     1     0     1     1
   NaN     1     0     0     1     1     0     1     0     1     1
   NaN     0     0     1     1     1     1     1     1     0     0
   NaN     1     0     0     1     1     1     0     1     0     1
   NaN     0     1     1     1     0     0     1     1     0     1
   NaN     1     0     1     1     0     0     0     1     1     1
   NaN     1     1     1     1     1     1     0     0     0     0
   NaN     0     0     0     0     1     1     1     1     1     1
   NaN     1     0     0     0     1     1     1     0     1     1
   NaN     0     0     1     1     1     0     0     1     1     1
   NaN     0     0     1     1     1     0     1     1     0     1
   NaN     1     1     0     1     1     0     0     1     1     0
   NaN     1     0     1     1     1     0     1     0     1     0
   NaN     1     0     0     1     0     1     1     1     1     0
   NaN     0     1     0     1     1     1     0     0     1     1
   NaN     1     1     1     1     1     0     0     1     0     0
   NaN     1     0     0     1     1     0     0     1     1     1
   NaN     1     1     0     1     1     1     0     0     1     0
   NaN     0     1     1     0     1     0     1     0     1     1
   NaN     0     0     1     0     1     1     1     1     1     0
   NaN     1     1     0     0     1     1     0     0     1     1
   NaN     1     0     1     0     1     1     0     1     1     0
   NaN     0     1     0     0     1     1     1     1     0     1
   NaN     0     1     1     1     0     0     1     0     1     1
   NaN     1     0     0     1     1     0     1     1     1     0
   NaN     1     1     1     0     1     0     0     1     0     1
   NaN     1     1     1     0     0     1     1     0     1     0
   NaN     1     0     1     1     1     0     0     1     0     1
   NaN     0     1     1     0     1     1     1     0     0     1
   NaN     0     1     1     1     1     0     0     0     1     1
   NaN     1     1     1     1     1     0     0     0     1     0
   NaN     0     1     0     1     1     0     1     1     1     0
   NaN     1     1     1     1     1     1     0     0     0     0
   NaN     0     1     1     1     1     1     1     0     0     0
   NaN     1     1     1     0     1     0     0     0     1     1
   NaN     1     0     1     0     0     1     1     0     1     1
   NaN     0     1     1     0     1     1     1     1     0     0
   NaN     1     0     1     0     1     0     1     0     1     1
   NaN     1     0     1     0     1     1     0     1     1     0
   NaN     1     0     0     0     0     1     1     1     1     1
   NaN     1     1     1     0     1     1     0     1     0     0
   NaN     1     1     0     1     1     1     0     0     0     1
   NaN     0     0     1     1     1     1     1     1     0     0
   NaN     0     1     1     0     1     0     0     1     1     1
   NaN     1     0     1     0     1     1     0     1     0     1
   NaN     1     0     1     1     0     0     0     1     1     1
   NaN     1     1     1     1     1     0     0     0     1     0
   NaN     1     1     1     0     0     0     0     1     1     1
   NaN     1     1     0     0     1     1     1     0     0     1
   NaN     1     1     0     1     1     0     1     0     1     0
   NaN     1     1     1     0     1     1     0     0     0     1
   NaN     1     1     1     0     1     1     0     0     1     0
   NaN     1     1     0     1     0     1     0     1     1     0
   NaN     1     1     0     0     1     1     0     0     1     1
   NaN     0     0     0     1     1     1     1     0     1     1
   NaN     0     1     1     1     1     0     0     0     1     1
   NaN     1     1     1     0     0     1     1     0     0     1
   NaN     1     0     0     0     1     1     1     1     0     1
   NaN     0     1     1     1     0     0     1     1     1     0
   NaN     0     1     1     0     0     1     1     0     1     1
   NaN     1     1     0     0     1     1     0     1     1     0
   NaN     1     1     0     1     1     0     1     1     0     0
   NaN     0     1     1     0     1     1     0     1     1     0
   NaN     1     1     1     0     1     0     0     0     1     1
   NaN     1     1     0     0     0     1     1     1     0     1
   NaN     1     1     0     1     1     0     0     0     1     1
   NaN     0     0     1     0     1     1     0     1     1     1
   NaN     1     0     1     0     0     1     1     1     1     0
   NaN     1     0     0     1     1     0     0     1     1     1
   NaN     1     0     1     1     0     1     1     1     0     0
   NaN     1     0     1     1     1     1     1     0     0     0
   NaN     0     1     0     0     1     1     1     1     1     0
   NaN     1     0     1     1     1     0     1     0     0     1
   NaN     1     0     0     0     1     1     1     1     0     1
   NaN     1     0     1     1     0     0     0     1     1     1
   NaN     1     0     1     0     1     1     1     0     0     1
   NaN     1     0     0     1     0     0     1     1     1     1
   NaN     0     1     1     0     1     0     1     0     1     1
   NaN     0     1     1     0     1     1     1     0     0     1
   NaN     1     0     1     0     1     0     1     1     1     0
   NaN     1     1     1     0     0     0     1     1     0     1
   NaN     0     1     1     0     0     1     1     1     1     0
   NaN     0     1     1     1     1     0     1     0     0     1
   NaN     1     1     0     0     1     0     1     1     0     1
   NaN     1     0     0     1     1     1     0     0     1     1
   NaN     1     0     1     0     1     0     1     1     0     1
   NaN     1     0     0     1     0     0     1     1     1     1
   NaN     0     0     1     1     1     0     1     0     1     1
   NaN     1     1     1     0     0     1     1     0     1     0
   NaN     1     0     0     0     1     0     1     1     1     1
   NaN     1     1     0     0     0     1     1     1     0     1
   NaN     0     1     0     1     1     1     0     1     0     1
   NaN     1     1     1     1     0     0     0     1     1     0
   NaN     1     0     1     0     1     0     0     1     1     1
   NaN     0     1     1     0     1     1     1     1     0     0
   NaN     1     0     0     1     1     1     0     1     1     0
   NaN     0     1     1     0     1     1     1     0     0     1
   NaN     0     1     1     0     1     0     1     1     0     1
   NaN     1     0     1     1     1     0     0     1     1     0
   NaN     1     1     1     0     0     0     1     1     1     0
   NaN     0     1     0     1     1     1     0     0     1     1
   NaN     0     1     0     0     1     1     0     1     1     1
   NaN     1     1     1     0     0     1     0     1     0     1
   NaN     1     0     0     0     1     1     1     1     0     1
   NaN     0     1     0     1     0     0     1     1     1     1
   NaN     0     1     0     0     1     1     1     1     0     1
   NaN     1     1     0     0     1     1     1     0     1     0
   NaN     0     1     0     1     0     1     0     1     1     1
   NaN     1     1     1     1     1     1     0     0     0     0
   NaN     1     0     1     1     0     0     0     1     1     1
   NaN     1     0     1     1     0     1     0     1     0     1
   NaN     1     1     1     1     0     1     0     1     0     0
   NaN     1     1     1     0     1     1     0     0     0     1
   NaN     1     1     0     0     1     0     1     0     1     1
   NaN     1     0     1     1     1     1     1     0     0     0
   NaN     0     0     1     0     1     1     1     0     1     1
   NaN     1     1     1     0     0     0     0     1     1     1
   NaN     1     0     1     0     1     0     1     0     1     1
   NaN     0     1     0     0     1     1     0     1     1     1
   NaN     0     1     1     0     1     0     1     1     1     0
   NaN     1     1     1     1     0     0     0     0     1     1
   NaN     0     1     0     1     1     1     1     0     1     0
   NaN     0     0     1     0     1     1     0     1     1     1
   NaN     1     1     0     1     1     0     0     1     1     0
   NaN     1     1     0     0     1     1     1     1     0     0
   NaN     1     0     1     0     0     1     1     1     0     1
   NaN     1     0     1     0     0     0     1     1     1     1
   NaN     0     1     1     0     0     1     1     1     1     0
   NaN     1     0     1     0     1     1     0     0     1     1
   NaN     0     1     1     0     1     0     1     1     1     0
   NaN     0     0     0     1     1     1     1     1     0     1
   NaN     1     0     1     1     0     1     1     1     0     0
   NaN     1     1     0     0     1     0     1     1     1     0
   NaN     1     1     1     1     1     0     1     0     0     0
   NaN     0     1     1     1     0     1     1     0     1     0
   NaN     0     1     0     1     1     1     1     0     0     1
   NaN     0     1     1     1     0     0     1     0     1     1
   NaN     1     1     1     1     1     0     0     0     0     1
   NaN     1     1     0     0     0     1     1     1     0     1
   NaN     0     1     0     1     1     0     1     0     1     1
   NaN     0     1     0     1     0     1     1     1     1     0
   NaN     0     0     1     0     1     1     1     0     1     1
   NaN     0     0     1     0     0     1     1     1     1     1
   NaN     0     1     0     1     1     0     0     1     1     1
   NaN     1     0     0     1     1     0     1     1     0     1
   NaN     0     1     1     1     1     1     0     0     1     0
   NaN     1     1     1     0     1     0     0     1     0     1
   NaN     1     0     0     1     0     1     1     1     1     0
   NaN     1     1     0     1     0     0     1     0     1     1
   NaN     1     1     0     0     1     1     0     1     1     0
   NaN     1     1     0     0     1     1     1     1     0     0
   NaN     1     1     0     0     1     1     1     0     1     0
   NaN     1     0     1     0     1     0     1     0     1     1
   NaN     1     0     1     1     0     1     1     0     0     1
   NaN     1     1     0     1     0     0     0     1     1     1
   NaN     1     0     1     0     1     1     0     1     1     0
   NaN     0     0     1     1     0     1     0     1     1     1
   NaN     0     1     0     0     1     1     1     1     0     1
   NaN     1     0     0     0     1     1     0     1     1     1
   NaN     1     0     0     1     1     0     1     1     1     0
   NaN     1     1     1     0     1     1     1     0     0     0
   NaN     1     1     1     1     0     0     1     0     1     0
   NaN     0     1     1     1     0     1     0     0     1     1
   ];

num_stim = size(stimuli,1);


%List of pre-configured parameters

%1 if paranmeters derived from empirical study
%2 if small test grid search
%3 if large grid search
% if config_list == 1;

%cols: prior alpha (learning rate), beta (inverse temperature)
parameters = [
                      0.17                       0.1                       7.4
                       0.5                         0                       150
                       0.2                         0                      1.43
                      0.48                      0.03                     19.68
                      0.14                         0                      1.02
                      0.91                      0.89                         1
                      0.32                      0.15                      4.51
                      0.12                      0.27                      2.99
                      0.47                         0                      1.09
                      0.43                      0.49                      3.64
                       0.5                         0                       150
                       0.5                         0                       150
                       0.5                      0.03                     34.28
                       0.5                         1                      4.41
                      0.49                      0.01                       150
                      0.56                         0                      1.27
                       0.5                      0.01                    133.24
                      0.49                         0                       150
                      0.57                      0.78                      1.13
                      0.23                      0.52                         1
                       0.5                         0                      1.09
                      0.82                         0                      1.13
                      0.43                      0.08                      4.95
                       0.5                         0                       150
                      0.52                      0.07                         5
                       0.5                         0                      40.8
                      0.44                         0                     15.87
                      0.45                         0                      1.59
                      0.58                      0.14                      4.34
                      0.33                      0.07                      6.59
                      0.48                      0.05                     13.13
                       0.5                         0                       150
                         0                      0.12                      3.21
                      0.52                      0.77                      1.86
                      0.37                      0.08                      5.39
                      0.46                      0.98                      5.05
                      0.76                      0.28                      2.12
                       0.5                         0                       150
                         0                      0.15                      5.05
                      0.37                         0                     16.58
                      0.64                      0.34                         1
                         0                         0                      1.16
                       0.5                         0                      1.15
                      0.59                      0.12                      1.54
                       0.5                         0                       150
                      0.43                      0.03                      5.77
                       0.5                         0                      1.03
                      0.34                      0.26                      2.13
                       0.5                         0                       150
                      0.48                      0.14                      6.15
                      0.42                         0                     14.47
                       0.2                      0.04                      5.49
                      0.01                      0.76                         1
                      0.57                      0.71                      1.45
                       0.5                         0                       150
                      0.46                      0.01                         1
                       0.5                         0                       150
                      0.35                      0.09                      1.95
                      0.46                      0.05                         1
                      0.51                         0                       150
                      0.98                         0                      1.82
                      0.45                      0.06                      3.35
                      0.17                      0.04                      3.62
                      0.71                         1                      2.15
                       0.8                      0.15                         1
                      0.48                      0.04                        12
                       0.5                         0                       150
                      0.08                      0.44                     46.68
                       0.5                         0                      1.76
                      0.08                      0.23                      1.12
                         1                         1                         1
                      0.32                      0.14                      3.29
                      0.06                      0.23                      1.16
                       0.5                         0                       1.5
                      0.44                      0.19                      4.21
                      0.38                      0.05                     10.05
                      0.18                         0                      2.18
                       0.5                         0                       150
                      0.48                      0.01                         1
                      0.52                      0.09                      3.47
                       0.5                         0                      1.13
                      0.56                      0.11                      3.12
                       0.5                         0                       150
                       0.5                         0                       150
                      0.52                      0.23                      2.31
                      0.06                         0                      1.01
                       0.5                         0                    149.93
                       0.5                         0                       150
                       0.5                      0.26                         1
                       0.5                         0                       150
                       0.5                         0                     60.69
                      0.57                      0.12                      3.58
                       0.5                         0                       150
                      0.35                      0.13                      1.42
                       0.5                         0                       150
                      0.48                      0.45                      1.14
                      0.69                         0                         1
                       0.5                         0                       150
                       0.5                         0                       150
                      0.19                      0.02                      1.02
                      0.73                      0.72                      1.28
                      0.51                      0.02                     28.18
                      0.65                         0                      4.56
                       0.5                      0.02                         1
                       0.5                         0                       150
                       0.5                         0                     145.5
                       0.5                         0                      1.06
                      0.55                      0.13                      2.26
                       0.5                         0                      1.14
                      0.41                      0.02                      4.12
                      0.12                      0.53                         1
                      0.23                      0.25                      1.32
                       0.5                         0                      1.01
                      0.53                      0.02                     35.58
                      0.48                         0                       3.1
                      0.59                       0.1                         1
                      0.51                       0.1                       3.8
                      0.33                      0.07                      9.46
                         0                      0.86                      1.56
                       0.5                         0                       150
                      0.31                      0.06                      2.31
                       0.5                         0                       150
                      0.42                      0.05                         1
                      0.45                         0                      1.84
                      0.26                         0                      2.12
                       0.5                         0                      1.02
                      0.02                      0.92                         1
                      0.49                         0                       150
                      0.54                      0.73                         1
                      0.42                         0                     24.58
                      0.45                         0                         1
                       0.5                         0                       150
                      0.49                         0                     37.89
                       0.5                         0                       150
                       0.5                         0                       150
                       0.5                         0                       150
                       0.5                      0.01                       150
                      0.28                         1                      1.21
                       0.5                         0                     18.32
                         1                      0.33                       1.5
                       0.5                         0                         1
                       0.5                         0                       150
                      0.45                         0                      1.65
                      0.65                      0.15                      2.29
                       0.5                         0                       150
                       0.5                         0                       150
                      0.02                         0                      2.28
                      0.41                         0                     25.47
                       0.5                         0                       150
                      0.39                      0.21                      1.73
                      0.36                      0.39                      2.04
                       0.5                         0                      1.66
                       0.5                         0                      1.22
                       0.5                         0                       150
                         0                         0                      2.88
                       0.3                         0                       2.1
                      0.14                       0.3                      2.91
                       0.7                         0                      1.02
                      0.46                      0.01                         1
                      0.57                         0                      1.08
                       0.5                         0                      1.01
                       0.5                         0                       150
                      0.57                         0                      1.27
                       0.5                         0                    145.71
                       0.5                         0                     144.5
                      0.54                         1                      2.06
                      0.43                      0.04                         1
                      0.96                      0.27                         1
                       0.5                      0.01                       150
                      0.42                      0.81                      1.11
                       0.5                         0                       150
                      0.05                         0                      2.82
                      0.16                         0                      1.13
                       0.5                      0.36                         1
                      0.12                      0.28                       1.8
                      0.58                         0                      1.02
                      0.49                      0.01                       150
                      0.49                      0.08                       6.7
                       0.5                         0                       150
                       0.5                         0                         1
                      0.56                      0.16                      5.84
                      0.69                      0.51                         1
                       0.5                         0                      1.05
                       0.5                      0.01                       150
                       0.5                         0                       150
                      0.54                         0                      5.44
                       0.5                         0                       150
                      0.47                      0.02                      11.2
                      0.53                      0.06                      1.17
                      0.49                         0                       150
                         0                      0.07                      3.56
                       0.5                         0                       150
                       0.5                         0                       150
                      0.59                       0.2                      3.13
                      0.49                      0.01                       150
                      0.54                      0.04                      3.26
                       0.5                         0                       150
                      0.46                      0.13                      7.56
                       0.5                         0                       150
                      0.49                      0.47                     31.98
                      0.37                         0                      2.27
                      0.57                         1                      5.11
                       0.5                         0                      1.53
                         1                      0.23                      1.37
                       0.5                         0                       150
                      0.23                         0                      2.01
                       0.5                      0.01                       150
                       0.5                         0                    148.13
                      0.55                      0.12                      2.43
                      0.91                      0.18                         1
                       0.5                         0                       150
                      0.57                         0                         1
                       0.5                         0                       150
                      0.51                       0.1                      2.44
                       0.5                         0                     125.2
                      0.49                         0                       150
                       0.5                         0                       150
                      0.24                      0.12                      2.14
                      0.32                         0                      2.79
                       0.5                         0                       150
                      0.43                      0.76                      1.54
                      0.49                      0.01                     31.84
                      0.69                      0.14                         1
                      0.38                         0                      1.34
                      0.71                      0.35                         1
                      0.51                      0.09                      1.34
                      0.85                         1                       4.6
                       0.5                         0                      1.01
                      0.47                      0.03                         1
                      0.58                         0                      1.03
                      0.48                      0.06                         1
                      0.48                      0.01                     16.28
                      0.49                      0.01                       150
                      0.51                      0.02                     18.88
                      0.49                         0                       150
                      0.35                      0.13                      1.63
                      0.64                         0                       1.5
                         0                      0.43                         1
                       0.5                         1                       1.1
                      0.77                      0.24                      1.01
                      0.69                         0                      1.27
                       0.5                         0                      1.73
                      0.48                         0                     126.5
                      0.56                      0.12                      6.21
                      0.51                      0.41                         1
                       0.5                         0                       150
                       0.4                      0.19                      2.95
                      0.43                         0                     34.34
                      0.82                      0.87                         1
                       0.5                         0                      1.03
                       0.5                         0                       150
                       0.5                         0                         1
                       0.5                         0                       150
                      0.55                      0.06                         1
                       0.5                         0                       150
                      0.63                      0.58                         1
                      0.73                      0.19                         1
                         0                      0.38                         5
                      0.49                         0                       150
                      0.41                      0.13                      4.68
                       0.5                         0                       150
                      0.94                      0.63                      1.39
                       0.5                         0                    149.48
                       0.5                         0                         1
                       0.5                         0                       150
                      0.52                       0.1                      2.92
                      0.01                         0                      1.02
                      0.49                         0                       150
                      0.54                         0                      1.06
                       0.5                         0                       150
                      0.53                      0.61                      2.77
                       0.5                         1                      5.11
                       0.5                         0                       150
                      0.46                      0.12                     11.03
                      0.53                      0.56                         1
                       0.5                         0                      1.01
                         1                         1                      1.39
                      0.43                      0.64                         1
                      0.51                      0.02                      7.45
                       0.5                         0                       150
                      0.52                         0                      6.49
                       0.6                      0.17                      1.63
                       0.5                         0                       150
                       0.5                      0.89                      1.27
                      0.68                       0.1                         1
                         0                         0                      1.38
                       0.5                         0                      1.03
                      0.56                       0.5                         1
                      0.71                         0                         1
                      0.52                         0                         1
                      0.52                      0.56                      1.13
                      0.39                      0.05                         1
                       0.5                         0                       150
                       0.5                      0.05                         1
                       0.4                      0.09                      3.42
                       0.5                         0                       150
                      0.47                      0.18                      1.95
                      0.55                      0.78                         1
                      0.46                      0.06                      1.59
                      0.76                         0                      1.13
                      0.41                      0.03                      2.12
                      0.12                      0.24                         1
                      0.44                         0                      7.27
                      0.28                         0                      1.56
                       0.4                      0.03                         1
                       0.5                         0                       150
                       0.5                         0                       150
                       0.4                         0                     14.55
                       0.5                         0                      1.15
                      0.14                      0.76                      1.34
                       0.5                         0                       150
                       0.5                         0                       150
                         0                      0.14                      1.33
                      0.58                      0.05                         1
                       0.5                         0                      1.02
                      0.64                      0.98                      1.48
                      0.49                         0                       1.7
                       0.4                         0                      1.04
                      0.27                      0.09                       2.3
                       0.5                         0                       150
                      0.02                         0                      1.44
                         0                         0                      1.09
                      0.64                      0.19                      3.32
                      0.71                      0.21                         1
                      0.47                      0.04                      6.67
                      0.36                      0.88                      2.26
                      0.57                      0.05                         1
                       0.5                         0                       150
                       0.5                         0                       150
                      0.46                         0                      1.11
                       0.5                         0                      1.01
                       0.5                         0                       150
                      0.14                      0.07                      3.08
                      0.65                       0.3                         1
                       0.8                      0.35                      2.21
                       0.5                         0                       150
                       0.5                         0                      1.75
                      0.68                         1                      1.13
                      0.26                       0.1                         1
                      0.22                         0                      3.26
                      0.77                      0.35                         1
                         0                      0.68                         1
                      0.52                         0                      1.59
                      0.97                         0                      1.48
                      0.37                      0.12                      2.26
                      0.88                         0                         1
                       0.1                      0.97                      5.77
                       0.5                         0                       150
                      0.04                         0                       2.1
                      0.45                      0.16                      1.93
                         1                         0                      1.44
                      0.71                       0.1                         1
                      0.69                      0.18                         1
                         0                      0.09                      1.57
                       0.5                         0                         1
                       0.5                         0                      1.03
                       0.5                         0                       150
                      0.11                         0                      5.95
                      0.28                      0.09                         1
                      0.49                      0.05                      7.17
                      0.47                      0.15                      4.02
                       0.6                      0.01                         2
                      0.51                      0.03                         1
                       0.5                         0                      1.79
                       0.5                      0.03                      5.74
                       0.5                         0                       150
                      0.69                      0.02                      1.15
                      0.55                       0.1                         1
                      0.26                         1                         1
                       0.3                      0.06                      1.69
                       0.6                      0.68                         1
                      0.49                      0.02                     24.07
                      0.28                       0.1                      2.01
                      0.63                      0.96                      3.89
                      0.55                      0.19                      1.29
                       0.5                         1                      1.45
                       0.5                         0                    144.31
                       0.5                         0                       150
                      0.67                      0.08                         1
                         0                      0.05                      3.34
                      0.41                      0.03                      5.62
                       0.5                         0                         1
                       0.5                         0                       150
                       0.5                         0                      1.27
                       0.5                         0                      1.02
                      0.56                         0                      1.09
                       0.5                         0                       150
                      0.19                      0.24                         1
                       0.5                         0                      1.19
                      0.32                         0                      2.31
                      0.43                         0                      1.17
                       0.5                         0                      1.03
                       0.5                         0                       150
                         0                      0.12                      1.35
                      0.06                      0.74                      1.04
                       0.5                         0                    100.65
                       0.5                         0                       150
                       0.3                         0                      7.03
                       0.7                      0.11                      1.28
                      0.41                       0.1                      3.67
                      0.27                       0.1                         1
                      0.58                      0.75                      1.35
                      0.47                      0.11                      3.62
                       0.5                         0                       150
                      0.29                      0.29                         1
                      0.53                      0.01                         1
                      0.94                         0                      1.01
                      0.47                      0.03                     11.82
                      0.49                      0.02                     13.11
                       0.4                         0                      1.21
                      0.49                         0                       150
                       0.5                         1                      2.08
                      0.57                      0.47                      1.05
                       0.7                      0.61                         1
                       0.5                         0                       150
                       0.4                      0.01                      3.28
                       0.5                         0                      1.36
                      0.32                         0                      1.89
                       0.5                         0                    126.28
                      0.61                      0.27                         1
                         1                         1                      1.03
                      0.07                      0.08                      1.61
                      0.66                      0.36                         1
                      0.03                         0                      2.73
                      0.16                         0                      1.94
                      0.76                      0.28                      1.44
                       0.5                         0                       150
                      0.45                         0                      6.76
                       0.5                         0                       150
                      0.68                      0.81                         1
                      0.59                      0.17                      2.06
                      0.96                         0                      1.48
                      0.38                         0                         1
                      0.48                         1                         1
                       0.5                         0                      1.27
                      0.48                      0.23                         1
                      0.45                      0.03                     16.03
                      0.55                         0                      1.01
                       0.5                         0                    149.25
                      0.03                         0                      4.88
                      0.71                      0.23                         1
                      0.38                      0.31                         1
                       0.7                         0                      1.62
                       0.5                         0                       150
                      0.59                      0.97                      5.67
                       0.5                         0                       150
                       0.5                         0                       150
                      0.74                         0                      1.03
                       0.5                         0                      1.19
                       0.5                         0                       150
                      0.51                         0                     41.78
                       0.5                         0                       150
                       0.5                         0                       150
                      0.18                         0                      1.24
                       0.5                         0                      1.31
                       0.5                         0                     99.82
                       0.5                         0                       150
                       0.5                         0                    110.23
                       0.5                         0                       150
                      0.61                         0                      1.01
                      0.39                      0.18                      1.38
                      0.58                      0.06                      4.14
                      0.05                      0.65                         1
                      0.75                         0                         1
                      0.51                      0.06                      4.74
                         0                       0.9                      1.12
                       0.5                         0                       150
                       0.5                         0                       150
                       0.5                         0                       150
                       0.5                      0.01                     12.73
                      0.52                      0.24                         1
                       0.5                         0                       150
                       0.5                         0                       150
                       0.5                         0                      1.21
                       0.5                         0                      1.34
                      0.18                      0.11                      1.71
                      0.21                         0                      5.37
                      0.27                      0.01                      4.98
                      0.51                      0.04                     13.48
                      0.48                      0.41                         1
                      0.02                         0                      1.79
                      0.15                      0.31                         1
                      0.09                         0                      1.22
                      0.33                      0.14                         1
                       0.5                         0                       150
                       0.5                      0.22                         1
                      0.45                         0                      1.53
                      0.51                         0                         1
                      0.64                      0.46                         1
                      0.09                         0                      3.78
                       0.5                         0                         2
                       0.5                         0                       150
                      0.26                         0                       1.6
                      0.25                      0.93                         1
                       0.5                         0                       150
                      0.48                      0.82                         1
                         1                      0.77                         1
                      0.05                         0                      1.26
                       0.5                         0                       150
                       0.5                         0                      1.02
                      0.48                         0                     42.03
                      0.12                      0.14                     87.97
                      0.59                      0.17                      1.74
                      0.66                      0.04                         1
                      0.88                         0                       1.3
                      0.64                      0.16                      1.98
                      0.49                         0                       150
                       0.5                      0.31                      1.99
                       0.5                         0                       150
                       0.5                         0                    148.15
                       0.4                      0.07                      4.33
                      0.54                      0.21                         1
                      0.83                      0.39                         1
                      0.77                         0                      1.03
                      0.99                         0                      2.47
                       0.5                         0                      1.04
                      0.68                         0                       1.1
                      0.38                      0.09                      3.24
                       0.5                      0.13                      6.53
                      0.13                         0                      2.34
                      0.49                      0.27                         1
                       0.5                         0                       150
                      0.85                         0                         1
                       0.5                         0                       129
                      0.63                      0.01                      2.33
                      0.53                      0.07                      3.14
                       0.6                      0.07                      1.16
                      0.51                       0.6                      3.04
                      0.62                      0.27                         1
                       0.5                         0                       150
                      0.72                         1                         1
                      0.46                       0.2                      1.82
                      0.27                      0.08                      2.94
                      0.49                         0                     50.98
                      0.51                      0.11                         1
                      0.54                         0                      1.31
                      0.49                         0                      1.01
                      0.52                      0.02                         1
                       0.5                         0                      1.13
                      0.64                      0.25                         1
                      0.06                         1                      4.41
                         1                      0.27                     30.85
                      0.49                      0.41                         1
                      0.08                       0.9                         1
                      0.25                         0                       1.7
                       0.5                         0                      1.05
                       0.5                         0                       150
                         0                      0.54                         1
                         0                      0.36                      1.66
                       0.5                         0                       150
                      0.55                      0.06                      3.48
                       0.5                         0                       150
                         0                      0.01                      1.27
                      0.48                      0.88                         1
                       0.2                      0.14                      3.19
                       0.5                         0                      1.03
                      0.02                       0.2                      1.77
                      0.41                      0.82                      1.57
                       0.5                         0                       150
                       0.5                      0.02                     15.57
                      0.49                      0.07                         1
                      0.26                         0                      1.31
                       0.1                      0.22                         1
                      0.36                         0                      2.94
                      0.53                      0.92                       6.2
                       0.5                         0                       150
                      0.51                         0                       150
                      0.32                      0.02                         1
                      0.39                         0                      5.44
                      0.46                      0.03                         1
                       0.5                         0                       150
                      0.49                         0                      89.7
                      0.65                      0.92                      1.29
                         0                      0.11                     16.02
                      0.41                         0                      4.75
                      0.56                      0.03                      3.23
                      0.54                         0                      1.18
                      0.56                      0.63                         1
                      0.59                      0.21                         1
                      0.82                      0.34                      1.02
                       0.5                         0                       150
                       0.5                      0.01                     41.46
                       0.5                         0                      96.9
                       0.5                         0                    147.53
                       0.5                         1                       1.1
                      0.37                      0.05                       5.1
                       0.5                         0                       150
                      0.34                         0                      1.19
                       0.5                      0.03                      5.38
                       0.5                         0                       150
                      0.49                      0.32                         1
                       0.5                         0                    118.08
                       0.5                         0                       150
                      0.53                      0.11                         1
                      0.42                      0.92                      1.41
                       0.5                         0                       150
                       0.5                         0                       150
                       0.5                         0                       150
                      0.79                      0.89                         1
    ];

%initial values and bounds: prior, alpha, beta
params = [0.5 0.5 12];
lower_bounds = [0 0.1 0];   %fitting will not try parameters outside these values
upper_bounds = [1 1 150];

num_combos = size(parameters,1);
for combo = 1:num_combos;

    this_sequence = stimuli(combo,:);
    these_params = parameters(combo,:);

    %simulate responses to this sequence of agent using this combo of parameters
    configured_probabilities = prob_guilt_delta(these_params,this_sequence);

    %Now fit the "participant" (combo) data you just created

    clear data_to_fit;
    data_to_fit.stimuli = this_sequence;
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
ticklabels = {'Prior', 'Alpha', 'Beta'};
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

    this_sequence = data_to_fit.stimuli(stimulus,:);

    %simulate responses to this sequence of agent using this combo of parameters
    fitted_probabilities(stimulus,:) = prob_guilt_delta(params,this_sequence);

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
        %          std_resid = 6.6/100; %in forensic_beads_study1_2024_v3.m I took all the fitted probabilities, subtracted from them the human probabilities, took the abs of these differences, then found the std over all trials in all participants. Then I divide by 100, as that program was in percentages and here we have proportions.
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



