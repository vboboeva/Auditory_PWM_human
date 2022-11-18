%% to produce stimuli for each session. only one modality

%% General session matrix parameters: INSERT YOUR VALUES IN THIS CELL!!!
%
%clear all; close all; clc;
%

data_dir = '/home/vizhe/Documents/Auditory_PWM_human/data/'; %ALL EXPERIMENTS subject folders will be saved here!
% expNum_folderName = name; %the subfolder you wish to create! All the subject folders and trial paramter files will be stored in this folder. 

sdi = .15;  % the formula for sdi =(sd1 - sd2)/(sd1 + sd2)
ratio = (1+sdi)/(1-sdi);
ratio1 =2.5;
alfa = (ratio1-1)/(1+ratio1);
R = 10*log10(ratio);
R1 = 10*log10(ratio1);
sd1 = 1; %define the highest value stimulus in the stimulus set
stim_levels = 6; % the number of base&comparison stimuli levels; if diamond=0, comparison stimuli levels will be base_stim_count+2;
%Number of points in SGM = stim_levels*2 - 2;

diamond = 1; %1=yes, 0=no: specify if you want a diamond shaped SGM or 'vertical-lines' SGM

%Number of trials: 
%Trials per session = (stim_levs*2)-2 X repsPerSes;
%Sessions in experiment = modality count X session count;
num_trials = 400;
repsPerSes = num_trials / (stim_levels - 1) / 2; %in every session each SGM data point(stim pair) is presented this many times
modalities = 1; %the modalities vector; 1:A-A; 2:T-T; 3:TA-TA; 4:A-T; 5:T-A;
session_count = 5; %sessions per modality; session modality picked randomly from modalities vector without replacement; when all chosen, restart with a full set;

split_exp = 1; % 1=don't split the expriment in 2 visits; or 2=split in half; if 2, the response rule for the second half of the sessions will get be the opposite; 
seedSet = 50; %how many pseudorandom noise patterns to use;
stim1_dur = 400; %ms, FIXED VALUE;
stim2_dur = 400; %ms, FIXED VALUE;
pre_delay = 250; %ms, FIXED VALUE; to make it variable, implement the code as used for post_delay
inter_delay = [2000 4000 6000]; %ms, FIXED VALUE; [500,200,100]; %
post_delay = [200 300]; %ms, the length of this vector must divide number of points in SGM * repsPerSes
%% Consistency check:
balanced_delay = num_trials/length(post_delay);
if mod(balanced_delay,1) 
   disp('Nothing was generated! Element count in "post delay" must divide "num_trials" without remainder. Add/remove values from "post_delay"!')
return
end
%% Generate stimuli values, show a plot of the SGM
%Stim values:
if diamond
    sd1 = (sd1+sd1*sdi)/(1-sdi);
end
stim_lev_count = stim_levels + 2; 
stimVek = zeros(stim_lev_count,1);
stimVek(stim_lev_count) = sd1;
for pp = stim_lev_count-1:-1:1
    sd2 = (sd1-sd1*sdi) / (sdi + 1);
    sd1 = sd2;
    stimVek(pp) = sd2;
end

ndx_lst = 2:stim_lev_count-1; %drop the first and the last value, they are for the stim_lev_2

%Create all the stimuli pairs (if you want diamond shaped SGM, exclude
%first and last points of stim_levs)
for rr = 1:size(ndx_lst,2)
    stim_levs(rr*2-1:rr*2,1) = [stimVek(ndx_lst(rr)), stimVek(ndx_lst(rr))];
    stim_levs(rr*2-1,2) = stimVek(ndx_lst(rr)-1);
    stim_levs(rr*2,2) = stimVek(ndx_lst(rr)+1);
end
if diamond
    stim_levs = stim_levs(2:end-1,:);
end
%figure
%loglog(stim_levs(:,1),stim_levs(:,2),'o')
%grid on

% number of pairs of stimuli
num_stim_pairs = size(stim_levs,1);
num_pair_pairs = num_stim_pairs/2;

switch distr_type
case 'Uniform'
    reweight = ones(num_pair_pairs);
case 'NegSkewed'
    p_ratio = 5.;
    lambda = log(p_ratio)/(num_pair_pairs - 1);
    reweight = exp(-lambda*[1:num_pair_pairs]);
    reweight = flip(reweight);
case 'PosSkewed'
    p_ratio = 5.;
    lambda = log(p_ratio)/(num_pair_pairs - 1);
    reweight = exp(-lambda*[1:num_pair_pairs]);
case 'Bimodal'
    lambda = 2; % increase this for "sharper" bimodality
    reweight = exp(-lambda*[1:num_pair_pairs]);
    reweight = reweight + flip(reweight);
otherwise
    message = sprintf('Invalid dist_type = %s: switching to uniform', distr_type);
    disp(message);
    reweight = ones(num_pair_pairs);
end

% generate data with weights from chosen distribution
% Replicate stimulus values to fill the matrix pseudorandomly

counts_pair_pairs = round(reweight/sum(reweight) * repsPerSes * num_pair_pairs);
replic_stims = [];
for pp = 1:num_pair_pairs
    replic_stims = [replic_stims; repmat(stim_levs(pp*2-1,:), counts_pair_pairs(pp), 1)];
    replic_stims = [replic_stims; repmat(stim_levs(pp*2,:), counts_pair_pairs(pp), 1)];
end
    
stim1_col = replic_stims(:,1);
stim2_col = replic_stims(:,2);

% save('data.mat')

% return

rep_count = size(replic_stims,1)/size(post_delay, 2);
replic_post_delay = repmat(post_delay',rep_count,1);
replic_inter_delay= repmat(inter_delay',rep_count, 1);

%% Set up variables before loop. In the loop they only get reshuffled
modality_count = size(modalities,2);
tot_ses = modality_count * session_count;
trial_count = size(replic_stims,1);
%%trl_mtx = zeros(trial_count,29);  %empty block matrix; columns count is fixed

fullSeedSet_times = floor(trial_count/seedSet);
fullSeedSet_rep = repmat((1:seedSet),1, fullSeedSet_times);
remainSeeds_count = rem(trial_count,seedSet);
remainSeeds_rand = randperm(seedSet);
concat_allSeeds = [fullSeedSet_rep remainSeeds_rand(1:remainSeeds_count)];

cond_distr = reshape((1:tot_ses),modality_count, session_count)';
cond_shuff = zeros(size(cond_distr));
for kk = 1:size(cond_shuff,1)
    modIndx = randperm(modality_count);
    rw = cond_distr(kk,:); %shuffle order in which modalities come
    cond_shuff(kk,:) = rw(modIndx); %maintain pseudorandom structure: each modality once, only then the next round!;
end
%% Generate complete folder tree with all session parameters for all subjects

curr_folderNum = cond_shuff;
%trl_mtx = cell(16,1);  %start with empty matrix for each session
trl_mtx = struct;
trl_mtx.folderNum = curr_folderNum; %session
trl_mtx.trialnum = (1:trial_count); %trial
trl_mtx.modalities= modalities; %modality
trl_mtx.a1_time = ones(trial_count,1)*stim1_dur; %stim1 and stim2 duration

stimRandIndx = randperm(trial_count);
trl_mtx.a1_sigma= stim1_col(stimRandIndx); %stim1 values

randSeedIndx1 = randperm(trial_count);
concat_allSeeds1 =concat_allSeeds(randSeedIndx1);
trl_mtx.a1_seed = concat_allSeeds1; %seed for stim1

trl_mtx.a1_seedh = concat_allSeeds1+1; %high freq seed for stim1
trl_mtx.a2_time = ones(trial_count,1)*stim2_dur; %stim2 duration
trl_mtx.a2_sigma = stim2_col(stimRandIndx); %stim2 (using same randomized index so pairs are preserved)

randSeedIndx2 = randperm(size(concat_allSeeds,2));
concat_allSeeds2 =concat_allSeeds(randSeedIndx2);

trl_mtx.a2_seed = concat_allSeeds2; %seed for stim2

trl_mtx.a2_seedh = concat_allSeeds2+1; %high freq seed for stim2
trl_mtx.prestim= ones(trial_count,1)*pre_delay; %pre stim delay
  
randDelayIndx = randperm(trial_count);
trl_mtx.delaydur= replic_inter_delay(randDelayIndx);%delay between the two stimuli

randDelayIndx = randperm(trial_count);
trl_mtx.poststim= replic_post_delay(randDelayIndx); %variable post delay


comp_stims = (stim2_col(stimRandIndx) < stim1_col(stimRandIndx));
trl_mtx.siderule= comp_stims;
trl_mtx.trialdur = (pre_delay+stim1_dur + stim2_dur + replic_post_delay(randDelayIndx) +replic_inter_delay(randDelayIndx));
trl_mtx.sessiondur = (3000 + sum(pre_delay+stim1_dur + stim2_dur + replic_post_delay(randDelayIndx) +replic_inter_delay(randDelayIndx)))./60000;

%Prepare to save
if ~exist (data_dir,'dir')
    mkdir(data_dir);
else
%     clear all;
%     disp('Folder already exists! Nothing was generated!')
%     return
end
save([data_dir,'trl_mtx.mat'],'trl_mtx');

%save('data.mat')

disp('Done building session matrices!')