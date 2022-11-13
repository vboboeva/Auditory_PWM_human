% these are the all the parameters expected to run PoissonEvents
%
% for each subject Myname, make another file Myname/Myname_startup.m that calls this
% file, then make any appropriate changes (such as subject contact,
% modality, etc.)

subject_contact = '';
command_to_run = 'PWMHuman';
modality = 'auditory'; % 'auditory' or 'visual' or 'colors' or 'orientation'
fname_string = 'testing';

commit_to_svn = 0;
send_automated_email = 0; % email session summary to experimenters and subject
send_session_analysis = 0; % email analysis of session to experimenters
experimenter_contacts = {'aakrami@princeton.edu'};

feedback_interval = 10; % feedback of mean accuracy every this many trials
display_session_summary = 1;
min_paid_session_length = 200;
a_payscale = [0         0.70      0.725     0.75      0.775     0.80      0.825     0.85      0.9; ... % accuracy threshold
              0.029     0.03      0.031     0.032     0.036     0.044     0.048     0.051     0.053]; % dollars per correct answer
o_payscale = [0         0.70      0.725     0.75      0.775     0.80      0.825     0.85      0.9; ... % accuracy threshold
              0.029     0.03      0.031     0.032     0.036     0.044     0.048     0.051     0.053]; % dollars per correct answer
          
          
fixed_fixation_dur = 2;
punish_early_responses = 1;

show_fixation_in_prepause = 1;

a_volume = 0.5;
c_contrast = 0.8; % should vary between [0.5 1], 1 being the most contrast, 0.5 being they're identical


