clear all
name = 'Test';
Default_startup_pwm;

command_to_run = 'PWMHuman';
min_paid_session_length = 300;
show_fixation_in_prepause = 0;

modality = 'auditory';
filter_type ='GAUS';
fcut = 110;
lfreq = 2000;
hfreq = 6000;

% pick a type of distribution (uniform is Athena's original code)
% distr_type = 'Uniform'
% distr_type = 'Bimodal';
distr_type = 'NegSkewed';
% distr_type = 'PosSkewed'

PWMHuman_gen_trial_mtx;

fname_string = ['Sub_' name '_' distr_type]
