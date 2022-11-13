name = 'Athena';

subject_contact = 'aakrami@princeton.edu';
mainpth='/home/vizhe/Documents/Auditory_PWM_human/';
data_dir = 'home/vizhe/Documents/Auditory_PWM_human/'; %ALL EXPERIMENTS subject folders will be saved here!
addpath(mainpth)

Default_startup_pwm;

command_to_run = 'PWMHuman';
min_paid_session_length = 300;
show_fixation_in_prepause = 0;

modality = 'auditory';
fname_string = 'auditory_new_data';
a_volume = 0.4;
PWMHuman_gen_trial_mtx;
