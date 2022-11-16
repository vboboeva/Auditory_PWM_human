% Extended version of the auditory PWM Task for human subjects

% There's a number of settings parameters that should be in a separate
% file, tailored for each subject.  We expect the setup file for subject
% MyName to be MyName/MyName_startup.m
%
% some of these settings are optional; most are not
%
% check out Default_startup.m for a list of expected settings

%% parameters of task
% general parameters
% auditory stimulus presentation
a_bg_pre_post = [0 0];  % duration in sec of pre and post background white noise 
a_balance = 0;
a_srate = 40000;        % sample rate for sound production
a_fixation_box_size = [80 40];

v_Rad = [100 100];           % [x y] size in pixel of the stimulus shapes, which are filled polygons
v_LocationLeft = [-0.35, 0]; % position of left visual stimulus relative to center (0,0), in fractional screen space
v_LocationRight = [0.35, 0];
v_fixation_box_size = [80 40];
c_target_rad = [100 100];
c_Rad = [100 100];          % [x y] size in pixel of the stimulus shapes, which are filled polygons
c_Location = [0,  0];   % position of color stimulus relative to center (0,0), in fractional screen space
c_fixation_box_size = [20 20];
c_target_location = [0.25 0]; % position of the colored target shapes relative to center (0,0, in fractional screen space
c_transparency = 150;

% runtime parameters
antibias_beta = 3;      % antibias beta
antibias_tau = 20;
antibias_min_trials = 30;
perf_tau = 10;      % # of trials time constant for exponential window of recent performance

max_prepause = 1.5;
min_prepause = 0.75;
postpause = 0.1;

error_pause = 3;      % delay in sec before next trial after an error
hit_pause = 1;      % delay in sec before next trial after a hit
timeout_pause = 2;    % delay in sec before next trial after a timeout

max_rt = 2;             % reaction times longer than this count as timeout, invalid trial (in sec)

% task-feedback sounds
readytone = sin((0:1/a_srate:0.2)*a_srate/8)/10;  % the ready sound, a 200 msec pure tone
readytone = repmat(readytone, 2, 1)*a_volume;

[dingling, srate2] = audioread('DingLing.wav');    % the hit sound, dingling!
dingling = dingling(1:2:end); % compress sound by half
dingling = interp(dingling, ceil(a_srate/srate2))*0.5; 
dingling = repmat(0.3*dingling',2,1)*a_volume;  

[buzz, srate2] = audioread('Buzz02.wav');          % the miss sound, buzz!
buzz = interp(buzz(:,1), ceil(a_srate/srate2))*0.5;
buzz = repmat(0.5*buzz',2,1)*a_volume;                        

timeouttone = randn(size(readytone))*0.1*a_volume;

% setup preliminaries

% Store the original code in a string, so we can reproduce it fully.
fp = fopen([mfilename '.m'], 'r');
source_code.protocol = fscanf(fp, '%c');
fclose(fp);
try
    fp = fopen('make_pwm.m', 'r');
    source_code.make_pwm = fscanf(fp, '%c');
    fclose(fp);
    fp = fopen([data_dir name '_startup.m'], 'r');
    source_code.startup = fscanf(fp, '%c');
    fclose(fp);
    fp = fopen('Default_startup_pwm.m', 'r');
    if fp > 0
        source_code.defaults = fscanf(fp, '%c');
        fclose(fp);
    end
catch
    source_code.make_pwm = [];
    source_code.startup = [];
    fprintf(2, 'PWMHuman.m: could not save source code for make_pwm.m and the startup file\n\n');
end

% Make sure the directory exists
if ~exist([distr_type '/' name], 'dir')
   mkdir([distr_type '/' name]);
end

% --- find a session filename in the format name_fname_string_ddd where ddd is a session number
data_dir = ['/home/vizhe/Documents/Auditory_PWM_human/' distr_type '/' name '/'];

u = dir([data_dir name '_' fname_string '*.mat']);
session_number = numel(u)+1;

session_string = num2str(session_number);
session_string = ['0'*ones(1, 3-numel(session_string)) session_string];
session_fname  = [name '_' fname_string '_' session_string];


%% this needs to be resolved.
Screen('Preference', 'SkipSyncTests', 0);

% Set volume to a standrad level
% system('osascript -e "set Volume 4.2"');
system('amixer -c 0 set Master 4.2DB');

% Perform basic initialization of the sound driver:
InitializePsychSound;


% Open the default audio device [], with default mode [] (==Only playback),
% and a required latencyclass of zero 0 == no low-latency mode, as well as
% a frequency of freq and nrchannels sound channels.
% This returns a handle to the audio device:
try
    % Try with the 'freq'uency we wanted:
    pahandle = PsychPortAudio('Open', [], [], 0, a_srate, 2);
catch %#ok<*CTCH>
    % Failed. Retry with default frequency as suggested by device:
    fprintf('\nCould not open device at wanted playback frequency of %i Hz. Will retry with device default frequency.\n', a_srate);
    fprintf('Sound may sound a bit out of tune, ...\n\n');

    psychlasterror('reset');
    pahandle = PsychPortAudio('Open', [], [], 0, [], 2);
end

%  Data structure in which exptl data will be stored:
pwm_history = struct('modality', {}, 'a1_time', {}, 'a2_time', {}, 'fcut', {}, ...
   'filter_type', {}, ...
   'a1_sigma',{},'a2_sigma',{},'delaydur',{},...
   'prestim',{},'poststim',{},...
   'lfreq', {}, 'hfreq',{}, 'rt', {}, 'hit', {}, ...
   'wentR', {}, 'selside', {}, 'siderule', {}, 'trialnum', {}, ... % selside stands for selected size
   'valid', {}, 'T', {}, 'trialtime', {}, 'prepause', {});

prestims_real = [];
delaydurs_real = [];
poststims_real = [];
a1_times_real = [];
a2_times_real = [];

% fields = fieldnames(pwm_history)
  
% store date and time of this session
sessiondate = datetime("now");%
% , 31);

% these are internal temporary variables, clear them so they are not
% contaminated by any previous external values
wentR = NaN;
selside = 'N';
side  = NaN;

% these keys are exempt from being counted as a premature/violation key
% stroke
violation_exempt_keys = {'space', 'q'};


%% text instructions

insReady = 'Ready? Press space bar to go \n (Press q to quit)';
insPress = '... Press any key to continue ...';

insHit = 'Yes!';
insMiss = 'sorry... :(';
insTimeout = 'TimeOut'; 
insViolation = 'Premature responses are invalid';

insByeBye = 'All done! Press any key to quit.';

%% main task
ListenChar(2); % disable keyboard input to matlab; if there's an error and you don't have the keyboard, press Ctrl-c

try
    % Setting this preference to 1 suppresses the printout of warnings.
    oldEnableFlag = Screen('Preference', 'SuppressAllWarnings');


    % Get the screen up, collect all the info we need about it
    HideCursor;
    [w, wRect]=Screen('OpenWindow',0,0);
    Screen('TextSize', w, 25);
    wWidth = wRect(3); wHeight = wRect(4);
    ifi = Screen('GetFlipInterval', w);
    white = WhiteIndex(w);
    priorityLevel=MaxPriority(w);
    Priority(priorityLevel);
    
    % store specs of the screen refresh
    screen_specs.ifi = ifi;
    
    % prepare the vertices and colors of the shapes of the center, left and right target
    color0 = [0 1 0] * 255;
    color1 = [1 0 0] * 255;
    color2 = [0 0 1] * 255;
    target_CenterShape = makeShape('square', c_Rad, c_Location, wWidth, wHeight);
    target_LeftShape = makeShape('square', c_target_rad, -c_target_location, wWidth, wHeight);
    target_RightShape = makeShape('square', c_target_rad, c_target_location, wWidth, wHeight);

    
    DrawFormattedText(w, ['Please adjust your chair so that you are comfortable. \n ' insPress], ...
                     'center', 'center', white);
    Screen('Flip', w);
    KbWait();
    
    DrawFormattedText(w, ['"s" is the Left response, "l" is the Right response \n ' insPress], ...
                     'center', 'center', white);
    Screen('Flip', w);
    WaitSecs(0.2);
    KbWait();
    
    %% get the session parameters
    [a1_times,a2_times,a1_sigmas,a2_sigmas,siderules]=extract_from_struct(trl_mtx,'a1_time','a2_time','a1_sigma','a2_sigma','siderule');
    [prestims,poststims,delaydurs]=extract_from_struct(trl_mtx,'prestim','poststim','delaydur');

    trialnum = 1;
    n_done_trials = 0;
    avg_hitfrac = 0;
    weighted_hitfrac = 0;
    
    fixed_fixation_durs = [];
    myTs = [];

    while(trialnum < numel(prestims))
        % which modality are we in for this trial?
        switch modality
            case 'auditory'
                fixation_color    = white; % color of fixation cross
                fixation_box_size = a_fixation_box_size;
                fixation_shape    = 'box_cross';
                generate_sound    = 1; % whether the sound will actually be generated in call to make_pwm.m
                Screen('BlendFunction', w, GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA); % enable transparency
            otherwise
                fprintf(2, 'PWMHuman.m: does not recognize modality %s, quitting.\n', modality);
                break;
        end
        
        [hits, valids, selsides] = extract_from_struct(pwm_history, 'hit', 'valid', 'selside');
        hits = hits(valids==1); selsides = selsides(valids==1);
        
        
        if numel(hits) > antibias_min_trials
            kernel = (exp((-(numel(hits)-1):0)/antibias_tau))';
            lefts  = find(selsides=='l');
            rights = find(selsides=='r');
            hitl   = sum((hits(lefts) .*kernel(lefts))) /sum(kernel(lefts));
            hitr   = sum((hits(rights).*kernel(rights)))/sum(kernel(rights));
            %probs  = probabilistic_trial_selector([hitr, hitl], [0.5, 0.5], antibias_beta);
        end
        
        % Make noise sounds
        a1_time = a1_times(trialnum)/1000;
        a2_time = a2_times(trialnum)/1000;
        a1_sigma = a1_sigmas(trialnum);
        a2_sigma = a2_sigmas(trialnum);
        prestim = prestims(trialnum)/1000;
        poststim = poststims(trialnum)/1000;
        delaydur = delaydurs(trialnum)/1000;
        siderule = siderules(trialnum); %1 s2>s1 left - 0 s1>s2 right
        
        Fs=a_srate;
        T=max(a2_time,a1_time);
        [rawA1, rawA2, normA1, normA2]=noisestim(1,1,T,fcut,Fs,filter_type);
        modulator=singlenoise(1,T,[lfreq hfreq],Fs,'BUTTER');
        a1=normA1(1:a1_time*Fs).*modulator(1:a1_time*Fs).*a1_sigma*5;
        a2=normA2(1:a2_time*Fs).*modulator(1:a2_time*Fs).*a2_sigma*5;
        snd1 = [a1';a1'];
        snd2 = [a2';a2'];
        snd1 = balance_knob(a_balance, snd1);  % adjust for volume and balance
        snd2 = balance_knob(a_balance, snd2);  % adjust for volume and balance
        snd1(snd1(:)> 0.995) = 0.995;
        snd1(snd1(:)<-0.995) = -0.995;
        snd2(snd2(:)> 0.995) = 0.995;
        snd2(snd2(:)<-0.995) = -0.995;  
        

        % PROBLEM
        % Once myT has reached the maximum possible value that it can take,
        % fixed_fixation_dur is set to that value + 0.5, and it will never
        % decrease below that value (because we will always have 
        % fixed_fixation_dur > myT)


        myT = a1_time + a2_time + delaydur + prestim + poststim;

        % if fixed_fixation_dur is shorter than longer stimulus, then
        % increase fixed_fixation_dur
        if fixed_fixation_dur>0 && fixed_fixation_dur <= myT
            fixed_fixation_dur = myT + 0.5;
        end
        
        if ~fixed_fixation_dur
            % Choose a prepause duration:
            prepause = rand(1)*(max_prepause-min_prepause)+min_prepause;
        else
            prepause = max(fixed_fixation_dur - myT, 0);
        end

        % save all stimulus info for this trial
        pwm_history(trialnum).fcut   = fcut;
        pwm_history(trialnum).lfreq   = lfreq;
        pwm_history(trialnum).hfreq   = hfreq;
        pwm_history(trialnum).filter_type   = filter_type;
        pwm_history(trialnum).modality   = modality;
        pwm_history(trialnum).T          = myT;
        pwm_history(trialnum).trialnum   = trialnum;
        pwm_history(trialnum).prepause   = prepause;
        pwm_history(trialnum).prestim   = prestim;
        pwm_history(trialnum).poststim   = poststim;
        pwm_history(trialnum).a1_time   = a1_time;
        pwm_history(trialnum).a2_time   = a2_time;
        pwm_history(trialnum).a1_sigma   = a1_sigma;
        pwm_history(trialnum).a2_sigma   = a2_sigma;
        pwm_history(trialnum).delaydur   = delaydur;
        
        
        % READY?
        %PsychPortAudio('FillBuffer', pahandle, readytone);
        DrawFormattedText(w, sprintf('Trial %i', trialnum), 'center', round(wHeight*0.1), white);
        if feedback_interval == 1
            DrawFormattedText(w, sprintf('%.3f', avg_hitfrac), round(wWidth*0.8), round(wHeight*0.1), white);
            DrawFormattedText(w, sprintf('%.3f', weighted_hitfrac), round(wWidth*0.85), round(wHeight*0.1), white);
        elseif mod(trialnum-1, feedback_interval) == 0
            DrawFormattedText(w, sprintf('Mean accuracy this session: %.3f', avg_hitfrac), round(wWidth*0.7), round(wHeight*0.1), white);
        end
            
        
        DrawFormattedText(w, insReady, 'center', 'center', white);
        Screen('Flip', w);
        
        %PsychPortAudio('Start', pahandle, 1, 0, 1); 
        WaitSecs(0.02);
        
        
        responded = 0;
        while ~responded
            [keyIsDown, t_start, keyCode] = KbCheck();
            keyname = KbName(keyCode);
            if keyIsDown && any(ismember(keyname, {'space', 'q'}))
                responded = 1;
            end
            WaitSecs(0.005);
        end        
        
            
        % timestamp of when this trial started
        pwm_history(trialnum).timestamp  = rem(now,1); % current time

        
        % if subject is done and tries to quit
        if any(strcmp(keyname, 'q'))
            break;
        end
        
        % put up fixation cross at center of screen, display trialnum
        DrawFormattedText(w, sprintf('Trial %i', trialnum), 'center', round(wHeight*0.1), white);
        if show_fixation_in_prepause
            DrawFixationCross(w, fixation_box_size, wWidth, wHeight, fixation_color, fixation_shape);
        end
        Screen('Flip', w);
        WaitSecs(0.1);
        
        premature_response = 0;
        t_start_fixation = GetSecs;
        switch modality
            case 'auditory'

                start = GetSecs;
                %% center target on for prestim dur
                DrawFormattedText(w, sprintf('Trial %i', trialnum), 'center', round(wHeight*0.1), white);
                Screen('FillPoly', w, color1, target_LeftShape);            
                Screen('Flip', w); 
                WaitSecs(prestim);
                prestims_real(end+1) = GetSecs - start;
                
                %% play the first sound
                start = GetSecs;
                DrawFormattedText(w, sprintf('Trial %i', trialnum), 'center', round(wHeight*0.1), white);
                Screen('FillPoly', w, color1, target_LeftShape);            
                Screen('Flip', w); 
                PsychPortAudio('FillBuffer', pahandle, snd1);
                a1LengthFrames = round(a1_time/ ifi);
                PsychPortAudio('Start', pahandle, 1, 0, 1);             
                for i = 1:a1LengthFrames
                    DrawFormattedText(w, sprintf('Trial %i', trialnum), 'center', round(wHeight*0.1), white);
                    Screen('FillPoly', w, color1, target_LeftShape);            
                    Screen('Flip', w);
                    if punish_early_responses
                        responded = 0;
                        % ketIsDown -> whether any key has been pressed
                        % t_premature_response -> time when key was pressed
                        % keyCode -> which key was it
                        [keyIsDown, t_premature_response, keyCode] = KbCheck();
                        keyname = KbName(keyCode);
                        if keyIsDown && ~any(ismember(keyname, violation_exempt_keys))
                            responded = 1;
                            premature_response = 1;
                            premature_keyname = KbName(keyCode);
                            PsychPortAudio('Stop', pahandle);

                            DrawFormattedText(w, sprintf('Trial %i', trialnum), 'center', round(wHeight*0.1), white);
                            DrawFormattedText(w, sprintf('You responded before the end of the stimulus!'), 'center', 'center', white);
                            Screen('Flip', w);
                        end
                    end
                end
                PsychPortAudio('Stop', pahandle);
                a1_times_real(end+1) = GetSecs - start;
               
                %% delay interval
                start = GetSecs;
                Screen('Flip', w);
                frame_durations = [];
                first_part_durations = [];
                second_part_durations = [];
                lhs = [];
                rhs = [];
                delayLengthFrames = round(delaydur/ ifi);
                for i = 1:delayLengthFrames
                    DrawFormattedText(w, sprintf('Trial %i', trialnum), 'center', round(wHeight*0.1), white);
                    DrawFormattedText(w, sprintf('WAIT!'), 'center', 'center', white);
                    Screen('Flip', w);
                    if punish_early_responses
                        responded = 0;
                        % ketIsDown -> whether any key has been pressed
                        % t_premature_response -> time when key was pressed
                        % keyCode -> which key was it
                        [keyIsDown, t_premature_response, keyCode] = KbCheck();
                        keyname = KbName(keyCode);
                        if keyIsDown && ~any(ismember(keyname, violation_exempt_keys))
                            responded = 1;
                            premature_response = 1;
                            premature_keyname = KbName(keyCode);
                            PsychPortAudio('Stop', pahandle);

                            DrawFormattedText(w, sprintf('Trial %i', trialnum), 'center', round(wHeight*0.1), white);
                            DrawFormattedText(w, sprintf('You responded before the end of the stimulus!'), 'center', 'center', white);
                            Screen('Flip', w);
                        end
                    end
                end
                delaydurs_real(end+1) = GetSecs - start;

                %% play the second sound
                start = GetSecs;
                DrawFormattedText(w, sprintf('Trial %i', trialnum), 'center', round(wHeight*0.1), white);
                Screen('FillPoly', w, color2, target_RightShape);            
                Screen('Flip', w); 
                PsychPortAudio('FillBuffer', pahandle, snd2);
                PsychPortAudio('Start', pahandle, 1, 0, 1);             

                a2LengthFrames = round(a2_time/ ifi);
                for i = 1:a2LengthFrames
                    DrawFormattedText(w, sprintf('Trial %i', trialnum), 'center', round(wHeight*0.1), white);
                    Screen('FillPoly', w, color2, target_RightShape);            
                    Screen('Flip', w);                 
                    if punish_early_responses
                        responded = 0;
                        % ketIsDown -> whether any key has been pressed
                        % t_premature_response -> time when key was pressed
                        % keyCode -> which key was it
                        [keyIsDown, t_premature_response, keyCode] = KbCheck();
                        keyname = KbName(keyCode);
                        if keyIsDown && ~any(ismember(keyname, violation_exempt_keys))
                            responded = 1;
                            premature_response = 1;
                            premature_keyname = KbName(keyCode);
                            PsychPortAudio('Stop', pahandle);

                            DrawFormattedText(w, sprintf('Trial %i', trialnum), 'center', round(wHeight*0.1), white);
                            DrawFormattedText(w, sprintf('You responded before the end of the stimulus!'), 'center', 'center', white);
                            Screen('Flip', w);
                        end
                    end               
                end
                PsychPortAudio('Stop', pahandle);
                a2_times_real(end+1) = GetSecs - start;

           
%                 %% play the second sound (ORIGINAL CODE)
%                 start1 = GetSecs;
%                 DrawFormattedText(w, sprintf('Trial %i', trialnum), 'center', round(wHeight*0.1), white);
%                 Screen('FillPoly', w, color2, target_RightShape);            
%                 Screen('Flip', w);
%                 PsychPortAudio('FillBuffer', pahandle, snd2);
%                 PsychPortAudio('Start', pahandle, 1, 0, 1);
%                 if punish_early_responses
%                     responded = 0;
%                     [keyIsDown, t_premature_response, keyCode] = KbCheck();
%                     keyname = KbName(keyCode);
%                     if keyIsDown && ~any(ismember(keyname, violation_exempt_keys))
%                         responded = 1;
%                         premature_response = 1;
%                         premature_keyname = KbName(keyCode);
%                         PsychPortAudio('Stop', pahandle);
%                         
%                         DrawFormattedText(w, sprintf('Trial %i', trialnum), 'center', round(wHeight*0.1), white);
%                         DrawFormattedText(w, sprintf('You responded before the end of the stimulus!'), 'center', 'center', white);
%                         Screen('Flip', w);
%                     end
%                 end
%                 a2_times_real(end+1) = GetSecs - start1;

                %% center target on for poststim dur
                start = GetSecs;
                DrawFormattedText(w, sprintf('Trial %i', trialnum), 'center', round(wHeight*0.1), white);
                Screen('FillPoly', w, color2, target_RightShape);            
                Screen('Flip', w); 
                WaitSecs(poststim);
                poststims_real(end+1) = GetSecs - start;
                 
            otherwise
                
                fprintf(2, 'PWMHuman.m: dont know about modality %s, quitting\n', modality);
                break;                
        end
        
        % stimulus presentation is done, take the fixation cross away
        DrawFormattedText(w, sprintf('Trial %i', trialnum), 'center', round(wHeight*0.1), white);
        Screen('Flip', w);
        
        % wait for response
        t_wait_for_response = GetSecs;
        if premature_response
            t_response = t_premature_response;
            response = 'INVALID';
        else
            responded = 0;
            while ~responded
                if GetSecs - t_wait_for_response > max_rt
                    % if we've passed the deadline, count the response as a
                    % fail
                    responded = 1;
                    response = 'INVALID';
                else
                    [keyIsDown, t_response, keyCode] = KbCheck();
                    keyname = KbName(keyCode);
                    if keyIsDown && any(ismember(keyname, {'r', 'a', 's', 'l', ';', 'd', 'i'}))
                        responded = 1;
                        if numel(keyname) > 1
                            % if more than one key was pressed, count the first one
                            keyname = keyname{1};
                        end

                        if ismember(keyname, {'a', 's', 'd'})
                            response = 'LEFT';
                            wentR = 0;
                            selside  = 'L';
                        elseif ismember(keyname, {'l', ';', 'i'})
                            response = 'RIGHT';
                            wentR = 1;
                            selside  = 'R';
                        elseif ismember(keyname, 'r')
                            response = 'INVALID';
                            wentR = NaN;
                            selside  = 'r';
                        end
                    end
                end
                WaitSecs(0.005); % we don't have to go as fast as possible
            end
        end
        
        % Store outcome for this trial
        pwm_history(trialnum).trialtime = GetSecs - t_start;
        pwm_history(trialnum).wentR     = wentR;
        pwm_history(trialnum).selside      = selside;
        pwm_history(trialnum).siderule      = siderule;
        pwm_history(trialnum).rt        = t_response - t_wait_for_response; % premature violation trials have rt<0
        pwm_history(trialnum).valid     = ~strcmp(response, 'INVALID'); % early and late response trials are invalid
        if pwm_history(trialnum).valid == 0
            if trialnum ==1
            pwm_history(trialnum).hit = 1;   
            else
            pwm_history(trialnum).hit = NaN;
            end
        else
            pwm_history(trialnum).hit  = double(xor(siderule, wentR));
        end
        
        Screen('Close');

        %sca;
        
        % save data after every trial
        n_done_trials = n_done_trials + 1;
        save([data_dir session_fname], 'pwm_history', 'sessiondate', 'source_code');
        

        % feedback to the subject
        if pwm_history(trialnum).hit == 1
            report = 'correct';
        elseif pwm_history(trialnum).hit == 0
            report = 'miss';
        elseif pwm_history(trialnum).hit == 0.5
            % if there's no right answer (equal events on left and right,
            % report randomly
            if rand(1) > 0.5 
                report = 'correct';
            else              
                report = 'miss'; end
        else
            % an invalid trial
            report = 'timeout';
        end
        
        % report mean hit fraction and moving average hit fraction for
        % valid trials
        x = cell(numel(pwm_history),1);
        [x{:}] = deal(pwm_history.hit); x = cell2mat(x);
        if sum(~isnan(x))>0
            avg_hitfrac = mean(x);
            eweights = exp(-(trialnum - (1:trialnum)')/perf_tau);
            eweights(isnan(x)) = 0;
            eweights = eweights/sum(eweights);
            weighted_hitfrac = sum(eweights.*x);
        else
            avg_hitfrac = 0;
            weighted_hitfrac = 0;
        end
        
        % report to the screen and provide auditory feedback
        DrawFormattedText(w, sprintf('Trial %i', trialnum), 'center', round(wHeight*0.1), white);
        if feedback_interval == 1
            DrawFormattedText(w, sprintf('%.3f', avg_hitfrac), round(wWidth*0.8), round(wHeight*0.1), white);
            DrawFormattedText(w, sprintf('%.3f', weighted_hitfrac), round(wWidth*0.85), round(wHeight*0.1), white);
        elseif mod(trialnum, feedback_interval) == 0
            DrawFormattedText(w, sprintf('Mean accuracy this session: %.3f', avg_hitfrac), round(wWidth*0.7), round(wHeight*0.1), white);
        end
        if pwm_history(trialnum).rt <= 0 % premature response
            DrawFormattedText(w, insViolation, 'center', 'center', white);
            Screen('Flip', w);
            PsychPortAudio('FillBuffer', pahandle, timeouttone);
            PsychPortAudio('Start', pahandle, 1, 0, 1);
            WaitSecs(timeout_pause);
        elseif pwm_history(trialnum).valid == 0 % timed out on response
            DrawFormattedText(w, insTimeout, 'center', 'center', white);
            Screen('Flip', w);
            PsychPortAudio('FillBuffer', pahandle, timeouttone);
            PsychPortAudio('Start', pahandle, 1, 0, 1);
            WaitSecs(timeout_pause);
        elseif strcmp(report, 'correct') % hit
            DrawFormattedText(w, insHit, 'center', 'center', white);
            Screen('Flip', w);
            %PsychPortAudio('FillBuffer', pahandle, dingling);
            %PsychPortAudio('Start', pahandle, 1, 0, 1);
            WaitSecs(max([hit_pause size(dingling,2)/a_srate]));            
        else % miss
            DrawFormattedText(w, insMiss, 'center', 'center', white);
            Screen('Flip', w);
            %PsychPortAudio('FillBuffer', pahandle, buzz);
            %PsychPortAudio('Start', pahandle, 1, 0, 1);
            WaitSecs(max([error_pause size(buzz,2)/a_srate]));
        end
            
        trialnum = trialnum + 1;

    prestims_real
    delaydurs_real
    poststims_real
    a1_times_real
    a2_times_real    

%         WaitSecs(2);
%         sca;

    end % end of trials loop

    % add ending time to session
    sessiondate(2,:) = datetime("now");%, 31);
    
    
    % end of session accounting
    [hits, timestamps] = extract_from_struct(pwm_history, 'hit', 'timestamp');
    timestamps(end+1) = rem(now,1);
    session_duration = (datetime(sessiondate(2,:)) - datetime(sessiondate(1,:)))*24; % in hours
    
    % how much time was spent working?
    t_trials = diff(timestamps)*24; % in hours
    working_duration = sum(t_trials(t_trials < 1/60)); % all time spent working, minus breaks of longer than 1 min
    
    
    % figure out the payscale for this session based on the mean hitfrac
    mypayscale = 0.0231;
    if strcmp(modality, 'auditory')
        payscale = a_payscale;
    elseif strcmp(modality, 'orientation')
        payscale = o_payscale;
    else
        payscale = [];
    end
    if ~isempty(payscale)
        th = find(payscale(1,:) < avg_hitfrac, 1, 'last');
        mypayscale = payscale(2, th);
    end
    
    % if we're displaying the session summary:
    net_pay = 0;
    if display_session_summary        
        insByeBye = 'Thank you for completing this session!\n'; 
        insByeBye = [insByeBye sprintf('\n\nThe session was %s min (%s working min) in duration and had a mean accuracy of %5.2f.\n', round(session_duration*60), round(working_duration*60), avg_hitfrac)];
        
        if trialnum < min_paid_session_length
            net_pay = 12*working_duration;
            insByeBye = [insByeBye sprintf('You completed less than %i trials, so will be paid $%5.2f for this session.', min_paid_session_length, net_pay)];
        else
            net_pay = sum(hits==1)*mypayscale;
            insByeBye = [insByeBye sprintf('You had %i correct answers, and at $%5.4f per answer, you will be paid $%5.2f for this session.\n', sum(hits==1), mypayscale, net_pay)];
            insByeBye = [insByeBye sprintf('You made approximately $%5.2f per hour.\n', net_pay/working_duration)];
        end
        insByeBye = [insByeBye '\n\n... Press any key to end this session....\n'];
    end
    
    DrawFormattedText(w, insByeBye, 'center', 'center', white);
    Screen('Flip', w);
    WaitSecs(0.2);
    KbWait();

    
    
    % ---- send automated email with summary of this session ---
    % assemblem the body of the email message
    msg = cell(0);
    msg{end+1} = sprintf('%s just completed a session of %s, in version %s,\n', name, command_to_run, modality);
    msg{end+1} = '';
    msg{end+1} = sprintf('Session started at %s and ended at %s.\n', sessiondate(1,:), sessiondate(2,:));
    msg{end+1} = '';
    msg{end+1} = sprintf('The total session duration was %s minutes, or %s minutes excluding breaks\n', round(session_duration*60), round(working_duration*60));
    msg{end+1} = sprintf('%i trials were done at an average hitfrac of %5.2f.\n', trialnum, avg_hitfrac);
    msg{end+1} = '';
    if trialnum < min_paid_session_length
        msg{end+1} = sprintf('%s completed less than %i trials, so will be paid $%5.2f for this session', name, min_paid_session_length, net_pay);
    else
        msg{end+1} = sprintf('%s had %i correct answers, and at $%5.4f per answer, the net pay for this session will be $%5.2f,', name, sum(hits==1), mypayscale, net_pay);
        msg{end+1} = sprintf('or approximately $%5.2f per hour.\n', net_pay/working_duration);
    end
    msg{end+1} = '';
    msg{end+1} = '';
    msg{end+1} = sprintf('The data file should be on svn: %s.mat\n', session_fname);
    if send_automated_email
        try            
            % send the email message
            setpref('Internet','SMTP_Server','sonnabend.princeton.edu');
            setpref('Internet','E_mail','PsychophysicsReport@Princeton.EDU');
            email_recipients = experimenter_contacts;
            if ~isempty(subject_contact) % add subject to email list
                email_recipients{end+1} = subject_contact;
            end
            sendmail(email_recipients, sprintf('Human Psychophysics Report: %s', name), msg);
            
            fprintf(1, '\n\n The automated email was sent successfully.\n');
            email_status = 'email sent';
        catch
            fprintf(1, '\n\n There was an error in sending the automated email.\n');
            email_status = 'error in sending email';
        end
    else
        fprintf(1, '\n\n Automated email was not sent.\n'); 
        email_status = 'email not sent';
    end
    % ----------
    
    
    % save data
    pwm_history = pwm_history(1:n_done_trials);
    save([data_dir session_fname], 'pwm_history', 'sessiondate', 'source_code', 'screen_specs', 'email_status', 'msg', 'net_pay');

    % Tidying up at the end
    PsychPortAudio('Close', pahandle); % Close the audio device:
    Screen('close');
    Screen('CloseAll');
    Priority(0);
    ShowCursor;
    
    % ---- send data file to svn ---
%     if commit_to_svn
%         [status, res] = system(sprintf('svn add %s', [data_dir session_fname '.mat'])); 
%         if status
%             fprintf(2, 'Error %s while adding data file to svn!  Please let Bing know.\n', res);
%         end
%         
%         [status, res] = system(sprintf('svn commit -m "" %s', [data_dir session_fname '.mat']));
%         if status
%             fprintf(2, 'Error %s while commiting data file to svn!  Please let Bing know.\n', res);
%         end
%         
%         [status, res] = system(sprintf('svn commit -m "" %s', [data_dir name '_startup.m']));
%         if status
%             fprintf(2, 'Error %s while commiting startup file to svn!  Please let Bing know.\n', res);
%         end
%     else
%         fprintf(1, '\n\n Data file has not been committed to svn\n\n\n');
%     end  
    % ----------
    
    % At the end of your code, it is a good idea to restore the old level.
    Screen('Preference','SuppressAllWarnings',oldEnableFlag);

catch   % If there's an error, save the data, quit and tell us about it
    ListenChar(0); % re-enable keyboard input to matlab    
    save([data_dir session_fname], 'pwm_history', 'sessiondate', 'source_code');
    PsychPortAudio('Close', pahandle); % Close the  audio device:
    Screen('CloseAll');
    Priority(0);
    psychrethrow(psychlasterror);

end

ListenChar(0); % re-enable keyboard input to matlab
