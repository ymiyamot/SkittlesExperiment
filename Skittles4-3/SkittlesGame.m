function SkittlesGame(subject_name, dayNum)

    %%% On Day 2, no baseline blocks %%%
    if dayNum == 1,
        dirName = subject_name;
        baseblocks = 2;
    else
        dirName = [subject_name, 'Day2']; 
        baseblocks = 0;
    end
    mkdir(dirName);
    totalblocks = 30;
    blockLen = 50;
    
    %%% Baseline blocks with no feedback %%%
    for bl = 1:baseblocks,
        fprintf('Block %d: baseline', bl);
        keyboard;
        Tester(subject_name, bl, blockLen, 1, dayNum, dirName);
    end
    
    %%% Training blocks with feedback %%%
    for bl=baseblocks + 1:totalblocks,
        fprintf('Block %d: training', bl);
        keyboard;
        Tester(subject_name, bl, blockLen, 0, dayNum, dirName);
    end
end

function Tester(subject_name, blocknum, total_trials, baseMode, dayNum, dirName)
try

    
    % Variable for recording data
    tabletRec = NaN(1000000, 4);
    tabletCt = 0;

    % Simulation parameters
    feedbackMode = 4;

    if baseMode == 1,
        targloc = [NaN NaN];
        feedbackOn = 0;
    else
        targloc = [1100 + 75 * sqrt(3), 270];
        feedbackOn = 1;
    end
    
    % Initialization
    init_t = getSecs(); % Session start time 
    LL = FreeLaunch(960 + 150, 540 + 250, 960 + 50 - 185.71, 540 + 250, targloc);
    S = Sample(LL.w(1), LL.resx, LL.resy, 1, init_t);
    trial = 1;
    simTrigger = false;
    simRun = false;
    waitStart = NaN;
    ITIStart = NaN;
    ITI = NaN;
    launchTime = NaN;
    delayFeedbackStart = NaN;
    numFails = 0;
    releaseTime = 0;
    feedbackGiven = 0;

    drawscore = 0;
    L = 12000; % pole height
    explodeTime = 2; % explosion time
    dt = 0.01; % time step size
    total_time = 4; % total time to calculate trajectory
    t = 0:dt:total_time;
    waitTime = 0; % This changes according to type of feedback
    
    % Load audio
    ting = load('audiofiles/ting');
    ting.ting = ting.ting(401:end, :)';
    
    % retrieve history for history-dependent reward
    histwin = 25;
    % No history for baseline blocks and initial block.
    if (dayNum == 1 && blocknum == 3) ...
            || baseMode == 1 ...
            || (dayNum == 2 && blocknum == 1),
        threshbuffer = NaN(50, 1);
        threshbuffer(50) = 200;
    else
        rec = csvimport([dirName, '/', subject_name, '_trialrecord_', ...
            char('a' + blocknum - 2), '.csv']);
        distHist = rec(2:end, 10);
        recVec = NaN(size(distHist));
        
        % Parsing for NaNs
        for dist_i = 1:length(distHist),
            if ~isnan(distHist{dist_i}) && all('string' == class(distHist{dist_i}))
                recVec(dist_i) = str2num(distHist{dist_i});
            elseif ~isnan(distHist{dist_i}) && all('double' == class(distHist{dist_i}))
                recVec(dist_i) = distHist{dist_i};
            end
        end
        threshbuffer = recVec;
    end

    % Storage for trial performance information
    header = {'release_x', 'release_y', 'release_vx', 'release_vy', ...
        'targ_x', 'targ_y', 'feedbackMode', 'explodeTime', ...
        'pendulumLength', 'distFromTarg', 'releaseTime', ...
        'releaseInd', 'feedbackStartTime', 'feedbackStartInd', ...
        'feedbackEndTime', 'feedbackEndInd', 'waitTime', 'rwd', ...
        'rwdThresh', 'historyWindowLength', 'numFails', ...
        'startTime', 'startInd'};
    record = NaN(total_trials, length(header));
    
    [~,~, mouse] = GetMouse; % Detect if mouse was pressed for aborting experiment
    while(~mouse(3) && trial <= total_trials)
        curr_t = getSecs() - init_t;
        
        S.Update;
        if(S.samp)
            for samp_i = 1:length(S.x)
                
                % Record measurements
                tabletCt = tabletCt + 1;
                tabletRec(tabletCt, :) = [S.x(samp_i), S.y(samp_i), S.t(samp_i), curr_t];
                
            end
            
            % Update the laucher's velocity samples
            LL.UpdateBuffer(S);
            
            % Update the state of the launcher
            if isnan(ITIStart) && isnan(waitStart) && LL.state == 0 && abs(S.x(samp_i) - (LL.startBarX)) <= LL.barwidStart / 2 && abs(S.y(samp_i) - LL.startBarY) <= LL.barheiStart / 2 && LL.v < LL.threshvelStart,
                
                % Turned cursor green and is ready for launch
                disp('Ready to launch. state 1')
                LL.state = 1;
                startTime = curr_t;
                startInd = tabletCt;
            elseif LL.state == 1 && (S.x(samp_i) < LL.midBarX && (abs(S.y(samp_i) - LL.startBarY) <= LL.barheiLaunch / 2))
                disp('Launched ball. state 2')
                
                LL.state = 2;
                
                % Time of release
                launchTime = getSecs;
                releaseTime = launchTime - init_t;
                releaseInd = tabletCt;

                simRun = true;
                simTrigger = true;

            elseif LL.state == 2 ...
                    && ((~(abs(S.x(samp_i) - LL.launchBarX) <= LL.barwidLaunch / 2 ...
                    && abs(S.y(samp_i) - LL.launchBarY) <= LL.barheiLaunch / 2) ...
                    && LL.v < LL.threshvelStop) ...
                    || (S.x(samp_i) < LL.launchBarX - LL.barwidLaunch / 2) ...
                    || abs(S.y(samp_i) - LL.launchBarY) > LL.barheiLaunch),
                
                Snd('Quiet');
                % If the hand doesn't end up at the launch bar, then abort trial.
                disp('Did not reach launch bar. state 0')
                disp((~(abs(S.x(samp_i) - LL.launchBarX) <= LL.barwidLaunch / 2 && abs(S.y(samp_i) - LL.launchBarY) <= LL.barheiLaunch / 2) && LL.v < LL.threshvelStop))
                disp(abs(S.x(samp_i) - LL.launchBarX) <= LL.barwidLaunch / 2)
                disp(abs(S.y(samp_i) - LL.launchBarY) <= LL.barheiLaunch / 2)
                disp(LL.v < LL.threshvelStop)
%                 disp((S.x(samp_i) < LL.launchBarX - LL.barwid * 4))
%                 disp(abs(S.y(samp_i) - LL.launchBarY) > LL.barhei)

                LL.state = 0;
                simTrigger = false;
                simRun = false;
                launchTime = NaN;
                delayFeedbackStart = NaN;
                
                % Clear trajectory off screen
                LL.Recopy(3, 2);
                LL.Recopy(2, 1);
                sound(0.1*sin(1:500));
                numFails = numFails + 1;
                
                feedbackGiven = 0; %This is for feedbackMode = 1
                waitStart = NaN;
                ITIStart = NaN;
            elseif LL.state == 1 && ~(abs(LL.buffer(end, 1) - LL.startBarX) <= LL.barwidStart / 2 && abs(LL.buffer(end, 2) - LL.startBarY) <= LL.barheiStart / 2) && LL.v < LL.threshvelStop,
                
                % Slowed down before launch so regress to pre-launch state
                disp('Go back to preparation. state 0')
                LL.state = 0;
                
            end
        end
        
        % Calculate outcome of simulation
        if simTrigger,
            sound(0.1*sin(1.5 * (1:500)));
            
            simTrigger = false; % only run this calculation once per launch
            simCt = 1;
            
            % variables that affect outcome
            vx = diff(LL.buffer(end - 1:end, 1)) / diff(LL.buffer(end - 1:end, 3)); % pixel units, maybe smooth this more.
            vy = -diff(LL.buffer(end - 1:end, 2)) / diff(LL.buffer(end - 1:end, 3)); % pixel units
            px = S.x(samp_i);
            py = S.y(samp_i);

            rwdthresh = nanmedian(threshbuffer(max(1, length(threshbuffer) - histwin + 1):end));
%             [x, y, d, score] = FB.SkittlesSimulation(xr, yr, vx, vy, speed, L, explode_time, wait_time, drawscore, rwdthresh, feedback_on);
%             %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

          % Rotate the coordinates to adjust for screen coordinates
            xr_cent = (px - LL.cx); % pixel units
            yr_cent = -(py - LL.cy); % pixel units
            z = -sqrt((L^2 - xr_cent^2 - yr_cent^2));
            
            disp([xr_cent, yr_cent, z]);
            disp([vx, vy, 0])
            
            % Calculate resulting trajectory of ball
            pos = PenDynamics([xr_cent, yr_cent, z], [vx, vy, 0], dt, total_time, L);

            % Convert coordinates
            x = pos(:, 1) + LL.cx; % pixel units (uncentered)
            y = -pos(:, 2) + LL.cy; % pixel units (uncentered)
            % Create another set of coordinates for detecting off-screen
            % trajectories
            xOff = pos(:, 1) * LL.scaleFac + LL.cx; % pixel units (uncentered)
            yOff = -pos(:, 2) * LL.scaleFac + LL.cy; % pixel units (uncentered)
            
            % Buzz sound if ball will fly off screen
            if any(xOff < 0) || any(xOff > 1920) || any(yOff < 0) || any(yOff > 1080)
               sound(0.75 * sin((1:3000)/12))
            end

            % Calculate distance of each point from target
            dists = sqrt(bsxfun(@minus, x, targloc(1)).^2 + bsxfun(@minus, y, targloc(2)).^2);
            
            % Distance of ball at explosion time
            explodeInd = explodeTime / dt + 1;
            d = dists(explodeInd);
            
            %%%%%%%%%%%%%%%%%%%%%%%%

            % Update history
            threshbuffer = [threshbuffer(2:end); d];
        end
        
        % Animate simulation
        if simRun,
            score = d;
            
            if feedbackMode == 1,
                %%%%%% Instant path feedback simulation %%%%%%
                
                if feedbackGiven == 0,
                    if feedbackOn,
                        for traj_i = 2:explodeInd,
                            Screen('DrawLine',  LL.w(2), LL.green, ...
                                LL.scalePresentX(x(traj_i - 1)), ...
                                LL.scalePresentY(y(traj_i - 1)), ...
                                LL.scalePresentX(x(traj_i)), ...
                                LL.scalePresentY(y(traj_i)));
                            Screen('DrawLine',  LL.w(2), LL.green, ...
                                LL.resx - LL.scalePresentX(x(traj_i - 1)), ...
                                LL.resy - LL.scalePresentY(y(traj_i - 1)), ...
                                LL.resx - LL.scalePresentX(x(traj_i)), ...
                                LL.resy - LL.scalePresentY(y(traj_i)));
                        end
                    else
                        for traj_i = 2:explodeInd,
                            Screen('DrawLine',  LL.w(2), LL.green, ...
                                LL.scalePresentX(x(traj_i - 1)), ...
                                LL.scalePresentY(y(traj_i - 1)), ...
                                LL.scalePresentX(x(traj_i)), ...
                                LL.scalePresentY(y(traj_i)));
                        end
                    end
                    Screen('FillOval', LL.w(2), LL.green, ...
                        CenterRectOnPoint([0, 0, 15, 15], ...
                        LL.scalePresentX(x(explodeInd)), ...
                        LL.scalePresentY(y(explodeInd))));
                    if feedbackOn,
                        Screen('FillOval', LL.w(2), LL.green, ...
                            CenterRectOnPoint([0, 0, 15, 15], ...
                            LL.resx - LL.scalePresentX(x(explodeInd)), ...
                            LL.resy - LL.scalePresentY(y(explodeInd))));
                    end
                    
                    disp(sprintf('time since launch: %f', getSecs - launchTime))
                    disp('Finished simulation. state 0')
                    
                    sound(0.1*sin(2 * (1:500)));
                    
                    waitStart = getSecs;
                    waitTime = 1;
                    ITIStart = getSecs;
                    ITI = 3;
                    feedbackInd = tabletCt;
%                     LL.state = 0;
                    launchTime = NaN;
                    feedbackGiven = 1;
                else % if feedback was given already, then detect bad movement endings and beep and indicate
                    if numFails ~= -1 && ((~(abs(S.x(samp_i) - LL.launchBarX) <= LL.barwidLaunch / 2 && abs(S.y(samp_i) - LL.launchBarY) <= LL.barheiLaunch / 2) && LL.v < LL.threshvelStop) ...
                    || (S.x(samp_i) < LL.launchBarX - LL.barwidLaunch / 2) ...
                    || abs(S.y(samp_i) - LL.launchBarY) > LL.barheiLaunch),
%                         sound(0.1*sin(1:500)); % Bad beep
%                         numFails = -1; % Failed trial
%                         
%                         % Clear trajectory off screen
%                         LL.Recopy(3, 2);
%                         LL.Recopy(2, 1);
                    elseif LL.v < LL.threshvelStop || getSecs - waitStart > waitTime, % Don't go into simRun if wait time is over.
                        simRun = false;
                        feedbackGiven = 0;
                    end
                end 
            elseif feedbackMode == 1.5 && (abs(S.x(samp_i) - LL.launchBarX) <= LL.barwidLaunch / 2 && abs(S.y(samp_i) - LL.launchBarY) <= LL.barheiLaunch / 2) && LL.v < LL.threshvelStop
                %%%%%% Instant path feedback simulation that waits until movement is stopped %%%%%%
                if feedbackOn,
                    for traj_i = 2:explodeInd,
                        Screen('DrawLine',  LL.w(2), LL.green, LL.scalePresentX(x(traj_i - 1)), LL.scalePresentY(y(traj_i - 1)), LL.scalePresentX(x(traj_i)), LL.scalePresentY(y(traj_i)));
                        Screen('DrawLine',  LL.w(2), LL.green, LL.resx - LL.scalePresentX(x(traj_i - 1)), LL.resy - LL.scalePresentY(y(traj_i - 1)), LL.resx - LL.scalePresentX(x(traj_i)), LL.resy - LL.scalePresentY(y(traj_i)));
                    end
                else
                    for traj_i = 2:explodeInd,
                        Screen('DrawLine',  LL.w(2), LL.green, LL.scalePresentX(x(traj_i - 1)), LL.scalePresentY(y(traj_i - 1)), LL.scalePresentX(x(traj_i)), LL.scalePresentY(y(traj_i)));
                    end
                end
                Screen('FillOval', LL.w(2), LL.green, CenterRectOnPoint([0, 0, 15, 15], LL.scalePresentX(x(explodeInd)), LL.scalePresentY(y(explodeInd))));
                if feedbackOn,
                    Screen('FillOval', LL.w(2), LL.green, CenterRectOnPoint([0, 0, 15, 15], LL.resx - LL.scalePresentX(x(explodeInd)), LL.resy - LL.scalePresentY(y(explodeInd))));
                end
                
                disp(sprintf('time since launch: %f', getSecs - launchTime))
                disp('Finished simulation. state 0')

                sound(0.1*sin(2 * (1:500)));

                LL.state = 0;
                simRun = false;
                launchTime = NaN;

                waitStart = getSecs;
                waitTime = 1;
                ITIStart = getSecs;
                ITI = 3;
                feedbackInd = tabletCt;
            elseif feedbackMode == 2,
                %%%%%% Delayed path feedback simulation %%%%%%
                if isnan(delayFeedbackStart),
                    delayFeedbackStart = launchTime;
                end
                if getSecs - delayFeedbackStart > 2,
                    delayFeedbackStart = NaN;
                    if feedbackOn,
                        for traj_i = 2:explodeInd,
                            Screen('DrawLine',  LL.w(2), LL.green, ...
                                LL.scalePresentX(x(traj_i - 1)), ...
                                LL.scalePresentY(y(traj_i - 1)), ...
                                LL.scalePresentX(x(traj_i)), ...
                                LL.scalePresentY(y(traj_i)));
                            Screen('DrawLine',  LL.w(2), LL.green, ...
                                LL.resx - LL.scalePresentX(x(traj_i - 1)), ...
                                LL.resy - LL.scalePresentY(y(traj_i - 1)), ...
                                LL.resx - LL.scalePresentX(x(traj_i)), ...
                                LL.resy - LL.scalePresentY(y(traj_i)));
                        end
                    else
                        for traj_i = 2:explodeInd,
                            Screen('DrawLine',  LL.w(2), LL.green, LL.scalePresentX(x(traj_i - 1)), LL.scalePresentY(y(traj_i - 1)), LL.scalePresentX(x(traj_i)), LL.scalePresentY(y(traj_i)));
                        end
                    end
                    Screen('FillOval', LL.w(2), LL.green, CenterRectOnPoint([0, 0, 15, 15], LL.scalePresentX(x(explodeInd)), LL.scalePresentY(y(explodeInd))));
                    if feedbackOn,
                        Screen('FillOval', LL.w(2), LL.green, CenterRectOnPoint([0, 0, 15, 15], LL.resx - LL.scalePresentX(x(explodeInd)), LL.resy - LL.scalePresentY(y(explodeInd))));
                    end
                    
                    disp(sprintf('time since launch: %f', getSecs - launchTime))
                    disp('Finished simulation. state 0')
                    
                    sound(0.1*sin(2 * (1:500)));
                    LL.state = 0;
                    simRun = false;
                    launchTime = NaN;
                    
                    waitStart = getSecs;
                    waitTime = 1;
                    ITIStart = getSecs;
                    ITI = 1;
                    feedbackInd = tabletCt;
                end
            elseif feedbackMode == 3,
                %%%%%% Natural motion simulation %%%%%%
                currTime = getSecs - launchTime; % Time elapsed since launch
                
                % Find the point until which we need to show the trajectory
                t2 = min(find(currTime >= t, 1, 'last'), explodeInd); 
                for traj_i = simCt:t2,
                    if feedbackOn,
                        if traj_i > 1,
                            % Experimenter screen
                            Screen('DrawLine',  LL.w(2), LL.green, LL.scalePresentX(x(traj_i - 1)), LL.scalePresentY(y(traj_i - 1)), LL.scalePresentX(x(traj_i)), LL.scalePresentY(y(traj_i)));
                            
                            % Subject screen
                            Screen('DrawLine',  LL.w(2), LL.green, LL.resx - LL.scalePresentX(x(traj_i - 1)), LL.resy - LL.scalePresentY(y(traj_i - 1)), LL.resx - LL.scalePresentX(x(traj_i)), LL.resy - LL.scalePresentY(y(traj_i)));
                        end
                    else
                        if traj_i > 1,
                            % Experimenter screen
                            Screen('DrawLine',  LL.w(2), LL.green, LL.scalePresentX(x(traj_i - 1)), LL.scalePresentY(y(traj_i - 1)), LL.scalePresentX(x(traj_i)), LL.scalePresentY(y(traj_i)));
                        end
                    end
                end
                
                % Draw current position of ball and refresh buffer
                LL.Recopy(2, 1);
                Screen('FillOval', LL.w(1), LL.green, CenterRectOnPoint([0, 0, 15, 15], LL.scalePresentX(x(t2)), LL.scalePresentY(y(t2))));
                if feedbackOn,
                    Screen('FillOval', LL.w(1), LL.green, CenterRectOnPoint([0, 0, 15, 15], LL.resx - LL.scalePresentX(x(t2)), LL.resy - LL.scalePresentY(y(t2))));
                end
                recopyFlag = 0;
                
                simCt = t2 + 1;
                if simCt > explodeInd,
                    disp(sprintf('time since launch: %f', getSecs - launchTime))
                    disp('Finished simulation. state 0')
                    
                    
                    sound(0.1*sin(2 * (1:500)));
                    LL.state = 0;
                    simRun = false;
                    launchTime = NaN;
                    
                    waitStart = getSecs;
                    feedbackInd = tabletCt;
                    waitTime = 1;
                    ITIStart = getSecs;
                    ITI = 1;
                end
            elseif feedbackMode == 4,
                %%%%%% Instant Reward feedback %%%%%%
                rwdthresh = nanmedian(threshbuffer(26:end));
                 if feedbackGiven == 0,
                     % Only show experimenter
                     for traj_i = 2:explodeInd,
                         Screen('DrawLine',  LL.w(2), LL.green, ...
                             LL.scalePresentX(x(traj_i - 1)), ...
                             LL.scalePresentY(y(traj_i - 1)), ...
                             LL.scalePresentX(x(traj_i)), ...
                             LL.scalePresentY(y(traj_i)));
                     end
                    Screen('FillOval', LL.w(2), LL.green, ...
                        CenterRectOnPoint([0, 0, 15, 15], ...
                        LL.scalePresentX(x(explodeInd)), ...
                        LL.scalePresentY(y(explodeInd))));

                    disp(sprintf('time since launch: %f', getSecs - launchTime))
                    disp('Finished simulation. state 0')
                    
                    if d < rwdthresh,
                        %%%%%% Color target %%%%%%
                        %     Experimenter screen
%                         Screen('FillOval', LL.w(2), LL.green, CenterRectOnPoint([0, 0, 5, 5], LL.scalePresentX(targloc(1)), LL.scalePresentY(targloc(2))));
%                         Screen('FrameOval', LL.w(2), LL.yellow, CenterRectOnPoint([0, 0, 15, 15], LL.scalePresentX(targloc(1)), LL.scalePresentY(targloc(2))));
                        Screen('FrameOval', LL.w(2), LL.green, CenterRectOnPoint([0, 0, 100, 100], LL.scalePresentX(targloc(1)), LL.scalePresentY(targloc(2))));
                        
                        % Subject screen
%                         Screen('FillOval', LL.w(2), LL.green, CenterRectOnPoint([0, 0, 5, 5], LL.resx - LL.scalePresentX(targloc(1)), LL.resy - LL.scalePresentY(targloc(2))));
%                         Screen('FrameOval', LL.w(2), LL.yellow, CenterRectOnPoint([0, 0, 15, 15], LL.resx - LL.scalePresentX(targloc(1)), LL.resy - LL.scalePresentY(targloc(2))));
                        Screen('FrameOval', LL.w(2), LL.green, CenterRectOnPoint([0, 0, 100, 100], LL.resx - LL.scalePresentX(targloc(1)), LL.resy - LL.scalePresentY(targloc(2))));

                        Snd('Play', ting.ting, ting.fs);
                    else
                        sound(0.1*sin(2 * (1:500)));
                    end
                    disp(sprintf('time since launch2: %f', getSecs - launchTime))
                    waitStart = getSecs;
                    waitTime = 1;
                    ITIStart = getSecs;
                    ITI = 3;
                    feedbackInd = tabletCt;
                    launchTime = NaN;
                    feedbackGiven = 1;
                else % if feedback was given already, then detect bad movement endings and beep and indicate
                    if numFails ~= -1 ...
                            && ((~(abs(S.x(samp_i) - LL.launchBarX) <= LL.barwidLaunch / 2 ...
                            && abs(S.y(samp_i) - LL.launchBarY) <= LL.barheiLaunch / 2) ...
                            && LL.v < LL.threshvelStop) ...
                    || (S.x(samp_i) < LL.launchBarX - LL.barwidLaunch / 2) ...
                    || abs(S.y(samp_i) - LL.launchBarY) > LL.barheiLaunch),

                    % Don't go into simRun if wait time is over.
                    elseif LL.v < LL.threshvelStop || getSecs - waitStart > waitTime, 
                        simRun = false;
                        feedbackGiven = 0;
                    end
                end 
            
            elseif feedbackMode == 5,
                %%%%%% Delayed Reward feedback %%%%%%
                if isnan(delayFeedbackStart),
                    delayFeedbackStart = launchTime;
                end
                if getSecs - delayFeedbackStart > 2,
                    delayFeedbackStart = NaN;
                    
                    % Only show feedback to experimenter
                    for traj_i = 2:explodeInd,
                        Screen('DrawLine',  LL.w(2), LL.green, ...
                            LL.scalePresentX(x(traj_i - 1)), ...
                            LL.scalePresentY(y(traj_i - 1)), ...
                            LL.scalePresentX(x(traj_i)), ...
                            LL.scalePresentY(y(traj_i)));
                    end
                    Screen('FillOval', LL.w(2), LL.green, ...
                        CenterRectOnPoint([0, 0, 15, 15], ...
                        LL.scalePresentX(x(explodeInd)), ...
                        LL.scalePresentY(y(explodeInd))));
                    
                    disp(sprintf('time since launch: %f', getSecs - launchTime))
                    disp('Finished simulation. state 0')
                    
                    if d < rwdthresh,
                        %%%%%% Color target %%%%%%
                        %     Experimenter screen
%                         Screen('FillOval', LL.w(2), LL.green, CenterRectOnPoint([0, 0, 5, 5], LL.scalePresentX(targloc(1)), LL.scalePresentY(targloc(2))));
%                         Screen('FrameOval', LL.w(2), LL.yellow, CenterRectOnPoint([0, 0, 15, 15], LL.scalePresentX(targloc(1)), LL.scalePresentY(targloc(2))));
                        Screen('FrameOval', LL.w(2), LL.green, CenterRectOnPoint([0, 0, 100, 100], LL.scalePresentX(targloc(1)), LL.scalePresentY(targloc(2))));
                        
                        % Subject screen
%                         Screen('FillOval', LL.w(2), LL.green, CenterRectOnPoint([0, 0, 5, 5], LL.resx - LL.scalePresentX(targloc(1)), LL.resy - LL.scalePresentY(targloc(2))));
%                         Screen('FrameOval', LL.w(2), LL.yellow, CenterRectOnPoint([0, 0, 15, 15], LL.resx - LL.scalePresentX(targloc(1)), LL.resy - LL.scalePresentY(targloc(2))));
                        Screen('FrameOval', LL.w(2), LL.green, CenterRectOnPoint([0, 0, 100, 100], LL.resx - LL.scalePresentX(targloc(1)), LL.resy - LL.scalePresentY(targloc(2))));

                        Snd('Play', ting.ting, ting.fs);
                    else
                        sound(0.1*sin(2 * (1:500)));
                    end
                    
                    LL.state = 0;
                    simRun = false;
                    launchTime = NaN;
                    
                    waitStart = getSecs;
                    waitTime = 1;
                    ITIStart = getSecs;
                    ITI = 1;
                    feedbackInd = tabletCt;
                end
            end
        else
            recopyFlag = 1;
        end
        
        % Don't let feedback disappear until after fixed delay
        if getSecs - waitStart > waitTime,
            disp(sprintf('Cleared trajectory after: %f', getSecs - waitStart))

            
            LL.state = 0;
            
            % Clear trajectory off screen
            LL.Recopy(3, 2);
            LL.Recopy(2, 1);
            
            % Save trial performance stats
            record(trial, :) = [px, py, vx, vy, targloc(1), targloc(2), ...
                feedbackMode, explodeTime, L, d, releaseTime, releaseInd, ...
                waitStart - init_t, feedbackInd, curr_t, tabletCt, ...
                waitTime, d < rwdthresh, rwdthresh, histwin, numFails, ...
                startTime, startInd];
            csvwrite_with_headers([dirName, '/', subject_name, ...
                '_trialrecord_', char('a' + blocknum - 1), '.csv'], ...
                record(1:trial, :), header);
            
            numFails = 0;
            trial = trial + 1;
            waitStart = NaN;
        end
        
        % Don't let next trial begin until after fixed delay
        if getSecs - ITIStart > ITI,
           ITIStart = NaN; 
        end
        LL.Redraw(S, recopyFlag);
        [~,~, mouse] = GetMouse;
    end
    
    % Save all tablet recordings
    tabletRec = tabletRec(1:tabletCt, :);
    tabletRecHeader = {'x', 'y', 'tabletTime', 'computerTime'};
    save([dirName, '/', subject_name, '_tabletRec_', ...
        char('a' + blocknum - 1)], 'tabletRec', 'tabletRecHeader')    

    LL.Finish;
    S.Finish;
 
%%
catch ME
    LL.Finish;
    S.Finish;
    ME.message
    ME.stack.line
    keyboard;
end %catch
end



