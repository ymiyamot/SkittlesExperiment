function SkittlesGame(subject_name)
    mkdir(subject_name); 
    totalblocks = 30;
    baseblocks = 2;
    blockLen = 50;
    
    for bl = 1:baseblocks,
        fprintf('Block %d: baseline', bl);
        keyboard;
        Tester(subject_name, bl, blockLen, 1);
    end
    
    for bl=baseblocks + 1:totalblocks,
%     for bl=10:totalblocks,
        fprintf('Block %d: training', bl);
        keyboard;
        Tester(subject_name, bl, blockLen, 0);
    end
end

function Tester(subject_name, blocknum, total_trials, baseMode)
try
    % Initialization 
    init_t = getSecs();
    
    % Variable for recording data
    tabletRec = NaN(1000000, 4);
    tabletCt = 0;
    
    % Conversion: 
%     nzlevels = [0, 0];
%     [targloc0, targloc] = deal([1200, 270]);
%     [targloc0, targloc] = deal(targloc);
        

    % Simulation parameters
    feedbackMode = 1;

    if baseMode == 1,
        targloc = [NaN NaN];
        feedbackOn = 0;
    else
        targloc = [1100 + 75 * sqrt(3), 270];
        feedbackOn = 1;
    end
    
    LL = FreeLaunch(960 + 150, 540 + 250, 960 + 50 - 185.71, 540 + 250, targloc);
%     FB = FBSimulator(targloc);
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
    


%     handpath = NaN(10000, 4);
%     all_handpath = cell(total_trials, 1);
%     recordNew = 1;
    
    % retrieve history for history-dependent reward
    histwin = 25;
    if blocknum == 3 || baseMode == 1, % Assuming two baseline blocks
        threshbuffer = NaN(50, 1);
        threshbuffer(50) = 200;
    else
        rec = csvimport([subject_name, '/', subject_name, '_trialrecord_', char('a' + blocknum - 2), '.csv']);
        tmp = rec(2:end, 10);
        recVec = NaN(size(tmp));
        for tmp_i = 1:length(tmp),
            if ~isnan(tmp{tmp_i}) && all('string' == class(tmp{tmp_i}))
                recVec(tmp_i) = str2num(tmp{tmp_i});
            elseif ~isnan(tmp{tmp_i}) && all('double' == class(tmp{tmp_i}))
                recVec(tmp_i) = tmp{tmp_i};
            end
        end
        threshbuffer = recVec;
    end

    % Storage for trial performance information
    record = NaN(total_trials, 21);
    header = {'release_x', 'release_y', 'release_vx', 'release_vy', 'targ_x', 'targ_y', 'feedbackMode', 'explodeTime', 'pendulumLength', 'distFromTarg', 'releaseTime', 'releaseInd', 'feedbackStartTime', 'feedbackStartInd', 'feedbackEndTime', 'feedbackEndInd', 'waitTime', 'rwd', 'rwdThresh', 'historyWindowLength', 'numFails'};
    
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
            elseif LL.state == 1 && (S.x(samp_i) < LL.midBarX && (abs(S.y(samp_i) - LL.startBarY) <= LL.barheiLaunch / 2))
                disp('Launched ball. state 2')
                
                LL.state = 2;
                launchTime = getSecs;
                simRun = true;
                simTrigger = true;

            elseif LL.state == 2 && ((~(abs(S.x(samp_i) - LL.launchBarX) <= LL.barwidLaunch / 2 && abs(S.y(samp_i) - LL.launchBarY) <= LL.barheiLaunch / 2) && LL.v < LL.threshvelStop) ...
                    || (S.x(samp_i) < LL.launchBarX - LL.barwidLaunch / 2) ...
                    || abs(S.y(samp_i) - LL.launchBarY) > LL.barheiLaunch),
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
            
            % Time of release
            releaseTime = curr_t;
            releaseInd = tabletCt;

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
                    ITI = 1;
                    feedbackInd = tabletCt;
                end
            elseif feedbackMode == 3,
                %%%%%% Natural motion simulation %%%%%%
                currTime = getSecs - launchTime; % Time elapsed since launch
                
                t2 = min(find(currTime >= t, 1, 'last'), explodeInd); % Find the point until which we need to show the trajectory
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
            record(trial, :) = [px, py, vx, vy, targloc(1), targloc(2), feedbackMode, explodeTime, L, d, releaseTime, releaseInd, waitStart - init_t, feedbackInd, curr_t, tabletCt, waitTime, d < rwdthresh, rwdthresh, histwin, numFails];
            csvwrite_with_headers([subject_name, '/', subject_name, '_trialrecord_', char('a' + blocknum - 1), '.csv'], record(1:trial, :), header);
            
            numFails = 0;
            trial = trial + 1;
            waitStart = NaN;
        end
        
        % Don't let next trial begin until after fixed delay
        if getSecs - ITIStart > ITI,
           ITIStart = NaN; 
        end

%         [~, ~, keyCode, ~] = kbCheck();

%         if LL.state == 2
%             trial = trial + 1;
%             vx = diff(LL.buffer(end - 1:end, 1)) / diff(LL.buffer(end - 1:end, 3)); % pixel units
%             vy = -diff(LL.buffer(end - 1:end, 2)) / diff(LL.buffer(end - 1:end, 3)); % pixel units
%             xr = LL.buffer(end, 1);
%             yr = LL.buffer(end, 2);
%             
% %             speed = 1;
% %             wait_time = 0.01;
%             drawscore = 0;
%             L = 12000;
%             explode_time = 2;
%             
%             rwdthresh = nanmedian(threshbuffer(26:end));
%             [x, y, d, score] = FB.SkittlesSimulation(xr, yr, vx, vy, speed, L, explode_time, wait_time, drawscore, rwdthresh, feedback_on);
%             
%             % Update history
%             threshbuffer = [threshbuffer(2:end); d];
%             
%             LL.state = 3;
%             record(trial, :) = [xr, yr, vx, vy, targloc(1), targloc(2), speed, explode_time, L, d, score, curr_t, wait_time, d < rwdthresh, rwdthresh, histwin];
%             
%             csvwrite_with_headers([subject_name, '/', subject_name, '_trialrecord_', char('a' + blocknum - 1), '.csv'], record(1:trial, :), header);
% 
% %             targloc = targloc0 + (rand(1, 2) - 0.5) * diag(nzlevels);
%             targloc = targloc0 + randn(1, 2) * diag(nzlevels);
%             FB.targloc = targloc;
%             FB.Draw_Target;
%             
%         end
%         
%         if LL.state == 3,
%             % Record handpath
%             disp(sprintf('State 3. Recording path for trial %d', trial));
%             all_handpath{trial} = handpath(1:sct, :);
%             
%             recordNew = 1;
%         end
        
        LL.Redraw(S, recopyFlag);
        [~,~, mouse] = GetMouse;
    end
    
    % Save all tablet recordings
    tabletRec = tabletRec(1:tabletCt, :);
    tabletRecHeader = {'x', 'y', 'tabletTime', 'computerTime'};
    save([subject_name, '/', subject_name, '_tabletRec_', char('a' + blocknum - 1)], 'tabletRec', 'tabletRecHeader')    
%     all_handpath{trial} = handpath(1:sct, :);
%     save([subject_name, '/', subject_name, '_handpaths_', char('a' + blocknum - 1)], 'all_handpath')
    
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



