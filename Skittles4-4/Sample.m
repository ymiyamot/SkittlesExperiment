classdef Sample < handle
    properties
        x
        y
        t
        samp 
        usetablet
    end
    properties (SetAccess = private)
        resx
        resy
        startTime
    end
    methods
        function S = Sample(wPtr, resx, resy, usetablet, startTime)
            % If mouse==1, then we're using a tablet.
            S.resx = resx;
            S.resy = resy;
            S.startTime = startTime;
            S.samp = 0;
            S.usetablet = usetablet;
            disp('Initializing Sampler...')
            if S.usetablet==1
                WinTabMex(0, wPtr); %Initialize tablet driver, connect it to 'wPtr'  
                WinTabMex(2); %Empties the packet queue in preparation for collecting actual data
            end
            S.x = NaN;
            S.y = NaN;
            S.t = NaN;
            disp('Finished Initializing Sampler')
        end
        function S = Update(S)
            if S.usetablet, % Are we using the tablet?
                pkt = WinTabMex(5);
                if isempty(pkt) % Does the tablet have a sample?
                    S.samp = 0; 
                    return; % No, so exit
                else
                    S.samp = 1; % We have a sample
                end
                
                % Collect current sample(s)
                [S.x S.y S.t] = deal(NaN(20,1));
                ct = 0;
                while(~isempty(pkt))
                    ct = ct+1;
                    [S.x(ct) S.y(ct)] = S.count2pixel(pkt(1),pkt(2));
                    S.t(ct) = pkt(6)/1000 - S.startTime; % Time of measurement in sec
                    pkt = WinTabMex(5); % Look at next packet.
                end
                S.x = S.x(1:ct);
                S.y = S.y(1:ct);
                S.t = S.t(1:ct);
            else % If we're using computer mouse
                [S.x S.y] = GetMouse;
                S.t = GetSecs - S.startTime;
                S.samp = 1; % We have a sample
            end
        end
        function [x y] = count2pixel(S, c_x, c_y)
            xy = [c_x c_y] / 25.4; % Change to pixels
            xy(1) = -xy(1) + 1920; % Flip x
            if(xy(2)>1080) % Tablet may be larger than screen
                xy(2) = 1080;
            end
            x = xy(1);
            y = xy(2);
        end
        function S = Finish(S)
%             WinTabMex(0)
            disp('Sampler Finish');
            if S.usetablet,
                WinTabMex(3); % Stop/Pause data acquisition.
                WinTabMex(1); % Shutdown driver.
            end
        end
    end %methods
end %classdef
            