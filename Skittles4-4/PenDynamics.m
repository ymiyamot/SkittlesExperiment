% We will use leapfrog numerical integration
function allpos = PenDynamics(pos0, vel0, dt, total_t, L)

numsteps = total_t / dt;

pos = pos0;
vel = vel0;

% Initial half step for velocity
vel = vel + (dt / 2) * rad_acc(pos, L);
allpos = NaN(numsteps, 3);
allpos(1, :) = pos0;
for i = 1:numsteps - 1,
    % Full step for position
    pos = pos + dt * vel;
    pos = pos / sqrt(sum(pos.^2)) * L; % normalize, to keep position on sphere
    allpos(i + 1, :) = pos;
    vel = diff(allpos(i:i + 1, :)) / dt; % correct velocity
    
    % Full step for velocity except last step
    if i ~= numsteps - 1,
        vel = vel + dt * rad_acc(pos, L);
    else
        vel = vel + (dt / 2) * rad_acc(pos, L);
    end
end


% figure; 
% cols = cool(size(allpos, 1));
% scatter3(allpos(:, 1), allpos(:, 2), allpos(:, 3), 4, cols(1:end, :), 'fill');
% view([0.5, 50])
