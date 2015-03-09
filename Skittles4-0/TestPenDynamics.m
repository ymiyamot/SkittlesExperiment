sfac =  100 / 2.54 * 100;
L = 0.13 * sfac;
total_t = 0.1;
dt = 0.01;
pos0 = [0, -L, 0];
pos0 = pos0 / sqrt(sum(pos0.^2));
vel0 = [-0.3 * sfac, 0, 0];

pos = PenDynamics(pos0, vel0, dt, total_t, L);



L = 1;
total_t = 100;
dt = 0.01;
pos0 = [-L/10, 0, 0];
pos0 = pos0 / sqrt(sum(pos0.^2));
vel0 = [0, 0, 0];

pos = PenDynamics(pos0, vel0, dt, total_t, L);
figure; plot(pos)

x = diff(pos(:, 1));
dt * diff(find(x(1:end-1) > 0 & x(2:end) < 0))

z = diff(pos(:, 3));
diff(find(z(1:end-1) > 0 & z(2:end) < 0))