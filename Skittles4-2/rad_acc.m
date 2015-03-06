function acc = rad_acc(pos, L)
g = -9.8 * 100 / 2.54 * 100; % m/sec^2
% keyboard;
r = sqrt(pos(1)^2 + pos(2)^2);
if r==0,
    a_i = [0, 0, 0];
else
    a_i = [-pos(1) * pos(3) / r, -pos(2) * pos(3) / r, r] / L;
end
a_t = g * r / L;
acc = a_t * a_i;

end