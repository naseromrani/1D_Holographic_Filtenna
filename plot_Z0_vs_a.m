freq = 16e9;
d = 2.8;
omega = 2 * pi * freq;
epsilon_0 = 8.85e-12;

a = 0.2:0.2:2.4;
phaseX = [59.023, 59.052, 59.117, 59.253, 59.471, 59.823, 60.38, 61.23, 62.522, 64.521, 67.428, 72.061];

betaX = phaseX / d;

Z0 = betaX / (omega * epsilon_0);

plot(a, Z0);
