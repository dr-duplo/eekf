% invoke example program and store output in eekf.log
system('../../build/examples/eekf_example > ../../build/examples/eekf_example.log');

% read eekf.log into D
D = dlmread('../../build/examples/eekf_example.log', ' ', 1, 0);

% plot data
subplot(3,1,1)
plot(D(:,[2 8 10]))
xlabel('k')
ylabel('position')
legend x p z;

subplot(3,1,2)
plot(D(:,[3 9]))
xlabel('k')
ylabel('velocity')
legend dx v;

subplot(3,1,3)
plot(D(:,10) - D(:,8), 'c'); hold on;
plot(D(:,2) - D(:,8), 'r'); hold off;
xlabel('k')
ylabel('error')
legend 'z - p' 'x - p';

