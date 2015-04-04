% invoke example program and store output in eekf.log
system('../../build/examples/eekf_example > ../../build/exampleseekf.log');

% read eekf.log into D
D = dlmread('../../build/exampleseekf.log', ' ', 1, 0);

% plot data

subplot(2,1,1)
plot(D(:,[2 8]))
xlabel('k')
ylabel('x')
legend x rx;

subplot(2,1,2)
plot(D(:,[3 9 10]))
xlabel('k')
ylabel('dx')
legend dx rdx z;


