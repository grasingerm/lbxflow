prof100 = csvread('../data/prof-at-100_step-10000.dsv');
prof250 = csvread('../data/prof-at-250_step-10000.dsv');
prof500 = csvread('../data/prof-at-500_step-10000.dsv');

% analytical solution

plot(prof100(:,1), prof100(:,2), 'r-', ... 
    prof250(:,1), prof250(:,2), 'b--', ...
    prof500(:,1), prof500(:,2), 'g.-');
legend('x=100','x=250','x=500');