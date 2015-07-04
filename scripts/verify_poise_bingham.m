data0b = csvread('data/poise/chen/explicit/000000/ux_profile.dsv');
data1b = csvread('data/poise/chen/explicit/000004/ux_profile.dsv');
data2b = csvread('data/poise/chen/explicit/000008/ux_profile.dsv');
data3b = csvread('data/poise/chen/explicit/000012/ux_profile.dsv');
x = data1b(:,1);

figure();
plot(x,data0b(:,2),'r-<',x,data1b(:,2),'b->',x,data2b(:,2),'g-^',x,data3b(:,2),'y-O');
legend('\tau_y = 0.00000','\tau_y = 0.00004','\tau_y = 0.00008','\tau_y = 0.00012');
xlabel('y / H');
ylabel('u / u_{max}');

data0d = csvread('data/poise/wang/explicit/0000/ux_profile.dsv');
data1d = csvread('data/poise/wang/explicit/0001/ux_profile.dsv');
data2d = csvread('data/poise/wang/explicit/0005/ux_profile.dsv');
data3d = csvread('data/poise/wang/explicit/0010/ux_profile.dsv');
x = data1d(:,1);

figure();
plot(x,data0d(:,2),'r-<',x,data1d(:,2),'b->',x,data2d(:,2),'g-^',x,data3d(:,2),'y-O');
legend('\tau_y = 0.000','\tau_y = 0.001','\tau_y = 0.005','\tau_y = 0.010');
xlabel('y / H');
ylabel('u');