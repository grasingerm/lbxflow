data0b = csvread('data/poise/chen/explicit/000000/ux_profile.dsv');
data1b = csvread('data/poise/chen/explicit/000004/ux_profile.dsv');
data2b = csvread('data/poise/chen/explicit/000008/ux_profile.dsv');
data3b = csvread('data/poise/chen/explicit/0000012/ux_profile.dsv');
x = data1b(:,1);

figure();
plot(x,data0b(:,2),'r-<',x,data1b(:,2),'b->',x,data2b(:,2),'g-^',x,data3b(:,2),'y-O');
legend('\tau_y = 0.00000','\tau_y = 0.00004','\tau_y = 0.00008','\tau_y = 0.00012');
xlabel('y / H');
ylabel('u / u_{max}');

data0c = csvread('data/poise/chen/implicit/000000/ux_profile.dsv');
data1c = csvread('data/poise/chen/implicit/000004/ux_profile.dsv');
data2c = csvread('data/poise/chen/implicit/000008/ux_profile.dsv');
data3c = csvread('data/poise/chen/implicit/0000012/ux_profile.dsv');
x = data1c(:,1);

figure();
plot(x,data0c(:,2),'r-<',x,data1c(:,2),'b->',x,data2c(:,2),'g-^',x,data3c(:,2),'y-O');
legend('\tau_y = 0.00000','\tau_y = 0.00004','\tau_y = 0.00008','\tau_y = 0.00012');
xlabel('y / H');
ylabel('u / u_{max}');

data0d = csvread('data/poise/chen/explicit/000000/ux_profile.dsv');
data1d = csvread('data/poise/chen/implicit/000004/ux_profile.dsv');
data2d = csvread('data/poise/chen/implicit/000008/ux_profile.dsv');
data3d = csvread('data/poise/chen/implicit/0000012/ux_profile.dsv');
x = data1d(:,1);

figure();
plot(x,data0d(:,2),'r-<',x,data1d(:,2),'b->',x,data2d(:,2),'g-^',x,data3d(:,2),'y-O');
legend('\tau_y = 0.00000','\tau_y = 0.00004','\tau_y = 0.00008','\tau_y = 0.00012');
xlabel('y / H');
ylabel('u / u_{max}');