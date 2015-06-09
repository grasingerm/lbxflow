data0a = csvread('data/newt_006_000000/ubar_profile.dsv');
data1a = csvread('data/bmrt_006_000000/ubar_profile.dsv');
data2a = csvread('data/bmrt_006_000004/ubar_profile.dsv');
data3a = csvread('data/bmrt_006_000008/ubar_profile.dsv');
data4a = csvread('data/bmrt_006_000012/ubar_profile.dsv');

data0b = csvread('data/newt_006_000000/ux_profile.dsv');
data1b = csvread('data/bmrt_006_000000/ux_profile.dsv');
data2b = csvread('data/bmrt_006_000004/ux_profile.dsv');
data3b = csvread('data/bmrt_006_000008/ux_profile.dsv');
data4b = csvread('data/bmrt_006_000012/ux_profile.dsv');

x = data1a(:,1);

figure();
plot(x,data1a(:,2),'r-',x,data0a(:,2),'b--');
legend('\tau_y = 0.00000','Newtonian');
xlabel('y / H');
ylabel('u / u_{max}');

figure();
plot(x,data1b(:,2),'r-',x,data0b(:,2),'b--');
legend('\tau_y = 0.00000','Newtonian');
xlabel('y / H');
ylabel('u (lat / s)');

figure();
plot(x,data1a(:,2),'r-<',x,data2a(:,2),'b->',x,data3a(:,2),'g-^',x,data4a(:,2),'y-O');
legend('\tau_y = 0.00000','\tau_y = 0.00004','\tau_y = 0.00008','\tau_y = 0.00012');
xlabel('y / H');
ylabel('u / u_{max}');

figure();
plot(x,data1b(:,2),'r<',x,data2b(:,2),'b>',x,data3b(:,2),'g^',x,data4b(:,2),'y-O');
legend('\tau_y = 0.00000','\tau_y = 0.00004','\tau_y = 0.00008','\tau_y = 0.00012');
xlabel('y / H');
ylabel('u (lat / s)');