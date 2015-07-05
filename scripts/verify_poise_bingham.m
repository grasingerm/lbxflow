data0b = csvread('data/poise/chen/explicit/000000/ux_profile.dsv');
data1b = csvread('data/poise/chen/explicit/000004/ux_profile.dsv');
data2b = csvread('data/poise/chen/explicit/000008/ux_profile.dsv');
data3b = csvread('data/poise/chen/explicit/000012/ux_profile.dsv');
x = data1b(:,1);

nj = size(x,1);
h = (nj-1)/2.0;

anal = zeros(nj,4);
pgrad = -5.6e-6;
mu = 0.2;
tauys = [0.0; 4e-5; 8e-5; 12e-5];
xs = (1:nj) - nj/2 - 0.5;

for j=1:4
    tau = tauys(j);
    y_tau = -tau / pgrad;
    for i=1:size(x,1)
        xi = xs(i);
        if abs(xi) <= y_tau
            anal(i,j) = -1.0 / (2.0 * mu) * pgrad * (h^2 - y_tau^2) - tau / mu * (h - y_tau);
        else
            anal(i,j) = -1.0 / (2.0 * mu) * pgrad * (h^2 - xi^2) - tau / mu * (h - abs(xi));
        end
    end
end
            

figure();
plot(x,data0b(:,2),'r-',x,anal(:,1),'r<',x,data1b(:,2),'b--',x,anal(:,2),'b>',x,data2b(:,2),'g:',x,anal(:,3),'g^',x,data3b(:,2),'y-.',x,anal(:,4),'yO');
legend('lbm, \tau_y = 0.00000','analytical, \tau_y = 0.00000','lbm, \tau_y = 0.00004','analytical, \tau_y = 0.00004','lbm, \tau_y = 0.00008','analytical, \tau_y = 0.00008','lbm, \tau_y = 0.00012','analytical, \tau_y = 0.00012');
xlabel('y / H');
ylabel('u');

figure();
plot(x,data0b(:,2)/max(data0b(:,2)),'r-<',x,data1b(:,2)/max(data1b(:,2)),'b->',x,data2b(:,2)/max(data2b(:,2)),'g-^',x,data3b(:,2)/max(data3b(:,2)),'y-O');
legend('\tau_y = 0.00000','\tau_y = 0.00004','\tau_y = 0.00008','\tau_y = 0.00012');
xlabel('y / H');
ylabel('u / u_{max}');