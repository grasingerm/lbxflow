data1 = csvread('data\sukop\ux_profile.dsv');
data1(:,1) = (data1(:,1) - 6.5) / 11.0;

grad_p = -1.102e-3
nu = 1/6
rhoo = 1.0
mu = nu * rhoo
h = 5.5

y = data1(:,1)
u = -1 / (2 * mu) * grad_p * (h^2 - (2*h*y).^2)
e = zeros(length(u), 1);

for i=1:length(u)
    if (u(i) == 0)
        e(i) = abs(u(i) - data1(i,2));
    else
        e(i) = abs(u(i) - data1(i,2)) / abs(u(i));
    end
end

sum(e)

figure(1);
plot(y, data1(:,2), 'r-', y, u, 'bo');
legend('LBM-Gou', 'analytical');
ylabel('u (lu / s)');
xlabel('y (lu / lu)');

figure(2);
plot(y, e);
title('Relative Error');