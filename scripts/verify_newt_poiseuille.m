data = csvread('data\newt-srt_ss\prof-midchan_step-5000.dsv');
data(:,1) = (data(:,1) - 13) / 24.0;

p_in = 0.9 / 3.
p_out = 0.885 / 3.
delta_p = p_out - p_in
ni = 250
grad_p = delta_p / ni
nu = 0.02
rhoo = 1.0
mu = nu * rhoo
h = 12

y = data(:,1)
u = -h^2 / (2 * mu) * grad_p * (1 - (24 * y/h).^2)

u_bar = u / max(u);
lbm_u_bar = data(:,2) / max(data(:,2))

plot(y, data(:,2), 'r-', y, u, 'bo');
legend('LBM','analytical');

