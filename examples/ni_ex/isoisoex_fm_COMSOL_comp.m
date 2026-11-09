clear; clc; close all

T = 500;

true_params = load("true_params.mat");
lnkf = reshape(log(polyval(true_params.gold_film.k, T)),1,1,[]);
lnCf = reshape(log(polyval(true_params.gold_film.c, T)) ...
     + log(polyval(true_params.gold_film.rho, T)), 1,1,[]);
lnaf = log(true_params.gold_film.alpha);
lnhf = log(true_params.gold_film.h);

lnks = log(polyval(true_params.YSZ.k, T));
lnCs = log(polyval(true_params.YSZ.c, T)) ...
         + log(polyval(true_params.YSZ.rho, T));
lnas = log(true_params.YSZ.alpha);

logitR = -inf;  % R = 0 (phase is independent of R)
lnRth = -inf;
sx = 2;
sy = 2;
P = 1000;
f = reshape(10.^(-5:-1),1,1,1,1,[]);

Df = exp(get_val(lnkf)-get_val(lnCf)); % mm^2/s
Ds = exp(get_val(lnks)-get_val(lnCs)); % mm^2/s
Lthf = sqrt(reshape(Df,1,1,[]) ./ pi ./ reshape(f,1,1,1,1,[])); % um
Lths = sqrt(reshape(Ds,1,1,[]) ./ pi ./ reshape(f,1,1,1,1,[])); % um
x_max = max(100, 2*Lths);

Nx = 160;

xprobe = linspace(0,30,101);
Xprobe = [xprobe(:), xprobe(:)];

fm = ForwardModel("iso", "iso", true);

T0tilde = fm.solve( ...
    {lnkf}, {}, lnCf, lnaf, logitR, lnhf, ...
    {lnks}, {}, lnCs, lnas, logitR, [], ...
    lnRth, sx^2, sy^2, 0, P, f, x_max, Nx, Xprobe);

colors = colororder;
COMSOL = load("isoisoex_COMSOL_output.mat");

tiledlayout(2,1,"TileSpacing","tight","Padding","tight");
nexttile; hold on;
plot(COMSOL.r, unwrap(angle(squeeze(COMSOL.T0tilde(:,3,1,1:4:end)))), Color=colors(1,:), LineWidth=1.5)
plot(sqrt(Xprobe(:,1).^2+Xprobe(:,2).^2), unwrap(squeeze(angle(T0tilde))), "--", Color=colors(2,:), LineWidth=1.5)
xlim([0,30])
xline([5,10,20])
ylabel("Phase Lag [rad]", Interpreter="latex")
text(21, 0.04, "$f=10$ Hz", Interpreter="latex", Rotation=-1, FontSize=8)
text(21, -0.55, "$f=100$ Hz", Interpreter="latex", Rotation=-4, FontSize=8)
text(21, -1.06, "$f=1$ kHz", Interpreter="latex", Rotation=-12, FontSize=8)
text(21, -2.04, "$f=10$ kHz", Interpreter="latex", Rotation=-22, FontSize=8)
text(21, -3.6, "$f=100$ kHz", Interpreter="latex", Rotation=-35, FontSize=8)
ylim([-5,0.3])

nexttile; hold on;
plt1 = plot(COMSOL.r, unwrap(abs(squeeze(COMSOL.T0tilde(:,3,1,1:4:end)))), Color=colors(1,:));
plt2 = plot(sqrt(Xprobe(:,1).^2+Xprobe(:,2).^2), unwrap(squeeze(2*abs(T0tilde))), "--", Color=colors(2,:));
xlim([0,30])
yscale log
xline([5,10,20])

lgd = legend([plt1(1), plt2(1)], ["COMSOL", "Forward Model"], Interpreter="latex", Orientation="horizontal");
lgd.Layout.Tile = 'north';
ylabel("Amplitude [K]", Interpreter="latex")
xlabel("$\left\|  [x,y] \right\|$ [micron]", Interpreter="latex")
text(21, 2.4, "$f=10$ Hz", Interpreter="latex", Rotation=-5, FontSize=8)
text(21, 0.75, "$f=100$ Hz", Interpreter="latex", Rotation=-7, FontSize=8)
text(21, .3, "$f=1$ kHz", Interpreter="latex", Rotation=-8.5, FontSize=8)
text(21, 4.5e-2, "$f=10$ kHz", Interpreter="latex", Rotation=-19, FontSize=8)
text(21, 1.25e-3, "$f=100$ kHz", Interpreter="latex", Rotation=-32.5, FontSize=8)
ylim([10^-4,1e1])

fig = gcf;
fig.Units = 'inches';
fig.Position = [1 1 4 6];