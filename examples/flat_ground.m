clear;
args.freq_hz = 5900e6;
args.d_x = 5;
args.x_max_m = 700;
args.n_y = 1024*3;
args.n_z = 1024*4;
args.y_output_filter = 3;
args.z_output_filter = 4;
args.dy_wl = 0.4;
args.dz_wl = 0.4;
args.z_output_filter = 1;
args.antenna_height_m = 4;

res_ssf = SSF3D(args);

[~, zo1_ssf] = min(abs(res_ssf.z_grid_m-1.5));
[~, yo1_ssf] = min(abs(res_ssf.y_grid_m-0));

figure;imagesc(res_ssf.x_grid_m-args.d_x/2, res_ssf.z_grid_m(end:-1:1), ...
    squeeze(20*log10(abs(res_ssf.field(:, yo1_ssf, end:-1:1)))).'-38, [-120 -30])
colormap(jet)
set(gca,'YDir','normal')
xlabel('Range (m)');
ylabel('Height (m)');
ylim([0 30])
grid on
set(gcf, 'Position', [0 0 600 300]);

figure;imagesc(res_ssf.x_grid_m, res_ssf.y_grid_m, ...
    squeeze(20*log10(abs(res_ssf.field(:, :, zo1_ssf)))).'-38, [-120 -30])
colormap(jet)
xlabel('Range (m)');
ylabel('y (m)');
ylim([-20 20])
grid on
set(gcf, 'Position', [0 0 600 300]);

figure;plot(res_ssf.x_grid_m, squeeze(20*log10(abs(res_ssf.field(:, yo1_ssf, zo1_ssf))))-38);
line([0 700],[-96 -96],'linestyle','--');
xlabel('Range (m)');
ylabel('Power (dBm)');
legend('3D parabolic equation');
grid on
set(gcf, 'Position', [0 0 600 300]);