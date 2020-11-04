clear;
args.freq_hz = 5900e6;
args.d_x = 2;
args.x_max_m = 300;
args.n_y = 1024*3;
args.n_z = 1024*4;
args.y_output_filter = 3;
args.z_output_filter = 4;
args.dy_wl = 1 / (2*sqrt(2));
args.dz_wl = 1 / (2*sqrt(2));
args.antenna_height_m = 4;
brick.x_min_m = 100;
brick.x_max_m = 110;
brick.y_min_m = -10;
brick.y_max_m = 10;
brick.height_m = 20;
args.bricks = [brick];
args.two_way = true;

res_ssf = SSF3D(args);

[~, zo1_ssf] = min(abs(res_ssf.z_grid_m-1.5));
[~, yo1_ssf] = min(abs(res_ssf.y_grid_m-0));

figure;imagesc(res_ssf.x_grid_m+args.d_x/2, res_ssf.z_grid_m(end:-1:1), ...
    squeeze(20*log10(abs(res_ssf.field(:, yo1_ssf, end:-1:1)))).'-38, [-120 -30])
set(gca,'YDir','normal')
ylim([0 30])
colormap(jet)
hold on;
rectangle('Position', [100 0 10 20], 'FaceColor', uint8([220 220 220]));
grid on
xlabel('Range (m)');
ylabel('Height (m)');
colorbar;
set(gcf, 'Position', [0 0 600 300]);

figure;imagesc(res_ssf.x_grid_m, res_ssf.y_grid_m, ...
    squeeze(20*log10(abs(res_ssf.field(:, :, zo1_ssf)))).'-38, [-120 -30])
ylim([-20 20])
colormap(jet)
hold on;
rectangle('Position', [100 -10 10 20], 'FaceColor', uint8([220 220 220]));
grid on
xlabel('Range (m)');
ylabel('y (m)');
colorbar;
set(gcf, 'Position', [0 0 600 300]);

figure;plot(res_ssf.x_grid_m, squeeze(20*log10(abs(res_ssf.field(:, yo1_ssf, zo1_ssf))))-38);
line([0 300],[-96 -96],'linestyle','--');
ylim([-120 -30])
xlabel('Range (m)');
ylabel('Power (dBm)');
legend('3D parabolic equation');
grid on
set(gcf, 'Position', [0 0 600 300]);
