function [res] = SSF3D(args)
%SSF3D 3D two-way split-step Fourier Parabolic equation
% Input:
%   args.freq_hz - operating frequency (Hz) (obligatory)
%   args.beam_width - Gauss antenna 3dB beam width (deg) (default is 15)
%   args.d_x - longitudal grid step (m) (obligatory)
%   args.dy_wl - y grid step (wavelength) (default is 1)
%   args.dz_wl - z grid step (wavelength) (default is 1)
%   args.x_max_m - max propagation range (m) (obligatory)
%   args.n_y - number of grid points by y (obligatory)
%   args.n_z - number of grid points by z (obligatory)
%   args.antenna_height_m - transmitting antenna height (m) (obligatory)
%   args.bricks array of PEC cuboids
%       brick.x_min_m (obligatory)
%       brick.x_max_m (obligatory)
%       brick.y_min_m (obligatory)
%       brick.y_max_m (obligatory)
%       brick.height_m (obligatory)
%   args.two_way - use two-way PE (default is false)
%   args.n_iters - number of iterations for two-way PE (default is 1)
%   args.x_output_filter - save only each x_output_filter value by x (default is 1)
%   args.y_output_filter - save only each y_output_filter value by y (default is 1)
%   args.z_output_filter - save only each z_output_filter value by z (default is 1)
%
% Output:
% res.x_grid_m - output grid by x (m)
% res.y_grid_m - output grid by y (m); 
% res.z_grid_m - output grid by z (m);
% res.filed - output complex field 3d matrix (x,y,z)


f = args.freq_hz;
dx = args.d_x;
max_range_m = args.x_max_m;
n_y = args.n_y;
n_z = args.n_z;
ant_height = args.antenna_height_m;
if isfield(args, 'bricks')
    bricks = args.bricks;
    for ind=1:length(bricks)
        bricks(ind).front_faset_u = [];
        bricks(ind).front_faset_lo = [];
        bricks(ind).back_faset_u = [];
        bricks(ind).back_faset_lo = [];
    end
else
    bricks = [];
end
if isfield(args, 'beam_width')
    beam_width = args.beam_width;
else
    beam_width = 15;
end
if isfield(args, 'two_way')
    two_way = args.two_way;
else
    two_way = false;
end
if isfield(args, 'n_iters')
    n_iters = args.n_iters;
else
    n_iters = 1;
end
if isfield(args, 'x_output_filter')
    x_output_filter = args.x_output_filter;
else
    x_output_filter = 1;
end
if isfield(args, 'y_output_filter')
    y_output_filter = args.y_output_filter;
else
    y_output_filter = 1;
end
if isfield(args, 'z_output_filter')
    z_output_filter = args.z_output_filter;
else
    z_output_filter = 1;
end
if isfield(args, 'dy_wl')
    dy_wl = args.dy_wl;
else
    dy_wl = 1;
end
if isfield(args, 'dz_wl')
    dz_wl = args.dz_wl;
else
    dz_wl = 1;
end
c = 3e8;
wl = c/f;
k0 = 2*pi/wl;
dz_m = dz_wl*wl;
dy_m = dy_wl*wl;
sigmaz = 1/(pi*deg2rad(beam_width)) / 1.2159 * wl;

max_angle = rad2deg(asin(wl / (2*sqrt(dy_m^2+dz_m^2))));
fprintf('max angle = %f\n', max_angle);

y_max = n_y*dy_m/2;
y1 = 0:dy_m:y_max;
y = [y1 zeros(1,n_y/2-1)];
z_max = n_z*dz_m/2;
z1 = 0:dz_m:z_max;
z = [z1 zeros(1,n_z/2-1)];

ky_max = pi/dy_m;
Dky = 2*ky_max/n_y;
ky1 = -ky_max:Dky:ky_max;
ky = [ky1(:,(n_y/2)+1:n_y+1) ky1(:,2:n_y/2)];
kz_max = pi/dz_m;
Dkz = 2*kz_max/n_z;
kz1 = -kz_max:Dkz:kz_max;
kz = [kz1(:,(n_z/2)+1:n_z+1) kz1(:,2:n_z/2)];

Ky = meshgrid(ky,1:n_z);
Kz = meshgrid(kz,1:n_y).';
kx = sqrt(k0^2-Ky.^2-Kz.^2);
clear Ky;
clear Kz;

Hya=[];
for t=0:n_y/2
    if (t >= 0 & t <= 3*n_y/8)
        hy=1;
    elseif (t >= 3*n_y/8 & t <= n_y/2)
        hy=(sin(4*pi*t/n_y))^2;
    end
    Hya=[Hya hy];
end
Yy = fliplr(Hya(:,2:n_y/2));
HY = [Hya Yy];
Hzb = [];
for t = 0:n_z/2
    if (t >= 0 & t <= 3*n_z/8)
        hz = 1;
    elseif (t >= 3*n_z/8 & t <= n_z/2)
        hz = (sin(4*pi*t/n_z))^2;
    end
    Hzb = [Hzb hz];
end
Yz = fliplr(Hzb(:,2:n_z/2));

HZ = [Hzb Yz]';
Hmy = meshgrid(HY,1:n_z);
Hmz = meshgrid(HZ,1:n_y)';
absorption_layer=(Hmy.*Hmz);
clear Hmy;
clear Hmz;

g_tilda = exp(-kz.^2*sigmaz^2/2);
p_ez_e0_tilda = g_tilda.*cos(kz*ant_height);
p_ez_e0_tilda = meshgrid(p_ez_e0_tilda,1:n_y).';
p_ez_0_tilda = p_ez_e0_tilda;
p_ez_tilda = p_ez_0_tilda.*absorption_layer;
propagator = exp(1i*kx*dx);

x_grid_m = 0:dx:max_range_m;
x_output_grid = x_grid_m(1:x_output_filter:end);
y_grid_m = (0:dy_m:(n_y-1)*dy_m) - (n_y)*dy_m/2;
y_output_grid = y_grid_m(1:y_output_filter:end);
z_grid_m = 0:dz_m:(n_z-1)*dz_m/2;
z_output_grid = z_grid_m(1:z_output_filter:end);
res.field = zeros(length(x_output_grid), length(y_output_grid), length(z_output_grid));

wbar = waitbar(0, 'Please wait...');
for iter=1:n_iters
    for ind=2:length(x_grid_m)
        x = x_grid_m(ind);
        waitbar(x / x_grid_m(end), wbar, sprintf('Forward propagation %i of %i', iter, n_iters));
        p_ez_tilda_Dx = p_ez_tilda.*propagator;
        p_ez = Dkz*Dky*(n_y*n_z*ifft2(p_ez_tilda_Dx))/(2*pi)^2;
        p_ez = p_ez.*absorption_layer;

        for brick_i = 1:length(bricks)
            brick = bricks(brick_i);
            if brick.x_min_m <= x && x < brick.x_max_m
                zi = round(brick.height_m / dz_m) + 1;
                y_mask = brick.y_min_m <= y_grid_m & y_grid_m <= brick.y_max_m;
                y_mask = circshift(y_mask, n_y / 2);
                if brick.x_min_m > x - dx
                    brick.back_faset_u = p_ez(1:zi, y_mask);
                    brick.back_faset_lo = p_ez(end-zi:end, y_mask);
                end
                p_ez(1:zi, y_mask) = 0;
                p_ez(end-zi:end, y_mask) = 0;
            else
                if brick.x_min_m <= x && x - dx < brick.x_max_m
                    zi = round(brick.height_m / dz_m) + 1;
                    y_mask = brick.y_min_m <= y_grid_m & y_grid_m <= brick.y_max_m;
                    y_mask = circshift(y_mask, n_y / 2);
                    if ~isempty(brick.front_faset_u)
                        p_ez(1:zi, y_mask) = brick.front_faset_u;
                        p_ez(end-zi:end, y_mask) = brick.front_faset_lo;
                    end
                end
            end
            bricks(brick_i) = brick;
        end

        if mod(ind-1, x_output_filter) == 0
            iii = (ind-1)/x_output_filter + 1;
            res.field(iii,:,:) = squeeze(res.field(iii,:,:)) + ...
                circshift(p_ez(1:z_output_filter:n_z/2, 1:y_output_filter:end), round(length(y_output_grid)/2), 2).';
        end

        p_ez_z1 = p_ez(1:n_z/2+1,:);
        p_ez_z0 = flipud(p_ez(2:n_z/2,:));

        p_ez_e = [p_ez_z1; p_ez_z0];

        p_ez_tilda = (dz_m*dy_m)*fft2(p_ez_e);
        p_ez_tilda_H = p_ez_tilda.*absorption_layer;

        p_ez_tilda_g = p_ez_tilda_H;
        p_ez_tilda1 = p_ez_tilda_g;
        p_ez_tilda = p_ez_tilda1.*absorption_layer;
    end

    if ~two_way
        break
    end

    p_ez_tilda = p_ez_tilda * 0;

    for ind=length(x_grid_m)-1:-1:1
        x = x_grid_m(ind);
        waitbar(1 - x / x_grid_m(end), wbar, sprintf('Backward propagation %i of %i', iter, n_iters));
        p_ez_tilda_Dx = p_ez_tilda.*propagator;
        p_ez = Dkz*Dky*(n_y*n_z*ifft2(p_ez_tilda_Dx))/(2*pi)^2;
        p_ez = p_ez.*absorption_layer;

        % process bricks
        for brick_i = 1:length(bricks)
            brick = bricks(brick_i);
            if brick.x_min_m <= x && x < brick.x_max_m
                zi = round(brick.height_m / dz_m) + 1;
                y_mask = brick.y_min_m <= y_grid_m & y_grid_m <= brick.y_max_m;
                y_mask = circshift(y_mask, n_y / 2);
                if brick.x_max_m <= x + dx
                    brick.front_faset_u = p_ez(1:zi, y_mask);
                    brick.front_faset_lo = p_ez(end-zi:end, y_mask);
                end
                p_ez(1:zi, y_mask) = 0;
                p_ez(end-zi:end, y_mask) = 0;
            else
                if brick.x_max_m > x && x + dx >= brick.x_min_m
                    zi = round(brick.height_m / dz_m) + 1;
                    y_mask = brick.y_min_m <= y_grid_m & y_grid_m <= brick.y_max_m;
                    y_mask = circshift(y_mask, n_y / 2);
                    if ~isempty(brick.back_faset_u)
                        p_ez(1:zi, y_mask) = brick.back_faset_u;
                        p_ez(end-zi:end, y_mask) = brick.back_faset_lo;
                    end
                end
            end
            bricks(brick_i) = brick;
        end

        if mod(ind-1, x_output_filter) == 0
            iii = (ind-1)/x_output_filter + 1;
            res.field(iii,:,:) = squeeze(res.field(iii,:,:)) + ...
                circshift(p_ez(1:z_output_filter:n_z/2, 1:y_output_filter:end), round(length(y_output_grid)/2), 2).';
        end

        p_ez_z1 = p_ez(1:n_z/2+1,:);
        p_ez_z0 = flipud(p_ez(2:n_z/2,:));

        p_ez_e = [p_ez_z1; p_ez_z0];

        p_ez_tilda = (dz_m*dy_m)*fft2(p_ez_e);
        p_ez_tilda_H = p_ez_tilda.*absorption_layer;

        p_ez_tilda_g = p_ez_tilda_H;
        p_ez_tilda1 = p_ez_tilda_g;
        p_ez_tilda = p_ez_tilda1.*absorption_layer;
    end

    p_ez_tilda = p_ez_tilda * 0;
end

res.x_grid_m = x_output_grid;
res.y_grid_m = y_output_grid;
res.z_grid_m = z_output_grid;

end

