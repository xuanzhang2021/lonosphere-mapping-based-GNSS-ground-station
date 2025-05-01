
function [distance_lat,distance_lon] = latlon_to_meters(latGT, lonGT, lat, lon)
    % Earth's Radius (meters)
    R = 6371000;  
    
    % Latitude difference and longitude difference
    latitude_diff = lat - latGT;
    longitude_diff = lon - lonGT;
    
    distance_lat = latitude_diff * (pi / 180) * R;  
    distance_lon = longitude_diff * (pi / 180) * R * cosd(latGT);  
    
end

gt = [22.328444770087565; 114.1713630049711];

[error_lat,error_lon] = latlon_to_meters(22.328444770087565, 114.1713630049711, navSolutions.latitude, navSolutions.longitude);

error_height = navSolutions.height - 3;

function PL = calculate_3D_PL(n_sat, sigma, P_fa, P_md)
    dof = n_sat - 4; % Degree of freedom
    T_threshold = chi2inv(1 - P_fa, dof) * sigma^2;
    
    fun = @(PL) ncx2cdf(T_threshold, dof, (PL^2)/sigma^2) - (1 - P_md);
    PL_guess = 50; % Initial guess value
    options = optimset('Display','off');
    PL = fzero(fun, PL_guess, options);
end

sigma = 3;          % Standard deviation of pseudo-range noise [m]
P_fa = 1e-2;        % False alarm probability
P_md = 1e-7;        % Probability of missed detection
AL = 50;            % 3D Alarm Limit Value [m]
n_sat = 5;          % The number of satellites
num_samples = 79;

PL_3D = 29.63;
fprintf('3D Protection Level: %.2f meters\n', PL_3D);

% Generate simulation data (VPE and VPL)
rng(0); % Fix the random seeds for reproduction
VPE = error_lat'; % Vertical positioning error (normal distribution)
VPL = PL_3D*ones(num_samples,1); % Vertical protection level (with random fluctuations)

x_edges = linspace(0, 60, 50);
y_edges = linspace(0, 60, 50);
[N, ~, ~] = histcounts2(VPE, VPL, x_edges, y_edges);

% Stanford Chart
figure;
h = pcolor(x_edges(1:end-1), y_edges(1:end-1), log10(N' + 1)); 
set(h, 'EdgeColor', 'none');
colormap(turbo); 
colorbar;
hold on;

plot(xlim, [AL AL], 'r--', 'LineWidth', 2, 'DisplayName', sprintf('AL=%.1fm', AL));

plot([0 max(xlim)], [0 max(ylim)], 'r --', 'LineWidth', 1.5, 'DisplayName', 'VPE = VPL');

plot([AL AL], ylim, 'r--', 'LineWidth', 2, 'DisplayName', sprintf('AL=%.1fm', PL_3D));

% fill([min(xlim) max(xlim) max(xlim) min(xlim)], ...
%      [0 0 AL AL], 'g', 'FaceAlpha', 0.5, 'EdgeColor', 'none', ...
%      'DisplayName', 'Safe Zone (VPL < AL)');

% fill([min(xlim) max(xlim) max(xlim) min(xlim)], ...
%      [AL AL max(ylim) max(ylim)], 'r', 'FaceAlpha', 0.8, 'EdgeColor', 'none', ...
%      'DisplayName', 'Hazard Zone (VPL > AL)');
% 
% fill([min(xlim) max(xlim) max(xlim) min(xlim)], ...
%      [AL AL max(xlim) max(xlim)], 'r', 'FaceAlpha', 0.8, 'EdgeColor', 'none', ...
%      'DisplayName', 'Hazard Zone (VPL > AL)');


xlabel('Vertical Position Error (VPE) [m]');
ylabel('Vertical Protection Level (VPL) [m]');
title(sprintf('Stanford Chart Analysis\nn=%d satellites, σ=%.1fm, AL=%.1fm', n_sat, sigma, AL));
legend('Location', 'northeast');
grid on;

text(0.05*max(xlim), 0.9*max(ylim), ...
     sprintf('P_{fa}=%.0e\nP_{md}=%.0e', P_fa, P_md), ...
     'BackgroundColor', 'white', 'FontSize', 10);

xlim([0 max(x_edges)]);
ylim([0 max(y_edges)]);

function [PL_3D, HPL, VPL] = dynamic_PL(n_sat, sigma, P_fa, P_md, A, W)
    % Input real-time parameters: Number of satellites n_sat, noise sigma, geometric matrix A, weight W
    dof = n_sat - 4;
    
    % 3D PL calculation
    T_thresh = chi2inv(1-P_fa, dof) * sigma^2;
    fun = @(PL) ncx2cdf(T_thresh, dof, (PL^2)/sigma^2) - (1 - P_md);
    PL_3D = fzero(fun, 50);
    
    % Calculate by direction PL
    S = inv(A'*W*A); 
    K_H = 5.33;      
    K_V = 5.33;      
    
    HPL = K_H * sqrt(max(S(1,1), S(2,2)));
    VPL = K_V * sqrt(S(3,3));
end


load('skymask_A1_urban.csv');
azimuth = skymask_A1_urban(:,1);
elevation = skymask_A1_urban(:,2);

% ---Draw the polar coordinate celestial potential map ---
figure;
polarplot(deg2rad(azimuth), 90 - elevation, 'b-', 'LineWidth', 1.5); 
hold on;

az0 = [57.3027107992585
139.826705356260
21.8318470521365
41.5483282207439]; % Satellite azimuth Angle
el0 = [73.8679671141010
22.2923826402396
59.9140197171908
46.9124815902099]; % Satellite elevation Angle

theta = deg2rad(az0);     
r = 90 - el0;          

ssl = polarplot(theta(1), r(1), 'go', 'MarkerSize', 10, 'LineWidth', 2);
ssn = polarplot(theta(2), r(2), 'ro', 'MarkerSize', 10, 'LineWidth', 2);
polarplot(theta(3), r(3), 'go', 'MarkerSize', 10, 'LineWidth', 2);
polarplot(theta(4), r(4), 'ro', 'MarkerSize', 10, 'LineWidth', 2);

title('Zenith map (polar coordinate)');
rlim([0 90]);           
set(gca, 'ThetaDir', 'clockwise');     
set(gca, 'ThetaZeroLocation', 'top');  
ax = gca;
ax.RAxis.Label.String = 'Zenith Angle (90 - EL)';
legend([ssn, ssl], {'LOS', 'NLOS'}, 'Location', 'best');set(gca, 'ThetaDir', 'clockwise'); 