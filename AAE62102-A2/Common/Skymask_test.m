function blocked_sat_idx=Skymask(filename,el,az,ifplot)

% Read CSV file
% filename=settings.skymask_filename;
data = readtable(filename);

% Extract azimuth and elevation
azimuth = data.Azimuth_angle_deg;
elevation = data.Elevation_angle_deg;



el_sat=el;
az_sat=az;


building_el_at_sat = interp1(azimuth, elevation, az_sat(:,end), 'linear', 'extrap');

% Compare satellite elevation with building elevation
is_visible = el_sat(:,end) > building_el_at_sat; % logical array

% Find visible and blocked satellites
% visible_sat_idx = find(is_visible);
blocked_sat_idx = find(~is_visible);

visible_sat_idx = find(is_visible);

if(ifplot)
    % Visualization with correct polar orientation
    % figure;
    polarplot(deg2rad(az_sat(visible_sat_idx)), el_sat(visible_sat_idx), 'go'); hold on;
    polarplot(deg2rad(az_sat(blocked_sat_idx)), el_sat(blocked_sat_idx), 'rx');hold on;
    polarplot(deg2rad(azimuth), elevation, 'b-');
    
    % Adjust polar axes
    ax = gca;
    ax.ThetaDir = 'clockwise';   % Azimuth increases clockwise
    ax.ThetaZeroLocation = 'top'; % 0 degrees at the top
    ax.RDir = 'reverse';          % Elevation increases inward (as in typical sky plots)
    rlim([0 90]);
    title('Satellite Visibility Check');
    legend('Visible Satellites', 'Blocked Satellites', 'Building Profile');
end