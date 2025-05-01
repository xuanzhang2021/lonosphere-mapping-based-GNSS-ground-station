function [X_k,P_k] = ekf(satPos,satVel,obs,dopplers,settings,X,P,Q)
%UNTITLED 此处显示有关此函数的摘要
%   此处显示详细说明


numSat = size(obs,2);

%% time update
dt = settings.navSolPeriod/1000;  % ms to s
F = eye(8); % State Transition Model
F(1, 4) = dt;  F(2, 5) = dt;  F(3, 6) = dt;
F(7, 8) = dt;

% prediction 
X_kk = F * X;
P_kk = F*P*F'+Q;

%% measurement update
% Measurement Noise Covariance
R_pseudo = 1000;    % Pseudorange noise (meters)
R_doppler = 10000;    % Doppler noise (m/s)

% Initialize measurement model
H = zeros( 2*numSat, size(X,1));
Z = zeros( 2*numSat, 1);
h_x = zeros( 2*numSat, 1);

for i = 1:numSat
    %--- Correct satellite position (do to earth rotation) --------
    Xs = satPos(1, i);  Ys = satPos(2, i);  Zs = satPos(3, i);
    VXs = satVel(1,i); VYs = satVel(2,i); VZs = satVel(3,i);
    dX = Xs - X_kk(1);
    dY = Ys - X_kk(2);
    dZ = Zs - X_kk(3);
    dVX = VXs - X_kk(4);
    dVY = VYs - X_kk(5);
    dVZ = VZs - X_kk(6);
    rho = sqrt(dX^2 + dY^2 + dZ^2); % Range

    traveltime = rho / settings.c ;
    Rot_X = e_r_corr(traveltime, satPos(:, i));

    %--- Find the elevation angel of the satellite ----------------
    % [az(i), el(i), ~] = topocent(X_kk(1:3), Rot_X - X_k(1:3));

    % if (settings.useTropCorr == 1)
    %     %--- Calculate tropospheric correction --------------------
    %     trop = tropo(sin(el(i) * dtr), ...
    %         0.0, 1013.0, 293.0, 50.0, 0.0, 0.0, 0.0);
    % else
    %     % Do not calculate or apply the tropospheric corrections
    %     trop = 0;
    % end
    % weight(i)=sin(el(i))^2;

    Xs = Rot_X(1);  Ys = Rot_X(2);  Zs = Rot_X(3);
    dX = Xs - X_kk(1);
    dY = Ys - X_kk(2);
    dZ = Zs - X_kk(3);
    rho = sqrt(dX^2 + dY^2 + dZ^2); % Range

    % Pseudorange Measurement
    H(i, :) = [-dX/rho, -dY/rho, -dZ/rho, 0,0,0,1,0];
    Z(i) = obs(i) ;
    h_x(i) = rho +  X_kk(7);


    % Doppler Measurement
    if i > length(dopplers)  % Ensure index does not exceed doppler length
        disp('Lack Doppler measurement.');
    else
        H(numSat + i, :) = [0, 0, 0, -dX/rho, -dY/rho, -dZ/rho, 0, 1];
        h_x(numSat + i) = (dX*dVX + dY*dVY + dZ*dVZ ) / rho +  X_kk(8);
        Z(numSat + i) = dopplers(i);
    end


end

% Measurement Noise Covariance Matrix
R = diag([ones(1, numSat) * R_pseudo, ones(1, numSat) * R_doppler]);

% Innovation Calculation
r = Z - h_x;
S = H * P_kk * H' + R;
K = P_kk * H' /S; % Kalman Gain

% Update State Estimate
X_k = X_kk + (K * r);
I = eye(size(X, 1));
P_k = (I - K * H) * P_kk * (I - K * H)' + K * R * K';


end

