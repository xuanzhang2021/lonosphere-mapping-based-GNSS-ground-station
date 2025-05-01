function [pos,el, az, dop, is_faultout, omc,C,A,PL,WSSE_sqrt,T_threshold] = leastSquarePos(satpos,obs,settings)
%Function calculates the Least Square Solution.
%
%[pos, el, az, dop] = leastSquarePos(satpos, obs, settings);
%
%   Inputs:
%       satpos      - Satellites positions (in ECEF system: [X; Y; Z;] -
%                   one column per satellite)
%       obs         - Observations - the pseudorange measurements to each
%                   satellite corrected by SV clock error
%                   (e.g. [20000000 21000000 .... .... .... .... ....]) 
%       settings    - receiver settings
%
%   Outputs:
%       pos         - receiver position and receiver clock error 
%                   (in ECEF system: [X, Y, Z, dt]) 
%       el          - Satellites elevation angles (degrees)
%       az          - Satellites azimuth angles (degrees)
%       dop         - Dilutions Of Precision ([GDOP PDOP HDOP VDOP TDOP])
%       v           -velocity of receiver in ECEF
%--------------------------------------------------------------------------
%                           SoftGNSS v3.0
%--------------------------------------------------------------------------
%Based on Kai Borre
%Copyright (c) by Kai Borre
%Updated by Darius Plausinaitis, Peter Rinder and Nicolaj Bertelsen
%
% CVS record:
% $Id: leastSquarePos.m,v 1.1.2.12 2006/08/22 13:45:59 dpl Exp $
%==========================================================================
is_faultout = false;  % 初始化为无故障
is_fault = false;  % 初始化为无故障

%=== Initialization =======================================================
nmbOfIterations = 10;

dtr     = pi/180;
pos     = zeros(4, 1);   % center of earth
X       = satpos;
nmbOfSatellites = size(satpos, 2);

A       = zeros(nmbOfSatellites, 4);
omc     = zeros(nmbOfSatellites, 1);
weight = ones(nmbOfSatellites, 1);
az      = zeros(1, nmbOfSatellites);
el      = az;


%=== Iteratively find receiver position ===================================
for iter = 1:nmbOfIterations

    for i = 1:nmbOfSatellites
        if iter == 1
            %--- Initialize variables at the first iteration --------------
            Rot_X = X(:, i);
            trop = 2; 
        else
            %--- Update equations -----------------------------------------
            rho2 = (X(1, i) - pos(1))^2 + (X(2, i) - pos(2))^2 + ...
                   (X(3, i) - pos(3))^2;
            traveltime = sqrt(rho2) / settings.c ;

            %--- Correct satellite position (do to earth rotation) --------
            % Convert SV position at signal transmitting time to position 
            % at signal receiving time. ECEF always changes with time as 
            % earth rotates.
            Rot_X = e_r_corr(traveltime, X(:, i));
            
            %--- Find the elevation angel of the satellite ----------------
            [az(i), el(i), ~] = topocent(pos(1:3, :), Rot_X - pos(1:3, :));

            if (settings.useTropCorr == 1)
                %--- Calculate tropospheric correction --------------------
                trop = tropo(sin(el(i) * dtr), ...
                             0.0, 1013.0, 293.0, 50.0, 0.0, 0.0, 0.0);
            else
                % Do not calculate or apply the tropospheric corrections
                trop = 0;
            end

            if(i==2)
                sm = 1;
            elseif(i==5)
                sm = 1;
            else 
                sm = 0.4;
            end

            %weight(i)=sin(el(i)*3.1415926/180)*sm;
            weight(i) = 1;
        end % if iter == 1 ... ... else 
        

        %--- Apply the corrections ----------------------------------------
        omc(i) = ( obs(i) - norm(Rot_X - pos(1:3), 'fro') - pos(4) - trop ); 

        %--- Construct the A matrix ---------------------------------------
        A(i, :) =  [ (-(Rot_X(1) - pos(1))) / norm(Rot_X - pos(1:3), 'fro') ...
                     (-(Rot_X(2) - pos(2))) / norm(Rot_X - pos(1:3), 'fro') ...
                     (-(Rot_X(3) - pos(3))) / norm(Rot_X - pos(1:3), 'fro') ...
                     1 ];
        
        % weight by hd
        
    end % for i = 1:nmbOfSatellites
    % These lines allow the code to exit gracefully in case of any errors
    if rank(A) ~= 4
        pos     = zeros(1, 4);
        dop     = inf(1, 5);
        fprintf('Cannot get a converged solotion! \n');
        return
    end

    %--- Find position update (in the least squares sense)-----------------
    %x   = A \ omc;

    %--- Find position update (for the weighted least square)
    weight = weight./sum(weight);

    W = diag(weight);
    C=W'*W;
    x=(A'*C*A)\(A'*C*omc);
    
    %omc is the observed value


    pos = pos + x;
    % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % %--- Apply position update --------------------------------------------
    % if norm(x)<1e-4
    %         %%%%%%%%%%%%%%%%%%%%%%% RAIM insert %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        [is_fault, excluded_idx] = raim_detection(A, omc, diag(C), settings);
        if is_fault == 1
            is_faultout = 1;
        end
        if is_fault && (nmbOfSatellites - 1 >= 4)  % Ensure that the remaining number of satellites is ≥ 4
            % Reinitialize the data after eliminating the faulty satellite
            X(:, excluded_idx) = [];
            obs(excluded_idx) = [];
            %doppler(excluded_idx) = [];
            nmbOfSatellites = nmbOfSatellites - 1;
            A(excluded_idx, :) = [];
            omc(excluded_idx) = [];
            weight(excluded_idx) = [];
            C = diag(weight.^2);  % Update the weight matrix
            % Reset the iteration to recalculate
            iter = 0;  % The next round of the loop will start from iter=1
            pos     = zeros(4, 1);   % center of earth
            continue;
        end
    end

end % for iter = 1:nmbOfIterations




%--- Fixing resulut -------------------------------------------------------
pos = pos';

%=== Calculate Dilution Of Precision ======================================
if nargout > 4
    %--- Initialize output ------------------------------------------------
    dop     = zeros(1, 5);
    
    %--- Calculate DOP ----------------------------------------------------
    Q       = inv(A'*A);
    
    dop(1)  = sqrt(trace(Q));                       % GDOP    
    dop(2)  = sqrt(Q(1,1) + Q(2,2) + Q(3,3));       % PDOP
    dop(3)  = sqrt(Q(1,1) + Q(2,2));                % HDOP
    dop(4)  = sqrt(Q(3,3));                         % VDOP
    dop(5)  = sqrt(Q(4,4));                         % TDOP
end  % if nargout  == 4

%===To calculate receiver velocity=====HD
% b=[];
% satvelocity=satvelocity';
% for i=1:nmbOfSatellites
%     b(i)=-doppler(i)-satvelocity(i,:)*(A(i,1:3))';
% end
% v=(A'*C*A)\(A'*C*b');
% 
% 
% if size(el,2) == 4 
%     el = [el,0];
% end
% if size(az,2) == 4 
%     az = [az,0];
% end
% if size(omc,2) == 4 
%     omc = [omc,0];
% end
% if size(C,2) == 4 
%     C(5,5) = 0;
% end
% if size(A,2) == 4 
%     A(5,5) = 0;
% end
    
PL = 1;
WSSE_sqrt = 1;
T_threshold = 1;
    Nr_sat = 5;
    I = eye(Nr_sat);   
    isolation_mat = ones(Nr_sat, 1);
    I = I * diag(isolation_mat);

    y = omc;
    S = inv(A'*C*A)*(A'*C);
    P = A*S;
    W = C;

    %WSSE_sqrt = sqrt(y'*W*(I-P)*y ./ (sum(isolation_mat)-4));
    WSSE_sqrt = sqrt(y'*W*(I-P)*y);

    SLPOE = sqrt((S(1,:).^2+S(2,:).^2+S(3,:).^2)./(1-diag(P)'));

    SSE = omc' * C *omc;

    P_fa = 1e-2;      
    n_satellites = 5;  
    dof = n_satellites - 4; 

    
    T_threshold = sqrt(chi2inv(1 - P_fa, dof));

    P_md = 1e-7;               
    K_md = norminv(1 - P_md); 
    sigma = 3;

    K = S;

    for i = 1 : Nr_sat
    % OLS
     Pslope(i) = sqrt(sum((K(1:3,i)).^2)) * sqrt(Nr_sat-4) / sqrt(1-P(i,i)); 
    % WLS
     %Pslope(i) = sqrt(sum((K(1:3,i)).^2)) * sqrt(1/W(i,i)) / sqrt(1-P(i,i));
    end
    SLP = sort(Pslope,'descend');
    PL = SLP(2) * T_threshold +  K_md * sigma;
