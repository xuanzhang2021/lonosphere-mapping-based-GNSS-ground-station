function plotNavigation(navSolutions, settings)
%Functions plots variations of coordinates over time and a 3D position
%plot. It plots receiver coordinates in UTM system or coordinate offsets if
%the true UTM receiver coordinates are provided.  
%
%plotNavigation(navSolutions, settings)
%
%   Inputs:
%       navSolutions    - Results from navigation solution function. It
%                       contains measured pseudoranges and receiver
%                       coordinates.
%       settings        - Receiver settings. The true receiver coordinates
%                       are contained in this structure.

%--------------------------------------------------------------------------
%                           SoftGNSS v3.0
% 
% Copyright (C) Darius Plausinaitis
% Written by Darius Plausinaitis
%--------------------------------------------------------------------------
%This program is free software; you can redistribute it and/or
%modify it under the terms of the GNU General Public License
%as published by the Free Software Foundation; either version 2
%of the License, or (at your option) any later version.
%
%This program is distributed in the hope that it will be useful,
%but WITHOUT ANY WARRANTY; without even the implied warranty of
%MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%GNU General Public License for more details.
%
%You should have received a copy of the GNU General Public License
%along with this program; if not, write to the Free Software
%Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301,
%USA.
%--------------------------------------------------------------------------

% CVS record:
% $Id: plotNavigation.m,v 1.1.2.25 2006/08/09 17:20:11 dpl Exp $

%% Plot results in the necessary data exists ==============================
if (~isempty(navSolutions))

    %% If reference position is not provided, then set reference position
    %% to the average postion
    if isnan(settings.truePosition.E) || isnan(settings.truePosition.N) ...
                                      || isnan(settings.truePosition.U)

        %=== Compute mean values ========================================== 
        % Remove NaN-s or the output of the function MEAN will be NaN.
        refCoord.E = mean(navSolutions.E(~isnan(navSolutions.E)));
        refCoord.N = mean(navSolutions.N(~isnan(navSolutions.N)));
        refCoord.U = mean(navSolutions.U(~isnan(navSolutions.U)));

        %Also convert geodetic coordinates to deg:min:sec vector format
        meanLongitude = dms2mat(deg2dms(...
            mean(navSolutions.longitude(~isnan(navSolutions.longitude)))), -5);
        meanLatitude  = dms2mat(deg2dms(...
            mean(navSolutions.latitude(~isnan(navSolutions.latitude)))), -5);


        refPointLgText = ['Mean Position\newline  Lat: ', ...
                            num2str(meanLatitude(1)), '{\circ}', ...
                            num2str(meanLatitude(2)), '{\prime}', ...
                            num2str(meanLatitude(3)), '{\prime}{\prime}', ...
                         '\newline Lng: ', ...
                            num2str(meanLongitude(1)), '{\circ}', ...
                            num2str(meanLongitude(2)), '{\prime}', ...
                            num2str(meanLongitude(3)), '{\prime}{\prime}', ...
                         '\newline Hgt: ', ...
                            num2str(mean(navSolutions.height(~isnan(navSolutions.height))), '%+6.1f')];
    else
        refPointLgText = 'Reference Position';
        refCoord.E = settings.truePosition.E;
        refCoord.N = settings.truePosition.N;
        refCoord.U = settings.truePosition.U;        
    end    
     
    figureNumber = 300;
    % The 300 is chosen for more convenient handling of the open
    % figure windows, when many figures are closed and reopened. Figures
    % drawn or opened by the user, will not be "overwritten" by this
    % function if the auto numbering is not used.
 
    %=== Select (or create) and clear the figure ==========================
    figure(figureNumber);
    clf   (figureNumber);
    set   (figureNumber, 'Name', 'Navigation solutions');
 
    %--- Draw axes --------------------------------------------------------
    handles(1, 1) = subplot(4, 2, 1 : 4);
    handles(3, 1) = subplot(4, 2, [5, 7]);
    handles(3, 2) = subplot(4, 2, [6, 8]);    
 


%% Plot all figures =======================================================
    %--- Coordinate differences in UTM system -----------------------------

    plot(handles(1, 1), [(navSolutions.E - refCoord.E)', ...
                         (navSolutions.N - refCoord.N)',...
                         (navSolutions.U - refCoord.U)']);

    title (handles(1, 1), 'Coordinates variations in UTM system');
    legend(handles(1, 1), 'E', 'N', 'U');
    xlabel(handles(1, 1), ['Measurement period: ', ...
                                    num2str(settings.navSolPeriod), 'ms']);
    ylabel(handles(1, 1), 'Variations (m)');
    grid  (handles(1, 1));
    axis  (handles(1, 1), 'tight');   
    
    
    
        %--- Coordinate differences in UTM system -----------------------------
    % plot(handles(1, 1), navSolutions.velocity);
    % 
    % title (handles(1, 1), 'Coordinates variations in UTM system');
    % legend(handles(1, 1), 'E', 'N', 'U');
    % xlabel(handles(1, 1), ['Measurement period: ', ...
    %                                 num2str(settings.navSolPeriod), 'ms']);
    % ylabel(handles(1, 1), 'Variations (m)');
    % grid  (handles(1, 1));
    % axis  (handles(1, 1), 'tight');    
 
    %--- Position plot in UTM system --------------------------------------
    plot3 (handles(3, 1), navSolutions.E - refCoord.E, ...
                          navSolutions.N - refCoord.N, ... 
                          navSolutions.U - refCoord.U, '+');
    hold  (handles(3, 1), 'on');
    %Plot the reference point
    plot3 (handles(3, 1), 0, 0, 0, 'r+', 'LineWidth', 1.5, 'MarkerSize', 10);
    hold  (handles(3, 1), 'off');
    
    view  (handles(3, 1), 0, 90);
    axis  (handles(3, 1), 'equal');
    grid  (handles(3, 1), 'minor');    
    
    legend(handles(3, 1), 'Measurements', refPointLgText);
 
    title (handles(3, 1), 'Positions in UTM system (3D plot)');
    xlabel(handles(3, 1), 'East (m)');
    ylabel(handles(3, 1), 'North (m)');
    zlabel(handles(3, 1), 'Upping (m)');
    
    %--- Satellite sky plot -----------------------------------------------
    skyPlot(handles(3, 2), ...
            navSolutions.az, ...
            navSolutions.el, ...
            navSolutions.PRN(:, 1));
        
    title (handles(3, 2), ['Sky plot (mean PDOP: ', ...
                               num2str(mean(navSolutions.DOP(2,:))), ')']);  
                           
else
    disp('plotNavigation: No navigation data to plot.');
end % if (~isempty(navSolutions))
