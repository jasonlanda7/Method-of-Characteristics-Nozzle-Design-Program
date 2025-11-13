
% =========================================================================
%  METHOD OF CHARACTERISTICS (MOC) FOR AXISYMMETRIC NOZZLE DESIGN
% =========================================================================
%
%  Author: [Jason Da Silva]
%  Credits: Characteristic-line formulation and base logic adapted from
%           VDEngineering on youtube! (2025), 
%
%  Date: [11/12/2025]
%  Version: 1.0
%
% -------------------------------------------------------------------------
%  PURPOSE:
%  --------
%  This MATLAB program designs a *shock-free, minimum-length, axisymmetric
%  supersonic nozzle* using the Method of Characteristics (MOC).
%
%  The code computes the internal expansion region by tracing C⁺ and C⁻
%  characteristic lines and determining their intersections to solve for
%  the local flow angle (θ), Prandtl–Meyer function (ν), and Mach number (M).
%  The resulting wall contour provides a perfectly expanded, shock-free
%  supersonic flow field from the throat (M=1) to the desired exit Mach number.
%
%  The flowfield is calculated under the axisymmetric compatibility relations:
%       d(θ + ν) = -tan(μ)/r * dr    along C⁺
%       d(θ - ν) = +tan(μ)/r * dr    along C⁻
%
%  where:
%       θ  = local flow turning angle
%       ν  = Prandtl–Meyer function
%       μ  = Mach angle = sin⁻¹(1/M)
%       r  = local radius
%
%  The program constructs the entire characteristic net (expansion fan and
%  reflection region) and the resulting nozzle contour. It can also attach
%  a parabolic bell section for practical rocket nozzle shaping.
%
% -------------------------------------------------------------------------
%  FEATURES:
%  ----------
%   • Computes both planar and axisymmetric MOC fields
%   • Generates all C⁺ and C⁻ characteristic lines
%   • Produces a smooth, shock-free wall contour
%   • Allows hybrid parabolic bell extension
%   • Fully visualized characteristic net (blue lines)
%
% -------------------------------------------------------------------------
%  OUTPUT:
%  -------
%   - Wall coordinates (x, y) normalized by throat radius Rₜ
%   - Characteristic lines (C⁺ and C⁻)
%   - Centerline and throat geometry
%   - Optional parabolic bell extension
%
% -------------------------------------------------------------------------
%  HOW TO RUN:
%  -----------
%   Example:
%       out = moc_axisymmetric_nozzle(3.0, 1.4, 25, true);
%
%   Inputs:
%       Me     = desired exit Mach number
%       gamma  = ratio of specific heats
%       N      = number of characteristic rays
%       doPlot = true/false (display results)
%
% -------------------------------------------------------------------------
%  NOTES:
%  ------
%   - The MOC produces an *ideal*, inviscid, isentropic design.
%   - True shock-free performance occurs only at the design condition.
%   - For off-design conditions, flow separation or shocks may appear.
%
% =========================================================================
%   - For the initial Phase of the nozzle simulation we are going to use the 
%   - isentropic relationship equations in order to find the exit Mach number
%   - which is specified as (Me).
%   - The variables we are given to calculate the exit Mach is as follows:
%   
%   Qdot  ->  Total amount of energy in the system 
%   Gamma ->  Specific Heat Ratio at Chamber
%   T0    ->  Stagnation Temperature
%   P0    ->  Stagnation Pressure (Chamber Pressure)
%   P0/Pe ->  Pressure Ratio
%
%   - Ideally we would want the ambient pressure Patm to match the exit pressure
%   - Pe at the nozzle exit.
%   - All units are in International Metric Units.
% =========================================================================
%% ISENTROPIC CALCULATIONS FOR EXIT MACH NUMBER.
function out = moc_bell_nozzle(P0, g, TR)
   




Pe = 101325; %Pa
Me = sqrt(((2/(g-1))*((P0/Pe)^((g-1)/g)-1))); 
 

%==========================================================================
%Throat Radius (mm)
DtoR = pi/180;
RtoD = 180/pi;
P = []; %x axis points
%% Prandyl Meyer Function P(M)
A = sqrt((g+1)/(g-1));
B = (g-1)/(g+1);
V_PM = @(x) A*atan(sqrt(B*(x^2-1))) - atan(sqrt(x^2-1));
%% Calculate T_MAX, Then Break into Divisions
T_max = 0.5*(V_PM(Me))*RtoD;
DT = (90-T_max) - fix(90-T_max);
T(1) = DT*DtoR;
n = T_max*2;

for m = 2:n+1
    T(m) = (DT+(m-1))*DtoR;
    %Mach Number from T(i) using T(i) - V_PM (False Position)
    x_int = [1 1.01*Me];
    func = @(x) T(m) - V_PM(x);
    M(m) = fzero(func,x_int);
    P(m) = 0 + TR*tan(T(m));%X-Axis Points
    %Right Running Slopes
    RR(m) = -TR/P(m);
    %Left Running Slopes
    LR(m) = tan(T(m)+asin(1/M(m)));
    SL(m) = -RR(m);
end
%%Plotting 
P(1) = [];
l = length(P);

for j = 1:l
    P1 = [0 TR];
    P2 = [P(j) 0];
    plot(P2,P1,'c')
    hold on
    xlabel('Centerline')
    ylabel('Radius')
end
hold on;
LR(1) = []; RR(1) = [];
SL(1) = [];
F = RR(m-1);

for c = 1:length(P)-1
    x(c) = (TR+SL(c)*P(c))/(SL(c)-F);
    y(c) = F*x(c)+TR;
    X_P = [P(c) x(c)];
    Y_P = [0 y(c)];
    plot(X_P,Y_P,'c');
end
hold on 

%% === Wall Section Generation ===
TM = T_max*DtoR;           % Wall angle (max turning angle)
xw(1) = (TR + SL(1)*P(1)) / (SL(1) - tan(TM));
yw(1) = tan(TM)*xw(1) + TR;

% Plot the first wall segment
plot([P(1) xw(1)], [0 yw(1)], 'c', 'LineWidth', 1.5)

% Initialize slope array
DTW = tan(TM) / (length(P) - 1);
s(1) = tan(TM);
b(1) = TR;

for k = 2:length(P)-1
    % Slope of wall line for this section
    s(k) = tan(TM) - (k - 1) * DTW;
    
    % Intercept of the wall line
    b(k) = yw(k - 1) - s(k) * xw(k - 1);
    
    % Intersection of wall line with C- characteristic from grid point
    xw(k) = (b(k) + SL(k) * P(k)) / (SL(k) - s(k));
    yw(k) = s(k) * xw(k) + b(k);
    
    % Plot the wall segment between previous and current wall points
    plot([xw(k - 1) xw(k)], [yw(k - 1) yw(k)], 'c', 'LineWidth', 2)
end

for k = 2:length(P)-1
    s(k) = tan(TM)-(k-1)*DTW; %Slope
    b(k) = yw(k-1) - s(k)*xw(k-1); %y-int
    xw(k) = (b(k)+SL(k)*P(k))/(SL(k)-s(k));
    yw(k) = s(k)*xw(k)+b(k);
    X_P3 = [x(k) xw(k)];
    Y_P3 = [y(k) yw(k)];
    plot(X_P3,Y_P3,'m');
end
hold on 

%Last Point
%xf = (b(length(b))+SL(length(SL))*P(length(P)))/SL(length(b));
%yf = b(length(b));
%X_F = [P(length(P)) xf];
%Y_F = [0 yf];
%plot(X_F,Y_F,'r');

%xw = [0 xw];
%yw = [TR yw];

% If you also want to save dataMatrix, use another sheet (to avoid overwriting)
writematrix(transpose(xw), 'PARAMS.xlsx', 'Sheet', 'PTS', 'Range', 'A1');
writematrix(transpose(yw), 'PARAMS.xlsx', 'Sheet', 'PTS', 'Range', 'B1');

%% === Determine nozzle exit from final characteristic ===
k_end = length(P) - 1;   % last valid characteristic intersection

x_exit = xw(k_end);
y_exit = yw(k_end);

plot([0 0], [TR -TR], 'w--', 'LineWidth', 2)

% Append exit to wall arrays if not already included
xw = xw(1:k_end);
yw = yw(1:k_end);
xw(end+1) = x_exit;
yw(end+1) = y_exit;
%% === Mirror the Nozzle for Axisymmetric Visualization ===

% --- Plot upper wall up to the final characteristic ---
plot(xw, yw, 'c', 'LineWidth', 2)
hold on

% --- Mirror the wall ---
xw_mirror = xw;
yw_mirror = -yw;
plot(xw_mirror, yw_mirror, 'c', 'LineWidth', 2)

% --- Draw vertical exit line at that last characteristic ---
plot([x_exit x_exit], [y_exit -y_exit], 'w--', 'LineWidth', 2)

% --- Centerline
yline(0, 'w--', 'LineWidth', 1.5)

%% === Plotting

% --- Plot settings
axis equal
grid on
xlabel('Axial Distance (mm)')
ylabel('Radius (mm)')
title('Full Axisymmetric Nozzle Contour (Method of Characteristics)')

% --- Optional: === Save full contour coordinates to Excel ===
full_x = [xw, fliplr(xw_mirror)];
full_y = [yw, fliplr(yw_mirror)];

writematrix([full_x(:), full_y(:)], 'PARAMS.xlsx', 'Sheet', 'FullContour', 'Range', 'A1');





