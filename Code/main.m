%% Assumptions
% - Earth, Jupiter, and Pluto are in the same orbital plane
% - Earth, Jupiter, and Pluto are in circular orbits
% - Lamberts transfers are done in a Heliocentric system
% - J200 centered inertial frame
close all; clear all; clc;
%% Constants
mu_S = 1.327124400E11;
M_s = 1.989e30;
% Time Constants
hour = 3600;
day = 24*hour;
year = 365*day;
month = year/12;
t0 = datenum(2018,1,2,0,0,0);

% Initial Earth -> Neptune transfer over 1 year time with 5-10 year TOF
tvecn = linspace(0,year,100);
TOFvecn = linspace(1*year,25*year,100);
[DelV_matn, TOF_matn, TOA_matn, time_matn] = deltaVmat(tvecn, TOFvecn, 'e', 'n', mu_S,M_s);
% Minimum Delta-V Calculation
minDelNep = min(min(DelV_matn));
for i = 1:100
    for j = 1:100
        time_matn(i,j) = t0 + ceil(time_matn(i,j)/86400);
    end
end

% Initial Earth->J_L1 transfer over 1 year time with 1-5 year TOF
tvec1 = linspace(0,year,100);
TOFvec1 = linspace(1*year,10*year,100);
[DelV_mat1, TOF_mat1, TOA_mat1, time_mat1] = deltaVmat(tvec1, TOFvec1, 'e', 'j', mu_S,M_s);

for i = 1:100
    for j = 1:100
        time_mat1(i,j) = t0 + ceil(time_mat1(i,j)/86400);
    end
end

% Locate all flights where Delta-V is less than Minimum Earth->Neptune
logDvec1 = DelV_mat1 < minDelNep;
DelV_vec1 = DelV_mat1(logDvec1);
TOF_vec1  =  TOF_mat1(logDvec1);
TOA_vec1  =  TOA_mat1(logDvec1);
time_vec1 = time_mat1(logDvec1);

% J_L1->Neptune transfer starting at all succesful previous transfers with
% 1-5 year TOF
[tvec2,Ii,~] = unique(TOA_vec1);
TOFvec2 = linspace(1*year,10*year,100);
[DelV_mat2, TOF_mat2, TOA_mat2, time_mat2] = deltaVmat(tvec2', TOFvec2, 'j', 'n', mu_S,M_s);

% Create Vector of all Flgith paths which satisfy maximum delta-V min after
% refueling.
logDvec2  = DelV_mat2 < minDelNep;
DelV_vec2 = DelV_mat2(logDvec2);
TOF_vec2  =  TOF_mat2(logDvec2);
TOA_vec2  =  TOA_mat2(logDvec2);
time_vec2 = time_mat2(logDvec2);

% Plots
figure
contourf(time_matn,TOF_matn/year,DelV_matn)
title('Direct Flight Path from Earth to Neptune')
c = colorbar;
xlabel('Departure Date in 2018')
ylabel('Time of Flight (years)')
c.Label.String = 'Delta-V (km/s)';
NumTicks = 13;
L = get(gca,'XLim');
set(gca,'XTick',linspace(L(1),L(2),NumTicks))
datetick('x','mmm','keeplimits', 'keepticks')
set(gca,'FontSize',20)
set(gcf, 'Position', [100, 100, 1300, 9000])
saveas(gcf,'DirectFull.png')

[time_vec1U,Ind] = unique(time_vec1(Ii));
tvec2_mat = repmat(time_vec1U,1,100);

TOA_mat2U = TOA_mat2(Ind,:);
DelV_mat2U = DelV_mat2(Ind,:);
LogFinal = DelV_mat2U>minDelNep;

figure
contourf(tvec2_mat,TOA_mat2U/year,DelV_mat2U,18.5015:1:45.0015)
title('Refueled Flight Path from Earth to Neptune')
c = colorbar;
xlabel('Departure Date in 2018')
ylabel('Time of Arrival after January 1st, 2018 (years)')
c.Label.String = 'Delta-V (km/s)';
NumTicks = 11;
L = get(gca,'XLim');
set(gca,'XTick',linspace(L(1),L(2),NumTicks))
datetick('x','mmm','keeplimits', 'keepticks')
set(gca,'FontSize',20)
set(gcf, 'Position', [100, 100, 1300, 9000])
saveas(gcf,'RefueledFull.png')

figure
contour(tvec2_mat,TOA_mat2U/year,DelV_mat2U,[18.5015:1:45.0015],'LineWidth',2)
hold on
DelV_mat2U(LogFinal) = -inf;
contourf(tvec2_mat,TOA_mat2U/year,DelV_mat2U,0.0015:21.0015)
title('Refueled Flight Path from Earth to Neptune')
c = colorbar;
xlabel('Departure Date in 2018')
ylabel('Time of Arrival after January 1st, 2018 (years)')
c.Label.String = 'Delta-V (km/s)';
NumTicks = 11;
L = get(gca,'XLim');
set(gca,'XTick',linspace(L(1),L(2),NumTicks))
datetick('x','mmm','keeplimits', 'keepticks')
set(gca,'FontSize',20)
set(gcf, 'Position', [100, 100, 1300, 9000])
hold off
saveas(gcf,'RefueledHighlight.png')


figure
contourf(time_mat1,TOF_mat1/year,DelV_mat1,11.0015:2:60.0015)
title('Refueled Flight Path from Earth to Jupiter')
c = colorbar;
xlabel('Departure Date in 2018')
ylabel('Time of Flight (years)')
c.Label.String = 'Delta-V (km/s)';
NumTicks = 13;
L = get(gca,'XLim');
set(gca,'XTick',linspace(L(1),L(2),NumTicks))
datetick('x','mmm','keeplimits', 'keepticks')
set(gca,'FontSize',20)
set(gcf, 'Position', [100, 100, 1300, 9000])
saveas(gcf,'RefueledFirst.png')


figure
LogFinal2 = DelV_mat1>minDelNep;
contour(time_mat1,TOF_mat1/year,DelV_mat1,11.0015:2:60.0015)
hold on
DelV_mat1(LogFinal2) = -inf;
contourf(time_mat1,TOF_mat1/year,DelV_mat1,11.0015:2:21.0015)
title('Refueled Flight Path from Earth to Jupiter')
c = colorbar;
xlabel('Departure Date in 2018')
ylabel('Time of Flight (years)')
c.Label.String = 'Delta-V (km/s)';
NumTicks = 13;
L = get(gca,'XLim');
set(gca,'XTick',linspace(L(1),L(2),NumTicks))
datetick('x','mmm','keeplimits', 'keepticks')
set(gca,'FontSize',20)
set(gcf, 'Position', [100, 100, 1300, 9000])
saveas(gcf,'RefueledFirstHighlighted.png')