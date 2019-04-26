function [R,V,M,mu] = planetLoc(planet,t)
% Mu Constants
M_S  = 1.989e30;  mu_S = 1.327124400E11;
M_e  = 5.9742E24; mu_e = 3.986004415E5;
M_j  = 1.898E27;  mu_j = 1.268E8;
M_n  = 1.0278E26; mu_n = 6.809E6;

%  Earth's location as of January 1, 2018
a_e = 149598023;            % km
ta_e = 356.935*pi/180;      % rad
Pe = 2*pi*sqrt(a_e^3/mu_S); % sec
% Jupiter's location as of January 1, 2018
a_j = 778298361;            % km
ta_j = 206.557*pi/180;      % rad
Pj = 2*pi*sqrt(a_j^3/mu_S); % sec
Ld = findL1(M_j,M_S,a_j);
% Neptune's location as of January 1, 2018
a_n = 4504449769;           % km
ta_n = 302.531*pi/180;      % rad
Pn = 2*pi*sqrt(a_n^3/mu_S); % sec

% Assign constants based on planet type
if planet == 'e'
    a = a_e; ta_0 = ta_e; P = Pe; M = M_e; mu = mu_e;
elseif planet == 'j'
    a = Ld; ta_0 = ta_j; P = Pj; M = M_j; mu = mu_j;
else
    a = a_n; ta_0 = ta_n; P = Pn; M = M_n; mu = mu_n;
end
% Calculate rate
rate = 2*pi/P;
ta = ta_0 + t*rate;
V0 = 2*pi*a/P;
% calculate Location and Velocity
R = [ cos(ta)*a;  sin(ta)*a;  0];
V = [-sin(ta)*V0; cos(ta)*V0; 0];