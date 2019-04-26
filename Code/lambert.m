function [delV1, delV2, a, et, theta1, theta2] = lambert(R1,R2,V1,V2,TOF0,mu,delt0)
%% Standard
r1 = norm(R1); r2 = norm(R2);

%% Minmum Calculations
delt = acos(dot(R1,R2)/(norm(R1)*norm(R2)));
if delt0 == 1
    delt = 2*pi - delt;
end

c    = sqrt(r1^2 + r2^2 - 2*r1*r2*cos(delt));
s    = (1/2)*(r1+r2+c);

a_min  = s/2; n_min  = sqrt(mu/a_min^3);
alpham = pi;  betam  = 2*asin(sqrt((s-c)/c));

TOF_min = (1/n_min)*((alpham-betam)-(sin(alpham)-sin(betam)));

%% Iteration
if TOF0 < TOF_min
    amax = a_min;
else
    amax = 5*a_min;
end
error = inf;
amin = 0;
a = amax;
while abs(error) >= 60 && a > 0
    % Alpha Calc
    if TOF0 < TOF_min; alpha =    0 + abs(2*asin(sqrt(( s )/(2*a))));
    else;              alpha = 2*pi - abs(2*asin(sqrt(( s )/(2*a))));
    end
    % Beta Calc
    if delt < pi;      beta  =  abs(2*asin(sqrt((s-c)/(2*a))));
    else;              beta  = -abs(2*asin(sqrt((s-c)/(2*a))));
    end
    % n Calc
    n = sqrt(mu/a^3);
    % TOF Guess
    TOF = (1/n)*((alpha-beta)-(sin(alpha)-sin(beta)));
    % Error Comparison
    error = (TOF - TOF0);
    
    if error > 0
        amax = a;
        a = (a+amin)/2;
    else
        amin = a;
        a = (a+amax)/2;
    end
end

%% e, p, and thetas
p1 = 4*a*(s-r1)*(s-r2)*(sin((alpha+beta)/2)/c)^2;
p2 = 4*a*(s-r1)*(s-r2)*(sin((alpha-beta)/2)/c)^2;
e1 = sqrt(1-p1/a); e2 = sqrt(1-p2/a);

if (TOF0 < TOF_min && delt < pi) || (TOF0 > TOF_min && delt > pi)
    e = min([e1,e2]);           p = max([p1,p2]);
    theta1 = acos((p/r1-1)/e);  theta2 = acos((p/r2-1)/e);
else
    e = max([e1,e2]);           p = min([p1,p2]);
    theta1 = acos((p/r1-1)/e);  theta2 = -acos((p/r2-1)/e);
end

%% V1 and V2 Calculation
h = sqrt(mu*a*(1-e^2));
V1rth = [mu/h*e*sin(theta1);mu/h*(1+e*cos(theta1));0];
V2rth = [mu/h*e*sin(theta2);mu/h*(1+e*cos(theta2));0];

%% eccentricity
ht = cross(R1,V1rth);                % Angular Momentum Vector (km^2/s)
et = cross(V1rth,ht)/mu-R1/norm(R1); % Eccentricity Vector  (N/A)


%% Delta V Calculation
h_hat = cross(R1,R2)/norm(cross(R1,R2));
if delt > pi; h_hat = -h_hat; end;

r1hat = R1/(norm(R1));          r2hat = R2/(norm(R2));
t1hat = cross(h_hat,r1hat);     t2hat = cross(h_hat,r2hat);
Q1    = [r1hat';t1hat';h_hat']; Q2    = [r2hat';t2hat';h_hat'];
V1t   = Q1'*V1rth;              V2t   = Q2'*V2rth;
delV1 = norm(V1t-V1);           delV2 = norm(V2t-V2);
