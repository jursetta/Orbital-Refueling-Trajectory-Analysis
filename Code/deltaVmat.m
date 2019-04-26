function [DelV_mat, TOF_mat, TOA_mat, t_mat] = deltaVmat(tvec, TOFvec, planet1, planet2, mu,M_s)
% Mat Size Allocations
rows = length(tvec);
cols = length(TOFvec);
% Matrix Allocation
DelV_mat = zeros(rows,cols);
TOF_mat  = zeros(rows,cols);
TOA_mat  = zeros(rows,cols);
t_mat    = zeros(rows,cols);
c1 = 0;
for t = tvec
    c2   = 0; c1 = c1 + 1; % Reset and Update
    for TOF = TOFvec
        TOA = TOF+t;
        c2 = c2 + 1; % Update
        % Planetary Location and Velocity Calcualtion
        [R1,V1,M1,mu1] = planetLoc(planet1, t    );
        [R2,V2,M2,mu2] = planetLoc(planet2, TOA);
        % Lamberts Problem Solve for + and - 180 transfer
        [delV1_1, delV2_1, a1, ~, ~, ~] = lambert(R1,R2,V1,V2,TOF,mu,0);
        [delV1_2, delV2_2, a2, ~, ~, ~] = lambert(R1,R2,V1,V2,TOF,mu,1);
        
        % Departure Patched Conics Calculation (Only if leaving from Earth)
        if planet1 == 'e'
            % Constant Definition
            mu_e   = mu1;  M_e = M1;
            re     = 7000; vc  = sqrt(mu_e/re);
            % Sphere of Influence Calc
            rSOI_e = norm(R1)*(M_e/M_s)^(2/5);
            % Energy and Periapsis Velocity from Hyperbolic Entry
            ESOI1   = delV1_1^2/2 - mu_e/rSOI_e; ESOI2   = delV1_2^2/2 - mu_e/rSOI_e;
            Vp1     = sqrt(2*(ESOI1+mu_e/re));   Vp2     = sqrt(2*(ESOI2+mu_e/re));
            % Change in Velocity from Hyperbolic to Circular orbit
            delV1_1 = abs(Vp1 - vc);
            delV1_2 = abs(Vp2 - vc);
        end
        
        % Arrival Patched Conics Calculation (Only if Arriving at Neptune)
        if planet2 == 'n'
            % Constant Definition
            mu_n    = mu2;   M_n = M2;
            rn      = 28000; vc = sqrt(mu_n/rn);
            % Sphere of Influence Calc
            rSOI_n  = norm(R2)*(M_n/M_s)^(2/5);
            % Energy and Periapsis Velocity from Hyperbolic Entry
            ESOI1   = delV2_1^2/2 - mu_n/rSOI_n; ESOI2   = delV2_2^2/2 - mu_n/rSOI_n;
            Vp1     = sqrt(2*(ESOI1+mu_n/rn));   Vp2     = sqrt(2*(ESOI2+mu_n/rn));
            % Change in Velocity from Hyperbolic to Circular orbit
            delV2_1 = abs(Vp1 - vc);
            delV2_2 = abs(Vp2 - vc);
        end
        % Finals Delta V calculation
        delV01 = delV1_1 + delV2_1;
        delV02 = delV1_2 + delV2_2;
        %(RETURN VALS) Return Values for Calculations
        DelV_mat(c1,c2) = min([delV01,delV02]);
        t_mat(c1,c2)    = t ;
        TOF_mat(c1,c2)  = TOF;
        TOA_mat(c1,c2)  = TOA;
    end
end