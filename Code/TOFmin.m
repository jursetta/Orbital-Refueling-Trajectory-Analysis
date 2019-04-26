function TOFp = TOFmin(R1,R2,mu,delt0)
r1 = norm(R1); r2 = norm(R2);

delt = acos(dot(R1,R2)/(norm(R1)*norm(R2)));
if delt0 == 1
    delt = 2*pi - delt;
end

c    = sqrt(r1^2 + r2^2 - 2*r1*r2*cos(delt));
s    = (1/2)*(r1+r2+c);
TOFp = 1/3*sqrt(2/mu)*(s^(3/2)-(s-c)^(3/2))/3600;