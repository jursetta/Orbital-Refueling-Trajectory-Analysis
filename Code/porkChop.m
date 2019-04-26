function porkChop(DelV, minDel, Departure, TOF)

Departurevec = linspace(min(Departure),max(Departure), 50);

MinArrival = min(Departure) + min(TOF);
MaxArrival = max(Departure) + max(TOF);
Arrivalvec = linspace(MinArrival, MaxArrival, 50);

DeltaVmat = inf(50,50);

[row,col] = size(DelV);
for i = 1:row
    for j= 1:col
        DepartTime = Departure(i);
        DepartInd  = find(DepartTime >= Departurevec,1,'last');
        AriveTime  = DepartTime + TOF(j);
        AriveInd   = find(AriveTime >= Arrivalvec, 1, 'last');
        
        if DelV(i,j) < DeltaVmat(DepartInd,AriveInd)
            DeltaVmat(DepartInd,AriveInd) = DelV(i,j);
        end
        
    end
end

figure
surf(Arrivalvec,Departurevec,DeltaVmat)

disp('Done')
