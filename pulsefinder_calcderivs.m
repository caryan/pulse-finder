function [derivs,goodness,zderivs] = pulsefinder_calcderivs(prop,Uwant,rhoin,rhogoal,RFmatts,pulse,HTOT,params,corrmat)

nbspins = log2(size(Uwant,1));

plength = size(pulse,1);

%Intitialize the derivatives
derivs = zeros(size(pulse));

%Now calculate forward propagation
Uforward = zeros(2^nbspins,2^nbspins,plength);
Uforward(:,:,1) = prop(:,:,1);
for ct = 2:1:plength
    Uforward(:,:,ct) = prop(:,:,ct)*Uforward(:,:,ct-1);
end

%And now the backwards
Uback = zeros(2^nbspins,2^nbspins,plength);

%Whether we are doing unitary or state this is a bit different
if(params.searchtype == 1)
    Uback(:,:,1) = -Uwant;
else
    Uback(:,:,1) = eye(2^nbspins);
end

for ct = 2:1:plength
    Uback(:,:,ct) = prop(:,:,plength-ct+2)'*Uback(:,:,ct-1);
end
Uback = flipdim(Uback,3);

%Finally the derivatives (again this is different for state to state or
%unitary)
if(params.searchtype == 1) %Do the derivatives for unitary
    tmpRFmat = zeros(size(RFmatts,3),2^(2*nbspins));
    
    for RFmattct = 1:1:size(RFmatts,3)
        tmpRFmat(RFmattct,:) = reshape(transpose(RFmatts(:,:,RFmattct)),[1 2^(2*nbspins)]);
    end
    
    for ct = 1:1:plength
        tmpmat = Uforward(:,:,ct)*Uback(:,:,ct)';
        tmpmult = conj(trace(tmpmat));
        tmpmatt = reshape(tmpmat,[1 2^(2*nbspins)]);
        
        for RFmattct = 1:1:size(RFmatts,3)
            derivs(ct,RFmattct+1) = 2*pulse(ct,1)*imag(sum(tmpRFmat(RFmattct,:).*tmpmatt)*tmpmult);
%            derivs(ct,RFmattct+1) = -2*real(trace(Uback(:,:,ct)'*i*pulse(ct,1)*RFmatts(:,:,RFmattct)*Uforward(:,:,ct))*trace(Uforward(:,:,ct)'*Uback(:,:,ct))); 
        end %RFmattct loop

        %If we want calculate the time step derivatives
        if(params.tstepflag)
            derivs(ct,1) = 2*imag(sum(reshape(transpose(HTOT(:,:,ct)),[1 2^(2*nbspins)]).*tmpmatt)*tmpmult);
        end

    end %ct loop

else %Do the same for state to state
    rhon = Uforward(:,:,plength)*rhoin*Uforward(:,:,plength)';
    for ct = 1:1:plength
        lambdaj = Uback(:,:,ct)*rhogoal*Uback(:,:,ct)';
        rhoj = Uforward(:,:,ct)*rhoin*Uforward(:,:,ct)';

        for RFmattct = 1:1:size(RFmatts,3)
            derivs(ct,RFmattct+1) = 2*pulse(ct,1)*imag(trace(lambdaj'*(RFmatts(:,:,RFmattct)*rhoj - rhoj*RFmatts(:,:,RFmattct)))*trace(rhon'*rhogoal));
        end %RFmattct loop

        if(tstepflag)
            derivs(ct,1) = 2*imag(trace(lambdaj'*(HTOT(:,:,ct)*rhoj - rhoj*HTOT(:,:,ct)))*trace(rhon'*rhogoal));
        end %tstepflag

    end %ct loop

end %params.searchtype if

%Calculate the trace squared fidelity of the propagator
if(params.searchtype == 1)
    goodness = (abs(trace(Uforward(:,:,plength)'*Uwant)))^2;
else
    goodness = (abs(trace(rhon'*rhogoal)))^2;
end

%If we are asking for the z derivatives calculate them with finite
%difference (fast enough because there isn't very many of them and
%the matrices are diagonal)

if(nargout == 3)
    zderivs = zeros(2,nbspins);

    %Switch depending on whether we are doing state to state or
    %unitary
    if(params.searchtype == 1)

        %Do the Zpre first
        tmpmat = diag(Uwant'*Uforward(:,:,plength));
        for ct = 1:1:nbspins
            newgoodness = (abs(sum(exp(-i*pi*1e-6*corrmat(:,ct)).*tmpmat)))^2;
            zderivs(1,ct) = 1e6*(newgoodness-goodness);
        end

        %And now Zpost
        tmpmat = diag(Uforward(:,:,plength)*Uwant');
        for ct = 1:1:nbspins
            newgoodness = (abs(sum(exp(-i*pi*1e-6*corrmat(:,ct)).*tmpmat)))^2;
            zderivs(2,ct) = 1e6*(newgoodness-goodness);
        end

    else %Do state to state
        %Do the Zpre first
        Usim = Uforward(:,:,plength);
        for ct = 1:1:nbspins
            newgoodness = (abs(trace(rhogoal'*Usim*diag(exp(-i*pi*1e-6*corrmat(:,ct)))*rhoin*diag(exp(i*pi*1e-6*corrmat(:,ct)))*Usim')))^2;
            zderivs(1,ct) = 1e6*(newgoodness-goodness);
        end

        %And now Zpost
        for ct = 1:1:nbspins
            newgoodness = (abs(trace(rhogoal'*diag(exp(-i*pi*1e-6*corrmat(:,ct)))*Usim*rhoin*Usim'*diag(exp(i*pi*1e-6*corrmat(:,ct))))))^2;
            zderivs(2,ct) = 1e6*(newgoodness-goodness);
        end
    end %params.searchtype if
end %nargout == 3 if