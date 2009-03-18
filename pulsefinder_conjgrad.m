function [gooddirec,goodzdirec,maxmult] = pulsefinder_conjgrad(pulse,zangles,derivs,oldderivs,zderivs,oldzderivs,olddirec,oldzdirec,betaresetflag,goodness,params,HNAT_sub,RFmatts_sub,Uwant_sub,rhoin_sub,rhogoal_sub,iterct,stepsize_in)

%Helper function to calcuate the conjugate gradient direction and do a
%line search along that direction

%Calculate the conjugate gradient direction we should have to reset
%every once and a while but we don't seem to have to.
if(iterct ~= 1 && ~betaresetflag)
    diffderivs = derivs - oldderivs;
    if(params.Zfreedomflag)
        diffzderivs = zderivs - oldzderivs;
        beta = (sum(sum(derivs.*diffderivs))+sum(sum(zderivs.*diffzderivs)))/(sum(sum(oldderivs.^2)) + sum(sum(oldzderivs.^2)));
    else
        beta = sum(sum(derivs.*diffderivs))/sum(sum(oldderivs.^2));

    end

else
    beta = 0;
end

%Do a sort of reset.  If we have really lost conjugacy then beta will be
%negative.  If we than reset beta to zero then we start with the
%steepest descent again.
beta = max(beta,0);

%Define the good direction as the linear combination
gooddirec = derivs + beta*olddirec;
goodzdirec = zderivs + beta*oldzdirec;

%Now find the best step size along this direction by calculating two more points along this direction and fitting to a quadratic
mults = [0 1 2];
tmpgoodness = zeros(3,1);
tmpgoodness(1) = goodness;
for ct = 2:1:3
    for ct2 = 1:1:length(params.subsystem)
        subsys = params.subsystem{ct2};
        tmpgoodness_sub  =  pulsefinder_evalpulse(pulse+mults(ct)*stepsize_in*gooddirec,HNAT_sub{ct2},RFmatts_sub{ct2},Uwant_sub{ct2},rhoin_sub{ct2},rhogoal_sub{ct2},zangles(:,subsys)+mults(ct)*stepsize_in*goodzdirec(:,subsys),params);
        tmpgoodness(ct) = tmpgoodness(ct) + params.subsys_weight(ct2)*tmpgoodness_sub;
    end

end

%We have three points to fit a quadratic to.  The matrix to obtain
%the [a b c] coordinates for fitting points 0,1,2 is
fitcoeffs = [0.5   -1    0.5; -1.5    2   -0.5; ...
    1.0000         0         0]*tmpgoodness;

%If the quadratic is negative this method did not work so just
%go for the maximum value
if(fitcoeffs(1) > 0)
    [maximum,maxindex] = max(tmpgoodness);
    maxmult = mults(maxindex);
    %Otherwise choose the maximum of the quadratic
else
    maxmult = -fitcoeffs(2)/fitcoeffs(1)/2;
end

%If the max looks like it is beyond 2X the stepsize, try to fit up
%to 4X the stepsize using the same steps
if(maxmult > 1.99)
    mults = [0 2 4];
    tmpgoodness(2) = tmpgoodness(3);
    tmpgoodness(3) = 0;
    for ct2 = 1:1:length(params.subsystem)
        subsys = params.subsystem{ct2};
        tmpgoodness_sub  =  pulsefinder_evalpulse(pulse+mults(3)*stepsize_in*gooddirec,HNAT_sub{ct2},RFmatts_sub{ct2},Uwant_sub{ct2},rhoin_sub{ct2},rhogoal_sub{ct2},zangles(:,subsys)+mults(3)*stepsize_in*goodzdirec(:,subsys),params);
        tmpgoodness(3) = tmpgoodness(3) + params.subsys_weight(ct2)*tmpgoodness_sub;
    end

    fitcoeffs = [0.125   -0.25    0.125; -0.75    1   -0.25; ...
        1.0000         0         0]*tmpgoodness;

    if(fitcoeffs(1) > 0)
        [maximum,maxindex] = max(tmpgoodness);
        maxmult = mults(maxindex);
    else
        maxmult = -fitcoeffs(2)/fitcoeffs(1)/2;
    end

end %maxmult >2 if

%Move by at least 0.1 and at most 4;
maxmult = min(max(maxmult,0.1),4);



