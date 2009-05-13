function [goodness,derivs,zderivs] = pulsefinder_evalpulse(pulse,HNAT,RFmatts,Uwant_in,rhoin_in,rhogoal_in,zangles,params)

%Helper function to evaluate a pulse

%Modify the pulse by the tstepscale
pulse(:,1) = params.tstepscale*pulse(:,1);

%Load the rf and Hamiltonian distribution
rfdist = params.rfdist;
Hamdist = params.Hamdist;

%Initialize the goodness variable
goodness  = 0;

%Calculate the number of spins we are dealing with
nbspins = log2(size(HNAT,1));

%Initialize the derivatives
derivs = zeros(size(pulse));
zderivs = zeros(2,nbspins);

%Initialize a matrix for Hamiltonian storage
HTOT = [];

%Initialize the correction matrix for the Z rotations
corrmat = [];

%If we are allowing Zfreedom, then modify the Uwant or
%rhoin/rhogoal by the Z rotations
if(params.Zfreedomflag)
    
    %Set up the corrmat: each column contains the diagonal of ZIII,IZII....
    corrmat = zeros(2^nbspins,nbspins);
    for ct = 1:1:nbspins
        reps = 2^(ct-1);
        corrmat(:,end-ct+1) =  repmat([zeros(reps,1);ones(reps,1)],2^nbspins/reps/2,1);
    end
    corrmat  = -2*corrmat +1;
    
    Zpre = 0; Zpost = 0;
    for ct = 1:1:size(corrmat,2)
        Zpre = Zpre + zangles(1,ct)*corrmat(:,ct);
        Zpost = Zpost + zangles(2,ct)*corrmat(:,ct);
    end
    
end  %params.Zfreedomflag if

%Loop over the Hamiltonian distribution
for Hamct = 1:1:length(Hamdist)
    
    %Load the matrix of Hamiltonian shifts
    Hamshift = 2*pi*params.Hammatts{Hamct};
    
    %Reset the Uwant and rhoin and rhogoal to the params values
    if(params.searchtype == 1)
        Uwant = Uwant_in;
        rhoin = []; rhogoal = [];
    elseif(params.searchtype == 2)
        Uwant = [];
        rhoin = rhoin_in;
        rhogoal = rhogoal_in;
    end
    
    %Modify the desired unitary or the goal states by the soft pulse
    %buffer free evolution
    if(params.softpulsebuffer ~= 0)
        if(params.searchtype == 1)
            Uwant = expm(1i*params.softpulsebuffer*(HNAT+Hamshift))*Uwant*expm(1i*params.softpulsebuffer*(HNAT+Hamshift));
        elseif(params.searchtype == 2)
            rhoin = expm(-1i*params.softpulsebuffer*(HNAT+Hamshift))*rhoin* ...
                expm(+1i*params.softpulsebuffer*(HNAT+Hamshift));
            rhogoal = expm(1i*params.softpulsebuffer*(HNAT+Hamshift))*rhogoal* ...
                expm(-1i*params.softpulsebuffer*(HNAT+Hamshift));
        end
    end
    
    if(params.Zfreedomflag)
        if(params.searchtype == 1)
            Uwant = diag(exp(1i*pi*Zpost))*Uwant*diag(exp(1i*pi*Zpre));
        else
            rhoin = diag(exp(-1i*pi*Zpre))*rhoin*diag(exp(1i*pi*Zpre));
            rhogoal = diag(exp(1i*pi*Zpost))*rhogoal*diag(exp(-1i*pi*Zpost));
        end %params.searchtype if
    end
    
    %Loop over r.f. dist
    for rfct = 1:size(rfdist,1)
        
        %If we are doing a line search we only need the value of the
        %fitness function
        if(nargout==1)
            
            %Calculate the total propgator
            Usim = pulsefinder_calcprop(pulse,HNAT+Hamshift,RFmatts,rfdist(rfct,2),1);
            
            %Calculate the trace squared fidelity of the propagator
            if(params.searchtype == 1)
                goodness = goodness + Hamdist(Hamct)*rfdist(rfct,1)*(abs(trace(Usim'*Uwant)))^2;
            elseif(params.searchtype == 2)
                goodness = goodness + Hamdist(Hamct)*rfdist(rfct,1)*(abs(trace(Usim*rhoin*Usim'*rhogoal')))^2;
            end
            
            %Otherwise we need to store each timestep and calculate the derivatives
        else %nargout == 1 if
            
            %Calculate the Hamiltonian for each step and the propagator (if we
            %are doing timestep derivatives, we also need to store the HTOT for
            %each step
            if(params.tstepflag)
                [prop,HTOT] = pulsefinder_calcprop(pulse,HNAT+Hamshift,RFmatts,rfdist(rfct,2),0);
            else
                prop = pulsefinder_calcprop(pulse,HNAT+Hamshift,RFmatts,rfdist(rfct,2),0);
            end
            
            %Now calculate the derivatives
            if(params.Zfreedomflag)
                [tmpderivs,tmpgoodness,tmpzderivs] = pulsefinder_calcderivs(prop,Uwant,rhoin,rhogoal,RFmatts,pulse,HTOT,params,corrmat);
            else
                [tmpderivs,tmpgoodness] = pulsefinder_calcderivs(prop,Uwant,rhoin,rhogoal,RFmatts,pulse,HTOT,params,corrmat);
            end
            
            goodness = goodness + Hamdist(Hamct)*rfdist(rfct,1)*tmpgoodness;
            derivs = derivs + Hamdist(Hamct)*rfdist(rfct,1)*rfdist(rfct,2)*tmpderivs;
            
            if(params.Zfreedomflag)
                zderivs = zderivs + Hamdist(Hamct)*rfdist(rfct,1)*tmpzderivs;
            end
            
        end %nargout == 1 if
        
    end %rfct
end %Hamct loop

%Scale the goodness and derivs
if(params.searchtype == 1)
    goodness =  goodness/2^(2*nbspins);
    derivs = derivs/2^(2*nbspins);
    if(params.Zfreedomflag)
        zderivs = zderivs/2^(2*nbspins);
    end
else
    goodness = goodness/abs(trace(rhoin^2)*trace(rhogoal^2));
    derivs = derivs/abs(trace(rhoin^2)*trace(rhogoal^2));
    zderivs = zderivs/abs(trace(rhoin^2)*trace(rhogoal^2));
end

%Call the penalty function helper function
[goodness,derivs] = pulsefinder_penaltyfunctions(pulse,goodness,derivs,RFmatts,params);


