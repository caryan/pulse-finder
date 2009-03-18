function [Uwant_sub,rhoin_sub,rhogoal_sub,RFmatts_sub,HNAT_sub,params] = pulsefinder_subsystems(params,spins)

%Helper pulsefinder function to handle the subsystem stuff

%First do some error checking
if(isempty(params.subsystem{1}))
    params.subsystem{1} = 1:1:spins.nb;
end

%Setup the RFmatts and HNAT for the subsystems
nbsub = length(params.subsystem);

HNAT_sub = cell(1,nbsub);
RFmatts_sub = cell(1,nbsub);

Uwant_sub = cell(1,nbsub);
rhoin_sub = cell(1,nbsub);
rhogoal_sub = cell(1,nbsub);

for ct = 1:1:nbsub
    %Sort the subsystem 
    params.subsystem{ct} = sort(params.subsystem{ct});
   
    %The spins we are taking out
    tspins = setxor(params.subsystem{ct},1:1:spins.nb);

    %The reduced natural Hamiltonian
    HNAT_sub{ct} = Partrace(params.HNAT,tspins)/2^length(tspins);

    %The reduced RFmatts
    for ct2 = 1:1:size(params.RFmatts,3)
        RFmatts_sub{ct}(:,:,ct2) =  params.rfmax(ct2)*Partrace(params.RFmatts(:,:,ct2),tspins)/2^length(tspins);
    end

    %The reduced Uwant 
    %To avoid a traceless unitary we apply a random unitary and then take the partial trace
    if(~iscell(params.Uwant))
    randU = 1;
    for ct3 = 1:1:spins.nb
      if(ismember(ct3,tspins))
        %Make a random 1 qubit unitary.  These aren't properly
        %distributed but oh well
        a = 0; b = 2*pi*rand; c = 2*pi*rand; d = 2*pi*rand;
        singleU = [exp(i*(a-b/2-d/2))*cos(c/2) -exp(i*(a-b/2+d/2))*sin(c/2); exp(i*(a+b/2-d/2))*sin(c/2) exp(i*(a+b/2+d/2))*cos(c/2)];
        randU = kron(randU,singleU);
      else
        randU = kron(randU,eye(2));
      end
    end
    if(params.searchtype == 1)
      Uwant_sub{ct} = Partrace(randU*params.Uwant,tspins);
      Uwant_sub{ct} = Uwant_sub{ct}/sqrt(trace(Uwant_sub{ct}*Uwant_sub{ct}'))*sqrt(2^length(params.subsystem{ct}));
    elseif(params.searchtype == 2)
      rhoin_sub{ct} = Partrace(randU*params.rhoin*randU',tspins);
      rhogoal_sub{ct} = Partrace(randU*params.rhogoal*randU',tspins);
    end %typeflag if
    
    else
      Uwant_sub{ct} = params.Uwant{ct};
    end %iscell if 
end %nbsub ct

%Normalize the sub-system weighting
params.subsys_weight = params.subsys_weight./sum(params.subsys_weight);
