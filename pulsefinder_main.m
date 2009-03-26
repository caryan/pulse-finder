function pulses_out = pulsefinder_main(paramsfile,outputfilename)

%This MATLAB file tries to find strongly modulating pulses based on the
%GRAPE method see JMR Vol 172 pgs 296-305
%
% Usage:  [pulses] = pulsefinder(paramsfile (string));
% If you ctrl-c out of the function you can get at the results so far
% through the global variable pulses
%
%Written by Colm Ryan 26 February, 2006
%
%Updated by Colm Ryan 25 July, 2006 - to handle state to state,
%search with conjugate gradients, time step derivatives and
%handle all the different possibilities (Zfreedom,Hamdist,rfdist) nicely
%
%Updated by Colm Ryan 7 Aug,2006 to add option of on/off ramps
%
%
% Updated by Colm Ryan 12 January, 2007 to fix a small bug and
% change the averaging so that each Hamiltonian is averaged over
% the r.f.
%
% Updated by Colm Ryan 11 January, 2009 to spread to code out into
% different files to make it more readable. Also added possiblity to run
% on a cluster with a pulsefinder_cluster object.


%Define the pulses variable as global so that it can be retrieved after
%kill the program
global pulses
pulses = [];

%Load the output file or set to standard out
if(nargin == 1)
    outputFID = 1;
else
    outputFID = fopen(outputfilename,'w');
end

%Load the params file
params = pulsefinder_loadparams(paramsfile);

%Read the nucleus file and remove ignored spins
spins = read_nucleus_file(params.nucleifile);
if(isfield(params,'spins_ignore'))
    spins.ignore = params.spins_ignore;
    if(length(spins.ignore)>0)
        spins = spins_ignore(spins);
    end
end

%Update the spins.freqs to the pulsing frequencies
spins.freqs = spins.freqs - params.pulsefreq;

%Create the natural Hamiltonian in the pulsing reference frame
params.HNAT = full(CreateHamiltonian(spins));

%Set the total number of spins
params.nbspins = log2(size(params.HNAT,1));

%Sort out the subsystems
[Uwant_sub,rhoin_sub,rhogoal_sub,RFmatts_sub,HNAT_sub,params] = pulsefinder_subsystems(params,spins);

%Initialize some of the optimization variables
goodness = 0;
tryct = 0;

%Randomize the state of the random number generator
rand('state',sum(100*clock));

%Start the big optimization loop
while(goodness < params.fidelity && tryct < params.numtry)

    %Increment the try counter
    tryct = tryct+1;

    %Create a new pulse guess or load the starting point from
    %params.pulseguess
    [pulse,zangles,params] = pulsefinder_newpulse(params);

    %Print out a header line
    fprintf(outputFID,' Goodness |  Improvement/step |  Step Size  |  Iterations  |  CPU Time\n');

    %Initialize some checks and counters and other assorted variables
    stepsize = params.stepsize;
    oldgoodness = -1;
    improvechk = ones(1,20);
    dispchk = -1;
    iterct = 0;

    oldderivs = zeros(size(pulse));
    olddirec  = zeros(size(pulse));
    oldpulse = zeros(size(pulse));

    goodzdirec = zeros(2,spins.nb);
    oldzderivs = zeros(2,spins.nb);
    oldzdirec = zeros(2,spins.nb);
    oldzangles = zeros(2,spins.nb);

    betaresetct = 0;

    %Start the clock
    t = cputime;
    tic;
    lasttime = 0;

    %Initialize the pulses structure
    pulses{tryct}.params = params;

    %Save the initialguess
    pulses{tryct}.initialguess = pulse;

    %Start the small optimization loop
    while(stepsize > params.minstepsize && mean(improvechk) > params.improvechk && goodness< ...
            params.fidelity && betaresetct < 20)

        %Update the iterct
        iterct = iterct+1;

        %Reset the beta reset flag
        betaresetflag = 0;

        %Call a helper function which evaluates the goodness of the
        %current pulse and calculates the approximate derivatives
        goodness = 0;
        derivs = zeros(size(pulse));
        zderivs = zeros(2,spins.nb);
        if(params.Zfreedomflag)
            for ct = 1:1:length(params.subsystem)
                [tmpgoodness,tmpderivs,tmpzderivs] = pulsefinder_evalpulse(pulse,HNAT_sub{ct},RFmatts_sub{ct},Uwant_sub{ct},rhoin_sub{ct},rhogoal_sub{ct},zangles(:,params.subsystem{ct}),params);
                pulses{tryct}.subgood(ct) = tmpgoodness;
                goodness = goodness + params.subsys_weight(ct)*tmpgoodness;
                derivs = derivs + params.subsys_weight(ct)*tmpderivs;
                zderivs(:,params.subsystem{ct}) = zderivs(:,params.subsystem{ct}) + params.subsys_weight(ct)*tmpzderivs;
            end
        else
            for ct = 1:1:length(params.subsystem)
                [tmpgoodness,tmpderivs] = pulsefinder_evalpulse(pulse,HNAT_sub{ct},RFmatts_sub{ct},Uwant_sub{ct},rhoin_sub{ct},rhogoal_sub{ct},[],params);
                pulses{tryct}.subgood(ct) = tmpgoodness;
                goodness = goodness + params.subsys_weight(ct)*tmpgoodness;
                derivs = derivs + params.subsys_weight(ct)*tmpderivs;
            end
        end

        %Shift the derivatives by half a timestep (I'm not sure why I have to do
        %this but it takes ~1/1000 of the time to evalpulse and it makes the
        %derivatives almost perfect.  However, this method also
        %filters out high frequency components so I'm not
        %convinced it is the best way.
        derivs(:,2:end) = filter([0.5 0.5],1,derivs(:,2:end));
        derivs(1,2:end) = 2*derivs(1,2:end);

        %If we didn't improve then the conjugate gradient or our fitting
        %screwed up so go back to the old pulse and reset the beta
        if(goodness < oldgoodness)
            fprintf(outputFID,'Had to reset conjugate gradients!\n');

            pulse = oldpulse;
            derivs = oldderivs;
            goodness = oldgoodness;

            if(params.Zfreedomflag)
                zangles = oldzangles;
                zderivs = oldzderivs;
            end

            betaresetflag = 1;
            betaresetct = betaresetct+1;
        end

        %Call a helper function to calculate the conjugate gradient
        %direction and how far to move along it
        [gooddirec,goodzdirec,maxmult] = pulsefinder_conjgrad(pulse,zangles,derivs,oldderivs,zderivs,oldzderivs,olddirec,oldzdirec,betaresetflag,goodness,params,HNAT_sub,RFmatts_sub,Uwant_sub,rhoin_sub,rhogoal_sub,iterct,stepsize);

        %Update the improvement array for calculating the average
        %improvement
        improvechk = [goodness-oldgoodness improvechk(1:19)];

        %Transfer the variables to the old ones
        oldpulse = pulse;
        oldgoodness = goodness;
        oldderivs = derivs;
        olddirec = gooddirec;

        %Now move uphill and update the pulse
        pulse = pulse + maxmult*stepsize*gooddirec;

        if(params.Zfreedomflag)
            oldzangles = zangles;
            oldzderivs = zderivs;
            oldzdirec = goodzdirec;
            zangles = zangles + maxmult*stepsize*goodzdirec;
        end

        %Make sure we are not over %100 power
        pulse(:,2:end) = min(max(pulse(:,2:end),-1),1);

        %If we are changing the time steps remove those
        %that are negative
        if(params.tstepflag)
            zerolines = pulse(:,1)<0;
            if(sum(zerolines) ~= 0)
                pulse(zerolines,:) = [];
                derivs(zerolines,:) = [];
                gooddirec(zerolines,:) = [];
            end
        end

        %Change the stepsize by the sqrt of the maxmult so that we move in
        %the general direction of the right stepsize but don't jump around
        %too much
        stepsize = sqrt(maxmult)*stepsize;

        %Save the pulse in case we need to ctrl-c out of the
        %function
        tmppulse = pulse;
        tmppulse(:,1) = params.tstepscale*tmppulse(:,1);
        for ct = 1:1:size(params.RFmatts,3)
            tmppulse(:,ct+1) = params.rfmax(ct)*tmppulse(:,ct+1);
        end
        pulses{tryct}.pulse = tmppulse;
        pulses{tryct}.zangles = zangles;
        pulses{tryct}.goodness = goodness;

        %Now display if we want (basically we display when we are 20% closer
        %to the goal than the last time we displayed)
        if(goodness>dispchk || (cputime-lasttime > params.dispevery))
            fprintf(outputFID,'  %6.4f         %5.3e       %5.3e      %5d         %5.2f \n',goodness,mean(improvechk),stepsize,iterct,toc/60);

            %Set the next display trigger point
            dispchk = max((params.fidelity-goodness)*0.2,1e-4) + goodness;

            %Set the lasttime
            lasttime = cputime;
            
            %Save the result if we are logging to a file
            if(outputFID ~= 1)
                save([outputfilename '.mat'],'pulses');
            end
        end

    end %small optimization loop

    %Save each try because it might be as good as we get
    tmppulse = pulse;
    tmppulse(:,1) = params.tstepscale*tmppulse(:,1);
    for ct = 1:1:size(params.RFmatts,3)
        tmppulse(:,ct+1) = params.rfmax(ct)*tmppulse(:,ct+1);
    end
    pulses{tryct}.pulse = tmppulse;
    pulses{tryct}.zangles = zangles;
    pulses{tryct}.goodness = goodness;

    %Print out some info about why we gave up on this one.
    if(stepsize < params.minstepsize)
        fprintf(outputFID,'Failed to find pulse because stepsize is less than params.minstepsize.\n');
    end
    if(mean(improvechk < params.improvechk))
        fprintf(outputFID,'Failed to find pulse because no longer improving by params.improvechk.\n');
    end

    if(betaresetct > 19)
        fprintf(outputFID,'Failed to find pulse because had to reset conjugate gradients too many times. \n Usually this is because pulse is hitting maximum power limits\n');
    end

end %big optimization loop

%Print out some information about whether we got a good enough pulse or not
%and where the best one we found was
if(goodness>params.fidelity)
    fprintf(outputFID,'\nFound pulse of desired fidelity!\n');
else
    bestpulse = 1;
    for ct = 2:1:length(pulses)
        if(pulses{ct}.goodness > pulses{bestpulse}.goodness)
            bestpulse = ct;
        end
    end

    fprintf(outputFID,'\nUnable to find decent pulse.  Best one found is number %d at a fidelity of %f',bestpulse,pulses{bestpulse}.goodness);
    fprintf(outputFID,'You might want to try fiddling with the timestep, number of points or the random guess paramaters.\n');
end

if(outputFID ~= 1)
    fclose(outputFID);
end

%Return the pulses cell array
pulses_out = pulses;

