function [pulse,zangles,params] = pulsefinder_newpulse(params)

%Load the pulse from the intial guess or create a new random guess
if(~isempty(params.pulseguess))
    disp(sprintf('\nLoading pulse from params.pulseguess'));
    pulse = params.pulseguess.pulse;
    
    %Scale the pulse with respect to the max power
    for ct = 1:1:size(pulse,2)-1
        pulse(:,ct+1) = pulse(:,ct+1)/params.rfmax(ct);
    end

    if(params.Zfreedomflag)
        zangles = params.pulseguess.zangles;
    else
        zangles = zeros(2,params.nbspins);
    end
    %There is only point in trying one
    params.numtry = 1;
else

    %Initialize the pulse i.e. create a skeleton of random points for each rf. field and then
    %fit the points with a interpolating cubic spline
    disp(sprintf('\nTrying new random guess....'));

    xold = [1:params.randevery:params.plength params.plength]';
    xnew = [1:1:params.plength]';
    skeletonpts = [];
    for ct = 1:1:size(params.RFmatts,3)
        skeletonpts = [skeletonpts params.randscale(ct)*(2*rand(size(xold,1),1)-1)];
    end
    skeletonpts(1,:) = 0; skeletonpts(end,:) = 0;
    pulse = interp1(xold,skeletonpts,xnew,'spline');
    
    %Add on the time 
    pulse = [params.timestep*ones(params.plength,1) pulse];

    %Initialize the zangles if we are allowing Zfreedom
    if(params.Zfreedomflag)
        zangles = rand(2,params.nbspins);
    else
        zangles = zeros(2,params.nbspins);
    end

end %load pulseguess if

%Redefine the timesteps in terms of the mean time step = 0.1 units
params.tstepscale = 10*mean(pulse(:,1));
pulse(:,1) = pulse(:,1)/params.tstepscale;