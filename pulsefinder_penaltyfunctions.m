function [goodness,derivs] = pulsefinder_penaltyfunctions(pulse,goodness_in,derivs_in,RFmatts,params)

%Copy things over in case we don't do anything and/or to fill out the
%derivs
goodness = goodness_in;
derivs = derivs_in;

%Helper function to add any penalty functions.  Feel free to add whatever
%you want here but be careful that your derivatives match the goodness,
%otherwise the conjugate gradients will fail.  

%If we are doing timestep derivatives then modify the goodness and
%derivatives with the penalty function for extra time
if(params.tstepflag)
    tottime = sum(pulse(:,1));
    goodness = goodness_in - 0.01*exp(10*(tottime/params.tpulsemax - 1));

    if(params.searchtype == 1)
        derivs(:,1) = derivs_in(:,1) - 0.01*exp(-10)*(10/params.tpulsemax)*exp(10*tottime/params.tpulsemax);
    elseif(params.searchtype == 2)
        derivs(:,1) = derivs_in(:,1) - 0.01*exp(-10)*(10/params.tpulsemax)*exp(10*tottime/params.tpulsemax);
    end

    derivs(:,1) = derivs(:,1)*tstepscale;

end %tstepflag if

%If we want on and off ramps then modify the goodness and the
%derivatives with the penalty function
if(params.onofframps)
 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% These variables adjust the penalty function for the onoff ramps
% They are a bit of a hack for now however roughly speaking
% numpts says how many points at beginning/end are penalized
% param1 is the strength of the penalty (larger = stronger push to zero)
% param2 determines how quickly the penalty decreases/increases as
% you move away from the end points (larger = decreases faster)
% The penalty function is some sort of cosh function
  numpts = params.onoff_numpts;
  param1 = params.onoff_param1;
  param2 = params.onoff_param2;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  

  for ct = 1:1:size(RFmatts,3)
    %First the on ramp
    goodness = goodness_in - sum(param1*exp(-10*[param2:param2:param2*numpts]').*(exp(10*pulse(1:numpts,ct+1)) + exp(-10*pulse(1:numpts,ct+1))))...
                                     + 2*param1*sum(exp(-10*[param2:param2:param2*numpts]));
    
    %Now the off ramp
    goodness = goodness - sum(param1*exp(-10*[param2*numpts:-param2:param2]').*(exp(10*pulse(end-numpts+1:end,ct+1)) ...
                                  + exp(-10*pulse(end-numpts+1:end,ct+1)))) + 2*param1*sum(exp(-10*[param2:param2:param2*numpts]));

    %Now adjust the derivatives
    derivs(1:numpts,ct+1) = derivs_in(1:numpts,ct+1) - ...
        10*param1*(exp(-10*[param2:param2:param2*numpts]').*(exp(10*pulse(1:numpts,ct+1)) - exp(-10*pulse(1:numpts,ct+1))));
    
    derivs(end-numpts+1:end,ct+1) = derivs_in(end-numpts+1:end,ct+1) - 10*param1*exp(-10*[param2*numpts:-param2:param2]').* ...
        (exp(10*pulse(end-numpts+1:end,ct+1)) - exp(-10*pulse(end-numpts+1:end,ct+1)));
    
  end
end