function pulseout = smoothpulse(pulsein,finaldisc,window,plotflag,phaserampflag)

%This function smoothes and rediscretizes a pulse.
%
% pulseout = smoothpulse(pulsein,findaldisc,window,plotflag,phaserampflag)
%
% pulsein: a pulse structure that you wish to smooth
% finaldisc: the final time discritization 
% window: the smoothing window 
% plotflag: 1 to plot, 0 otherwise
% phaserampflag: 1 to convert Z controls into phase ramps

%Load the pulse
chunkypulse = pulsein.pulse;

%Now discritize the pulse in 50ns periods
maxsize = fix(sum(chunkypulse(:,1))/50e-9 + 1/50);
discpulse = zeros(maxsize,size(chunkypulse,2)-1);
pointer = 1;
for ct = 1:1:size(chunkypulse,1)
  numrep = fix(chunkypulse(ct,1)/50e-9 + 1/50);
  
  discpulse(pointer:pointer+numrep-1,:) = repmat(chunkypulse(ct,2: ...
                                                    end),numrep,1);
  pointer = pointer+numrep;
end

%Remove any extra points
discpulse(pointer:end,:) = [];

%Sample the points
samprate = fix(finaldisc/50e-9 + 1/50);
samppulse = discpulse(fix(samprate/2):samprate:end,:);

%Smooth it out
smoothedpulse = filtfilt(ones(1,window)/window,1,samppulse);

%If we want, plot the result
if(plotflag)
  time = finaldisc*[1:1:size(samppulse,1)]';
  
  for ct = 1:1:size(samppulse,2)
  figure
  plot(time,samppulse(:,ct));
  hold on
  plot(time,smoothedpulse(:,ct),'r');
  end

end

%If we have the Z rf handle try to turn it into a phase ramp
if(phaserampflag == 1)
  pulseout = pulsein;
  for channelct = 1:1:2
    phasetrack = 0;
    for ct = 1:1:size(samppulse,1)
      amp = abs(samppulse(ct,3*(channelct-1)+1) + i*samppulse(ct,3*(channelct-1)+2));
      phase = angle(samppulse(ct,3*(channelct-1)+1) + i*samppulse(ct,3*(channelct-1)+2));

      phase = phase + phasetrack;
      
      samppulsebis(ct,2*(channelct-1)+1) = amp*cos(phase);
      samppulsebis(ct,2*(channelct-1)+2) = amp*sin(phase);
      
      phasetrack = phasetrack - finaldisc*samppulse(ct,3*(channelct-1)+3);
    end

   if(channelct == 1)
    pulseout.zangles(2,1:3) = pulseout.zangles(2,1:3)-phasetrack/2/pi;
   else
      pulseout.zangles(2,4:7) = pulseout.zangles(2,4:7)-phasetrack/2/pi;
  end
    
    
  end
  
  pulseout.pulse = [finaldisc*ones(size(samppulsebis,1),1) samppulsebis];
  pulseout.params.plength = size(pulseout.pulse,1);
  

  
else

pulseout = pulsein;
pulseout.pulse = [finaldisc*ones(size(samppulse,1),1) smoothedpulse];
pulseout.params.plength = size(pulseout.pulse,1);
end    
    

return
