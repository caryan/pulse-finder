function make_bruker_shape(pulsestrucin,calib,outputfile,channel)

%This function converts x-y amplitude from pulsefinder into amplitude and phase for spectrometer/simulator
%
%Call  make_bruker_shape(pulsestrucin,calib,outputfile,channel)
% pulse is the pulse you want to convert
% calib is the calibration value in Hz
% outputfile is the outputfilename
% channel is optional input for which channel you want to convert

if(nargin<4)
  channel = 1;
end

pulse = pulsestrucin.pulse(:,2*channel:2*channel+1);

%Some additional paramaters
%Bruker or soft flag
bruker = 0;

%How many pts per period to write out
ptsper = 1;

%User name
user = 'caryan';

shpfile = fopen(outputfile,'w');

amp =  abs(pulse(:,1)+i*pulse(:,2))*100/calib/2/pi;
phase = mod((180/pi)*angle(pulse(:,1)+i*pulse(:,2)),360);

if(bruker)
%First the starting stuff
%Write the first few lines
fprintf(shpfile,'##TITLE= %s\n',outputfile);
fprintf(shpfile,'##JCAMP-DX= 5.00 Bruker JCAMP library\n');
fprintf(shpfile,'##DATA TYPE= Shape Data\n');
fprintf(shpfile,'##ORIGIN= Colm''s GRAPE Pulses \n');
fprintf(shpfile,'##OWNER= %s\n',user);
fprintf(shpfile,'##DATE= %s\n',date);
time = clock;
fprintf(shpfile,'##TIME= %d:%d\n',fix(time(4)),fix(time(5)));
fprintf(shpfile,'##MINX= %7.6e\n',min(amp));
fprintf(shpfile,'##MAXX= %7.6e\n',max(amp));
fprintf(shpfile,'##MINY= %7.6e\n',min(phase));
fprintf(shpfile,'##MAXY= %7.6e\n',max(phase));
fprintf(shpfile,'##$SHAPE_EXMODE= None\n');
fprintf(shpfile,'##$SHAPE_TOTROT= %7.6e\n',90);
fprintf(shpfile,'##$SHAPE_BWFAC= %7.6e\n',1);
fprintf(shpfile,'##$SHAPE_INTEGFAC= %7.6e\n',1);
fprintf(shpfile,'##$SHAPE_MODE= 1\n');
fprintf(shpfile,'##NPOINTS= %d\n',length(amp));
fprintf(shpfile,'##XYPOINTS= (XY..XY)\n');
end

for ct = 1:1:length(amp)
     for ct2 = 1:1:ptsper
fprintf(shpfile,'  %7.6e,  %7.6e\n',amp(ct),phase(ct));
end
end


if(bruker)
     fprintf(shpfile,'##END=\n');
end

fclose(shpfile);

return
