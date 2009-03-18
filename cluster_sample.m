%Sample script for submitting jobs to sharcnet
pfobj = pulsefinder_cluster;

%Add three pulsefinding jobs with the same params file (could be different
%for varying parameters)
for ct = 1:1:3
    pfobj.addpulse('sampleparams1',ct);
end

%Submit the jobs
pfobj.submit;

%Check the results
pfobj.dispresults;

%Look at the details for pulses 1 and 3
pfobj.displogfile([1 3]);

%Kill pulsefinder 2
%pfobj.kill(2);

