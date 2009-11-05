%Sample script for submitting jobs to Debbie's cluster
pfobj = pulsefinder_cluster_sge;

%Add three pulsefinding jobs with the same params file (could be different
%for varying parameters)
for ct = 1:1:3
    pfobj.addpulse('sampleparams1.m',ct);
end

%Submit the jobs
pfobj.submit(1:3);

%Check the results
pfobj.dispresults;

%Look at the details for pulses 1 and 3
pfobj.displogfile([1 3]);

%Get the current pulses for task number 2
pulse2 = pfobj.getresults(2)

%Kill pulsefinder 2
%pfobj.kill(2);

