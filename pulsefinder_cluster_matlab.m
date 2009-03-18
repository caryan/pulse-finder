classdef pulsefinder_cluster < handle

    properties
        %IP address of the cluster with username (note need passwordless
        %ssh log in
        clusterHost = 'c4ryan@whale.sharcnet.ca';
        %Directory on remote host to work out of
        remotedata = '/work/c4ryan/pulsefinding';
        %Scheduler object to submit jobs to
        scheduler
        %Job object to hold the pulsefinding tasks
        currentjob;
        %Location of files needed on local computer
        neededfiles = {'~/NMR/Programs/pulsefinder/dev','~/NMR/Programs/auxiliary_files/'};
        %Location of working directory on the local computer
        localdata = '/Users/caryan/tmp';
        %Array of pulse numbers
        pulsenums = [];
        %Array of jobids
        jobids = [];
        %Array of pulsefinding tasks
        tasks;
        
    end

    methods
        %Constructor function
        function pf_cluster_obj = pulsefinder_cluster()

            %Setup the scheduler
            pf_cluster_obj.scheduler = findResource('scheduler', 'type', 'generic');
            set(pf_cluster_obj.scheduler,'DataLocation',pf_cluster_obj.localdata);
            set(pf_cluster_obj.scheduler, 'HasSharedFilesystem', true);
            set(pf_cluster_obj.scheduler, 'ClusterOsType', 'unix');
            set(pf_cluster_obj.scheduler, 'SubmitFcn', {@lsfNonSharedSimpleSubmitFcn, pf_cluster_obj.clusterHost, pf_cluster_obj.remotedata});

            % Use the correct one below for R2007b or R2008b
            % for R2007b_start
            set(pf_cluster_obj.scheduler, 'ClusterMatlabRoot', '/opt/sharcnet/matlab/current');
            % for R2007b_end

            % for R2008b_start
            % set(pf_cluster_obj.scheduler, 'ClusterMatlabRoot', '/opt/sharcnet/matlab/R2008b');
            % set(pf_cluster_obj.scheduler, 'DestroyJobFcn', @lsfDestroyJob);
            % set(pf_cluster_obj.scheduler, 'GetJobStateFcn', @lsfGetJobState);
            % for R2008b_end

            %Create a job object:
            pf_cluster_obj.currentjob = createJob(pf_cluster_obj.scheduler);
            set(pf_cluster_obj.currentjob,'FileDependencies',pf_cluster_obj.neededfiles);

        end %Constructor function

        %Method to add tasks
        function addpulse(pf_cluster_obj,paramsfile,pulsenum)
            pf_cluster_obj.pulsenums(end+1) = pulsenum;

            %Add the pulsefinding task
            createTask(pf_cluster_obj.currentjob, @pulsefinder_main, 1, {paramsfile,sprintf('%s/pulsefinder%d.log',pf_cluster_obj.remotedata,pulsenum)});
        
            %Update the tasks property
            pf_cluster_obj.tasks = get(pf_cluster_obj.currentjob,'Tasks');
            
            
        end %addpulse method

        %Method to submit job to cluster
        function submit(pf_cluster_obj)
    
            %Remove all the previous log files for these pulse numbers
            for pulsect = pf_cluster_obj.pulsenums
                unix(sprintf('rm -f %s/pulsefinder%d.log',pf_cluster_obj.localdata,pulsect));
                runCmdOnCluster(sprintf('rm -f %s/pulsefinder%d.log',pf_cluster_obj.remotedata,pulsect),pf_cluster_obj.clusterHost);
            end

            %Submit the job
            submit(pf_cluster_obj.currentjob);

            %Load the jobid numbers from the local log file
            [status,result] = unix(sprintf('tail -%d %s/cluster_stdout.log',2*length(pf_cluster_obj.pulsenums),tempdir));
            jobids = regexp(result,'Job <(\d+)> is submitted to queue <matlab>.','tokens');
            for regexpct = 1:1:length(pf_cluster_obj.pulsenums)
                pf_cluster_obj.jobids(pf_cluster_obj.pulsenums(regexpct)) = str2double(jobids{regexpct}{1});
            end

        end %submit method

        %Method to find the best result so far
        function dispresults(pf_cluster_obj)

            numrunning = 0;
            numfinished = 0;
            numqueued = 0;
            bestfidelity = 0;
            bestpulsenum = 0;
            
            for pulsect = pf_cluster_obj.pulsenums

                %Copy the log file back to the local machine
                unix(sprintf('scp -q %s:%s/pulsefinder%d.log %s/.',pf_cluster_obj.clusterHost,pf_cluster_obj.remotedata,pulsect,pf_cluster_obj.localdata));

                %See if the file is there
                tmplogname = sprintf('%s/pulsefinder%d.log',pf_cluster_obj.localdata,pulsect);
                if(exist(tmplogname,'file'))

                    %Use awk to load the final fidelity
                    [status,result] = unix(sprintf('awk ''/  0./{print $1};'' %s | tail -1',tmplogname));
                    curfidelity = str2double(result);
                    if(curfidelity > bestfidelity)
                        bestfidelity = curfidelity;
                        bestpulsenum = pulsect;
                    end
                end

                tasknum = find(pf_cluster_obj.pulsenums == pulsect);
                
                if(strcmp(get(pf_cluster_obj.tasks(tasknum),'State'),'queued'))
                    numqueued = numqueued +1;
                    fprintf('Pulse %d is queued and there are no results yet.\n',pulsect);
                elseif(strcmp(get(pf_cluster_obj.tasks(tasknum),'State'),'running'))
                    numrunning = numrunning + 1;
                    fprintf('Pulse %d is running and the current fidelity is %f.\n',pulsect,curfidelity);
                elseif(strcmp(get(pf_cluster_obj.currentjob,'State'),'finished'))
                    numfinished = numfinished + 1;
                    fprintf('Pulse %d is finished and the final fidelity is %f.\n',pulsect,curfidelity);
                end

            end

            disp(sprintf('%d pulses are running and pulse number %d is the best with a fidelity of %f.',numrunning,bestpulsenum,bestfidelity));

        end %bestresult method

        %Method to look at log files
        function displogfile(pf_cluster_obj,pulsenums)

            for pulsect = pulsenums
                tasknum = find(pf_cluster_obj.pulsenums == pulsect);
                disp(sprintf('Pulse number %d:',pulsect));
                %See if the job is even running yet
                if(strcmp(get(pf_cluster_obj.tasks(tasknum),'State'),'queued'))
                    disp('Job is still queued. No results yet.');
                else
                    %Copy the log file back to the local machine
                    unix(sprintf('scp -q %s:%s/pulsefinder%d.log %s/.',pf_cluster_obj.clusterHost,pf_cluster_obj.remotedata,pulsect,pf_cluster_obj.localdata));

                    %Display the file
                    type(sprintf('%s/pulsefinder%d.log',pf_cluster_obj.localdata,pulsect));
                end
            end
        end %displogfile method

        %Method to kill jobs on cluster
        %Currently has serious environment issues. There must be a nicer
        %way to get the same environment that the interactive login has
        function kill(pf_cluster_obj,pulsenums)
            tmpenv = 'export LSF_SERVERDIR=/opt/hptc/lsf/top/6.2/linux2.6-glibc2.3-x86_64/etc;export LSF_LIBDIR=/opt/hptc/lsf/top/6.2/linux2.6-glibc2.3-x86_64/lib;export LD_LIBRARY_PATH=/opt/hptc/lib:/opt/hptc/lsf/top/6.2/linux2.6-glibc2.3-x86_64/lib;export PATH=/opt/sharcnet/compile/bin:/opt/hptc/bin:/opt/hptc/lsf/top/6.2/linux2.6-glibc2.3-x86_64/etc:/opt/hptc/lsf/top/6.2/linux2.6-glibc2.3-x86_64/bin:/usr/kerberos/bin:/usr/bin:/bin:/usr/sbin:/sbin:/usr/X11R6/bin:/opt/sharcnet/amber/current/exe:/opt/sharcnet/archive_tools/public:/usr/local/pgsql/bin:/opt/sharcnet/local/gaussian/bin:/opt/sharcnet/gaussian/bin:/opt/sharcnet/gpc/current/bin:/opt/sharcnet/gsl/current/bin:/opt/sharcnet/lammps/current/bin:/opt/sharcnet/matlab/current/bin:/opt/sharcnet/mpiblast/current/bin:/opt/sharcnet/octave/current:/opt/sharcnet/octave/current/bin:/opt/sharcnet/pathscale/current/bin:/opt/sharcnet/pgi/pgi-6.1/linux86-64/6.1/bin:/opt/sharcnet/r/current/bin:/opt/sharcnet/sharutils/current/bin:/opt/sharcnet/sq/bin:/opt/sharcnet/nwchem/current/bin;export LSF_BINDIR=/opt/hptc/lsf/top/6.2/linux2.6-glibc2.3-x86_64/bin;export LSF_ENVDIR=/opt/hptc/lsf/top/conf';
            for pulsect = pulsenums
                tasknum = find(pf_cluster_obj.pulsenums == pulsect);
                unix(sprintf('ssh %s "%s;sqkill %d"',pf_cluster_obj.clusterHost,tmpenv,pf_cluster_obj.jobids(pulsect)));
                %Destroy the job
                destroy(pf_cluster_obj.tasks(tasknum));
            end

        end %kill method

        %Method to get result
        function results = getresults(pf_cluster_obj,pulsenums)
            results = cell(length(pulsenums));
            tmpct = 0;
            for pulsect = pulsenums
                tmpct = tmpct+1;
                tasknum = find(pf_cluster_obj.pulsenums == pulsect);
                tmpresult = get(pf_cluster_obj.tasks(tasknum),'OutputArguments');
                results{tmpct} = tmpresult{1};
            end

        end

    end %Methods section

end %classdef

