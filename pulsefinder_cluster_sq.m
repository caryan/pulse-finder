classdef pulsefinder_cluster_sq < handle

    properties
        %IP address of the cluster with username (note need passwordless
        %ssh log in
        clusterHost = 'c4ryan@whale.sharcnet.ca';
        %Directory on remote host to work out of
        remotedata = '/work/c4ryan/pulsefinding';
        %Location of compiled puslefinder executable on cluster
        pulsefinder_exec = '/work/c4ryan/pulsefinding/pulsefinder';
        %Location of files needed on local computer
        neededfiles = {'~/NMR/Programs/pulsefinder/dev'};
        %Location of working directory on the local computer
        localdata = '/Users/caryan/tmp';
        %Array of pulse numbers
        pulsenums = [];
        %Cell array of pulsefinding tasks
        tasks = cell();
        
    end

    methods
        %Constructor function
        function pf_cluster_obj = pulsefinder_cluster_sq()


        end %Constructor function

        %Method to add pulsefinding jobs
        function addpulse(pf_cluster_obj,paramsfile,pulsenum,est_time)
            pf_cluster_obj.pulsenums(end+1) = pulsenum;

            %Add the pulsefinding task
            pf_cluster_obj.tasks{pulsenum}.paramsfile = paramsfile;
            
            %Set the status to loaded
            pf_cluster_obj.tasks{pulsenum}.status = 'loaded';
            
            %Set the estimated time
            pf_cluster_obj.tasks{pulsenum}.est_time = est_time;
            
        end %addpulse method

        %Method to submit jobs to cluster
        function submit(pf_cluster_obj,pulsenums)
    
            %Loop through each pulse
            for pulsect = pulsenums
                %Remove all the previous log files
                unix(sprintf('rm -f %s/pulsefinder%d.log',pf_cluster_obj.localdata,pulsect));
                pf_cluster_obj.runCmdOnCluster(sprintf('rm -f %s/pulsefinder%d.log',pf_cluster_obj.remotedata,pulsect),pf_cluster_obj.clusterHost);
                pf_cluster_obj.runCmdOnCluster(sprintf('rm -f %s/pulsefinder%d.log.mat',pf_cluster_obj.remotedata,pulsect),pf_cluster_obj.clusterHost);
            
            
                %Write the shell script to run the file
                
            
            end
           
            
            
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
        function kill(pf_cluster_obj,pulsenums)
            %Kill each pulse using the sqkill command
            for pulsect = pulsenums
                tasknum = find(pf_cluster_obj.pulsenums == pulsect);
                unix(sprintf('ssh %s "%s;sqkill %d"',pf_cluster_obj.clusterHost,tmpenv,pf_cluster_obj.jobids(pulsect)));
                %Destroy the job
                destroy(pf_cluster_obj.tasks(tasknum));
            end

        end %kill method

        %Method to get to most recent results
        function results = getresults(pf_cluster_obj,pulsenums)
            results = cell(length(pulsenums));
            tmpct = 0;
            for pulsect = pulsenums
                tmpct = tmpct+1;
                tasknum = find(pf_cluster_obj.pulsenums == pulsect);
                tmpresult = get(pf_cluster_obj.tasks(tasknum),'OutputArguments');
                results{tmpct} = tmpresult{1};
            end

        end %getresults method
        
        %Method to run a command on the cluster and write the response to
        %the log file (code mostly taken from matlab version)
        function runCmdOnCluster(pf_cluster_obj,command)
        
            %Use ssh to run the command
            cmdForCluster = sprintf('ssh %s ''%s''', pf_cluster_obj.clusterHost, command);

            [s, r] = unix(cmdForCluster);
            if s ~= 0
                error(['Failed to run the command\n' ...
                    '"%s"\n"' ...
                    'on the host "%s".\n' ...
                    'Command Output:\n' ...
                    '"%s"\n' ...
                    ], command, pf_cluster_obj.clusterHost, r);
            else
                fprintf('%s\n', r);

                %Also write it to the temporary log file
                tmpFID = fopen([tempdir 'cluster_stdout.log'],'a');
                fprintf(tmpFID,'%s\n',r);
                fclose(tmpFID);

            end
            
        end 


    end %Methods section

end %classdef

