classdef pulsefinder_cluster_sq < handle

    properties
        %IP address of the cluster with username (note need passwordless
        %ssh log in)
        clusterHost = 'c4ryan@whale.sharcnet.ca';
        %Directory on remote host to work out of
        remotedata = '/work/c4ryan/pulsefinding';
        %Location of compiled puslefinder executable on cluster
        pulsefinder_exec = '/work/c4ryan/pulsefinding/pulsefinder_whale';
        %Location of working directory on the local computer
        localdata = '/Users/caryan/tmp';
        %Array of pulse numbers
        pulsenums = [];
        %Array of pulsefinding task structures
        tasks = [];
        %Location of MCR on cluster
        MCRROOT = '/opt/sharcnet/local/matlab/matlab_2007b/';

    end

    methods
        %Constructor function
        function pf_cluster_obj = pulsefinder_cluster_sq()


        end %Constructor function

        %Method to add pulsefinding jobs
        function addpulse(pf_cluster_obj,paramsfile,pulsenum,est_time)
            pf_cluster_obj.pulsenums(end+1) = pulsenum;

            %Add the pulsefinding task
            pf_cluster_obj.tasks(pulsenum).paramsfile = paramsfile;

            %Set the status to loaded
            pf_cluster_obj.tasks(pulsenum).status = 'loaded';

            %Set the estimated time
            pf_cluster_obj.tasks(pulsenum).est_time = est_time;

        end %addpulse method

        %Method to submit jobs to cluster
        function submit(pf_cluster_obj,pulsenums)

            %Loop through each pulse
            for pulsect = pulsenums
                %Remove all the previous log files
                unix(sprintf('rm -f %s/pulsefinder%d.log',pf_cluster_obj.localdata,pulsect));
                pf_cluster_obj.runCmdOnCluster(sprintf('rm -f %s/pulsefinder%d.log %s/pulsefinder%d.mat.log %s/run_pulsefinder%d.sh',pf_cluster_obj.remotedata,pulsect,pf_cluster_obj.remotedata,pulsect,pf_cluster_obj.remotedata,pulsect));

                %Write the shell script to run the file
                pf_cluster_obj.writescript(pulsect);

                %Copy the shell scipt and params file over
                unix(sprintf('scp -q %s/run_pulsefinder%d.sh %s:%s/.',pf_cluster_obj.localdata,pulsect,pf_cluster_obj.clusterHost,pf_cluster_obj.remotedata));
                unix(sprintf('scp -q %s.m %s:%s/.',pf_cluster_obj.tasks(pulsect).paramsfile,pf_cluster_obj.clusterHost,pf_cluster_obj.remotedata));

                %Submit the job using sqsub
                sqsub_result = pf_cluster_obj.runCmdOnCluster(sprintf('bash --login -c "sqsub -q serial -r %s %s/run_pulsefinder%d.sh"',pf_cluster_obj.tasks(pulsect).est_time,pf_cluster_obj.remotedata,pulsect));

                %Load the jobid numbers from the local log file
                tmpjobid = regexp(sqsub_result,'submitted as jobid (\d+)','tokens','once');
                pf_cluster_obj.tasks(pulsect).jobid = str2double(tmpjobid{1});

                %Update the task status
                pf_cluster_obj.tasks(pulsect).status = 'submitted';

            end

        end %submit method

        %Method to find the best result so far
        function dispresults(pf_cluster_obj)

            %Update the status of the pulses
            pf_cluster_obj.updatestatus;

            bestfidelity = 0;
            bestpulsenum = 0;

            for pulsect = pf_cluster_obj.pulsenums


                if(strcmp(pf_cluster_obj.tasks(pulsect).status,'queued'))
                    fprintf('Pulse %d is queued and there are no results yet.\n',pulsect);
                else
                    curfidelity = 0;
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


                    if(strcmp(pf_cluster_obj.tasks(pulsect).status,'running'))
                        fprintf('Pulse %d is running and the current fidelity is %f.\n',pulsect,curfidelity);
                    elseif(strcmp(pf_cluster_obj.tasks(pulsect).status,'finished'))
                        fprintf('Pulse %d is finished and the final fidelity is %f.\n',pulsect,curfidelity);
                    end
                end

            end

            if(bestpulsenum == 0)
                disp('No results yet.');
            else
                disp(sprintf('Pulse number %d is the best with a fidelity of %f.',bestpulsenum,bestfidelity));
            end
        end %bestresult method

        %Method to look at log files
        function displogfile(pf_cluster_obj,pulsenums)

            %Update the status of the pulses
            pf_cluster_obj.updatestatus;

            for pulsect = pulsenums
                disp(sprintf('Pulse number %d:',pulsect));
                %See if the job is even running yet
                if(strcmp(pf_cluster_obj.tasks(pulsect).status,'queued'))
                    disp('Job is still queued. No results yet.');
                elseif(strcmp(pf_cluster_obj.tasks(pulsect).status,'loaded'))
                    disp('Job is not submitted yet.');
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
                %Use sqkill to kill the job if it is running
                if(strcmp(pf_cluster_obj.tasks(pulsect).status,'running') || strcmp(pf_cluster_obj.tasks(pulsect).status,'queued'))
                    %Use sqkill to kill the job
                    pf_cluster_obj.runCmdOnCluster(sprintf('bash --login -c "sqkill %d"',pf_cluster_obj.tasks(pulsect).jobid));
                    %Change the status
                    pf_cluster_obj.tasks(pulsect).status = 'killed';
                    disp(sprintf('Killed pulse number %d.',pulsect));
                end
            end

        end %kill method

        %Method to get to most recent results
        function results = getresults(pf_cluster_obj,pulsenums)
            results = cell(length(pulsenums));
            tmpct = 0;
            for pulsect = pulsenums
                tmpct = tmpct+1;

                %Copy the mat file back
                unix(sprintf('scp -q %s:%s/pulsefinder%d.log.mat %s/.',pf_cluster_obj.clusterHost,pf_cluster_obj.remotedata,pulsect,pf_cluster_obj.localdata));

                %Load the file
                tmppulse = load(sprintf('%s/pulsefinder%d.log.mat',pf_cluster_obj.localdata,pulsect));
                results{tmpct} = tmppulse.pulses;

            end

        end %getresults method

        %Method to run a command on the cluster and write the response to
        %the log file (code mostly taken from matlab version)
        function r = runCmdOnCluster(pf_cluster_obj,command)

            %Use ssh to run the command
            cmdForCluster = sprintf('ssh %s ''%s''', pf_cluster_obj.clusterHost, command);

            [s, r] = unix(cmdForCluster);
            if s ~= 0
                warning(['Failed to run the command\n' ...
                    '"%s"\n"' ...
                    'on the host "%s".\n' ...
                    'Command Output:\n' ...
                    '"%s"\n' ...
                    ], command, pf_cluster_obj.clusterHost, r);
            end

        end %runCmdOnCluster method

        %Method to write the script to correctly setup the environment and
        %run the pulsefinder
        function writescript(pf_cluster_obj,pulsenum)
            %Open the file
            scriptFID = fopen(sprintf('%s/run_pulsefinder%d.sh',pf_cluster_obj.localdata,pulsenum),'w');

            %Write the opening lines
            fprintf(scriptFID,'#/bin/sh\n');
            fprintf(scriptFID,'#Script to run a pulsefinding job\n\n');
            %Setup the environment
            fprintf(scriptFID,'MCRROOT=%s\n',pf_cluster_obj.MCRROOT);
            fprintf(scriptFID,'\n');
            fprintf(scriptFID,'LD_LIBRARY_PATH=.:${MCRROOT}/runtime/glnxa64;\n');
            fprintf(scriptFID,'LD_LIBRARY_PATH=${LD_LIBRARY_PATH}:${MCRROOT}/bin/glnxa64;\n');
            fprintf(scriptFID,'LD_LIBRARY_PATH=${LD_LIBRARY_PATH}:${MCRROOT}/sys/os/glnxa64;\n');
            fprintf(scriptFID,'\n');
            fprintf(scriptFID,'MCRJREVER=`cat ${MCRROOT}/sys/java/jre/glnxa64/jre.cfg`;\n');
            fprintf(scriptFID,'MCRJRE=${MCRROOT}/sys/java/jre/glnxa64/jre${MCRJREVER}/lib/amd64;\n');
            fprintf(scriptFID,'LD_LIBRARY_PATH=${LD_LIBRARY_PATH}:${MCRJRE}/native_threads;\n');
            fprintf(scriptFID,'LD_LIBRARY_PATH=${LD_LIBRARY_PATH}:${MCRJRE}/server;\n');
            fprintf(scriptFID,'LD_LIBRARY_PATH=${LD_LIBRARY_PATH}:${MCRJRE}/client;\n');
            fprintf(scriptFID,'LD_LIBRARY_PATH=${LD_LIBRARY_PATH}:${MCRJRE};\n');
            fprintf(scriptFID,'\n');
            fprintf(scriptFID,'XAPPLRESDIR=${MCRROOT}/X11/app-defaults;\n');
            fprintf(scriptFID,'\n');
            fprintf(scriptFID,'export LD_LIBRARY_PATH;\n');
            fprintf(scriptFID,'export XAPPLRESDIR;\n');
            fprintf(scriptFID,'\n');

            %Run the pulseprogram
            fprintf(scriptFID,'%s %s/%s %s/pulsefinder%d.log\n',pf_cluster_obj.pulsefinder_exec,pf_cluster_obj.remotedata,pf_cluster_obj.tasks(pulsenum).paramsfile,pf_cluster_obj.remotedata,pulsenum);
            fprintf(scriptFID,'\n');
            fprintf(scriptFID,'exit\n');

            %Close the file
            fclose(scriptFID);

            %Make the script executable
            unix(sprintf('chmod u+x %s/run_pulsefinder%d.sh',pf_cluster_obj.localdata,pulsenum));

        end %writescript method

        %Method to update the status of pulses
        function updatestatus(pf_cluster_obj)

            %Use sqjobs to get the status of submitted jobs
            tmpstr = pf_cluster_obj.runCmdOnCluster('bash --login -c "sqjobs"');
            tokens = regexp(tmpstr,'(\d+)\s+s\s+(R|Q)\s+\d\s+(\d+)\s+(\d+\w)','tokens');

            jobids = zeros(size(tokens,2),1);
            for ct = 1:1:size(tokens,2)
                jobids(ct) = str2double(tokens{ct}{1});
            end

            %Loop through the jobs
            for pulsect = pf_cluster_obj.pulsenums
                %Find which line it corresponds to
                linenum = find(jobids == pf_cluster_obj.tasks(pulsect).jobid);

                curstatus = pf_cluster_obj.tasks(pulsect).status;
                %If it wasn't there and was previously running then assume
                %it is finished
                if(isempty(linenum))
                    if(strcmp(curstatus,'running') || strcmp(curstatus,'queued'))
                        pf_cluster_obj.tasks(pulsect).status = 'finished';
                    end
                    %Otherwise update the status
                elseif(strcmp(tokens{linenum}{2},'R'))
                    pf_cluster_obj.tasks(pulsect).status = 'running';
                elseif(strcmp(tokens{linenum}{2},'Q'))
                    pf_cluster_obj.tasks(pulsect).status = 'queued';
                end

            end


        end %updatestaus method


    end %Methods section

end %classdef

