global pulseguess
global Uwantin

%Reset the params structure
params = [];

%Location of the nuclei file (full path name)
params.nucleifile = 'molecule_sample1.def';

%Number of timesteps
params.plength = 100;

%Length of each time step
params.timestep = 10e-6;

%Initial stepsize (this is reasonably important - run some tries
%and choose a value which is close to what the program is choosing
%after 50 or so iterations)
params.stepsize = 2e-2;

%Desired unitary
%H90
%params.Uwant = expm(-i*(pi/2)*(full(mkstate('+1IIX',0))));

Had = 1/sqrt(2)*[1 1;1 -1];
%Phase = [1 0;0 i];

I = eye(2);
params.Uwant = kron(kron(Had,I),I);

CNOT = [1 0 0 0;0 1 0 0;0 0 0 1;0 0 1 0];
%swap = [1 0 0 0;0 0 1 0;0 1 0 0;0 0 0 1];
%params.Uwant = kron(I,swap)*kron(CNOT,I)*kron(I,swap);

%params.Uwant = kron(I,CNOT);

%params.Uwant = Uwantin;

params.subsystem{1} = [1 2 3];

params.subsys_weight = [1];

%Desired fidelity for the unitary (this is the trace squared fidelity F = abs(Ugoal^dagger*Usim)^2/N^2)
params.fidelity = 0.999;

%RF distribution to optimize over (will slow down search and convergence dramatically)
%Two dimensional array first column is percentage of sample; second
%column is percentage of rf strength it sees. 
params.rfdist = [0.3 0.97;.4 1.00;0.3 1.03];

%params.rfdist = [1 1];

%Hamiltonian distribution to optimize over
%params.Hamdist = (1/8)*ones(1,8);

params.Hamdist = [1];

%Matrices for robustness to Hamiltonian distributions.  These will be
%multiplied by 2PI and added to the natural Hamiltonian.
params.Hammatts{1} = 0;

%Spins in nuclei file to ignore in search
params.spins_ignore = {};

%RF control fields (3D array of control Hamiltonians - as many as you like)
params.RFmatts(:,:,1) = (1/2)*full(mkstate('+1XII',0));
params.RFmatts(:,:,2) = (1/2)*full(mkstate('+1YII',0));
%params.RFmatts(:,:,3) = (1/2)*full(mkstate('+1ZII',0));
params.RFmatts(:,:,3) = (1/2)*full(mkstate('+1IXI+1IIX',0));
params.RFmatts(:,:,4) = (1/2)*full(mkstate('+1IYI+1IIY',0));
%params.RFmatts(:,:,6) = (1/2)*full(mkstate('+1IZI+1IIZ',0));


%The maximum rf power for each rf field in rad/s
params.rfmax = 2*pi*[25e3 25e3 16.7e3 16.7e3 15e3 15e3];

%Some parameters for the random guess 
%Scale of the random guess for each RFmatt (between 0 and 1)
params.randscale = [0.05 0.05 0.05 0.05 0.05 0.05];

%Choose every randevery points at random (rest will be fit to cubic spline)
params.randevery = 25;

%Tolerance for improving i.e. if over 20 tries we are not improving by an average of at least this, 
%we will try a new random starting point
params.improvechk = 1e-7;

%Minimum stepsize (if we're not moving anywhere we should stop searching)
params.minstepsize = 1e-7;

%Number of random guesses to try before giving up
params.numtry = 1;

%A vector of the pulsing frequency for each spin (this is defined
%with respect to the same frequency your nuclei file is)
params.pulsefreq = [0 0 0];	

%Soft pulse buffering delay (delays required before and after soft pulses.
%  (Since our time periods are usually greater than 350ns we could probably
    %try to use fast shapes)
params.softpulsebuffer = 10e-6;

%If there is a starting guess for the pulse load it in here (should
%be a structure
params.pulseguess = [pulseguess];

%The type of pulse we are searching for (1 for unitary, 2 for state
%to state)
params.searchtype = 1;

%Flag for whether you want to allow the time steps to vary
params.tstepflag = 0;

%Maximum length of pulse (rather soft boundary)
params.tpulsemax = 10e-3;

%Input and goal states for state to state
params.rhoin = mkstate('+1XII',1);
params.rhogoal = mkstate('+1YII',1);

%Allow Zfreedom or not
params.Zfreedomflag = 0;

%Parameter to force the beginnings and end of the pulse to zero
%with a penalty function
params.onofframps = 1;

%Display every x seconds irrespective of how things are improved
params.dispevery = 60;
