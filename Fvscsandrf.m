function Fvscsandrf(pulsestruc);

global spins corrmat params

%This function calculates the r.f. and cs profile of a pulse

rfrange = [0.9:0.005:1.1];
csrange = [-10:1:10];

%rfrange = 1;
%csrange = 0;

%Load the parameters
params = pulsestruc.params;

pulse = pulsestruc.pulse;

Zfreedomflag = params.Zfreedomflag;
if(Zfreedomflag)
    zangles = pulsestruc.zangles;
end

%Read the nucleus file and remove ignored spins
spins = read_nucleus_file(params.nucleifile);
if(isfield(params,'spins_ignore'))
    spins.ignore = params.spins_ignore;
    if(length(spins.ignore)>0)
        spins = spins_ignore(spins);
    end
end

%Update the spins.freqs to the pulsing frequencies
spins.freqs = spins.freqs - params.pulsefreq;

%Load the control matts
RFmatts = params.RFmatts;

%Create the natural Hamiltonian in the pulsing reference frame
HNAT = full(CreateHamiltonian);

%If we allow Z freedom calculate a matrix of the corrections (each
%column contains the diagonal of ZIII,IZII....
if(Zfreedomflag)
    corrmat = zeros(2^spins.nb,spins.nb);
    for ct = 1:1:2^spins.nb
        tmpstr = dec2bin(ct-1,spins.nb);
        for ct2 = 1:1:spins.nb
            corrmat(ct,ct2) = str2num(tmpstr(ct2));
        end
    end
    corrmat = -2*corrmat + 1;
end

%Now first do the rfdist
fprintf('Calculating fidelity over r.f. profile........');
for rfct = 1:1:length(rfrange);

    %Calculate the goodness
    if(Zfreedomflag)
        [rfgoodness(rfct),rfworstcase(rfct)] = evalpulse(pulse,HNAT,rfrange(rfct)*RFmatts,zangles);
    else
        [rfgoodness(rfct),rfworstcase(rfct)] = evalpulse(pulse,HNAT,rfrange(rfct)*RFmatts);
    end
    %Display some info about how much we have finished
    fprintf('\b\b\b\b\b%3.0f%% ',100*rfct/length(rfrange));
end

%Now the csdist
fprintf('\nCalculating fidelity over chemical shift profile........');
for csct = 1:1:length(csrange);

    %Modify the Hamiltonian
    spins.freqs = spins.freqs + csrange(csct);

    HNAT = full(CreateHamiltonian);

    %Calculate the goodness
    if(Zfreedomflag)
        [csgoodness(csct),csworstcase(csct)] = evalpulse(pulse,HNAT,RFmatts,zangles);
    else
        [csgoodness(csct),csworstcase(csct)] = evalpulse(pulse,HNAT,RFmatts);
    end

    spins.freqs = spins.freqs - csrange(csct);
    fprintf('\b\b\b\b%3.0f%%',100*csct/length(csrange));
end

fprintf('\n');

%Now plot the results
figure
subplot(2,1,1)
plot(rfrange,rfgoodness);
if(params.searchtype == 1)
    hold on
    plot(rfrange,rfworstcase,'r');
    legend('Average Fidelity','Worst Case Fidelity');
end

xlim([rfrange(1) rfrange(end)]);
xlabel('RF Power Multiplier'); ylabel('Fidelity'); title('Fidelity over RF Distribution');
subplot(2,1,2);
plot(csrange,csgoodness);
if(params.searchtype == 1)
    hold on
    plot(csrange,csworstcase,'r');
    legend('Average Fidelity','Worst Case Fidelity');
end
xlabel('Chemical Shift in Hz'); ylabel('Fidelity'); title('Fidelity over Chemical Shift Distribution');


return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [goodness,worstcase] = evalpulse(pulse,HNAT,RFmatts,zangles);

%Subfunction to cacluate goodness of pulse and approximate derivatives
global spins corrmat
global params

plength = size(pulse,1);

%Calculate the Hamiltonian for each step and the propagator
prop = eye(2^spins.nb);
for ct = 1:1:plength
    HTOT = HNAT;
    for RFmattct = 1:1:size(RFmatts,3)
        HTOT = HTOT + pulse(ct,RFmattct+1)*RFmatts(:,:,RFmattct);
    end
    prop = expm(-i*pulse(ct,1)*HTOT)*prop;
end

%If we have allowed Zfreedom, include that
if(nargin == 4)
    Zpre = 0; Zpost = 0;
    for ct = 1:1:size(corrmat,2)
        Zpre = Zpre + zangles(1,ct)*corrmat(:,ct);
        Zpost = Zpost + zangles(2,ct)*corrmat(:,ct);
    end

    prop = diag(exp(-i*pi*Zpost))*prop*diag(exp(-i*pi*Zpre));
end

%Modify by the soft pulse buffer
prop = expm(-i*params.softpulsebuffer*HNAT)*prop*expm(-i*params.softpulsebuffer*HNAT);


%Calculate the fidelity of the propagator

if(params.searchtype == 1)

    goodness = (abs(trace(prop'*params.Uwant)))^2;
    worstcase = wfidel(prop,params.Uwant);

    %Scale the propagator to between 0 and 1
    goodness = goodness/2^(2*spins.nb);
else

    %Or state fidelity for a certain input state
    goodness = (abs(trace(params.rhogoal'*prop*params.rhoin*prop')))^2;
    worstcase = 0;
end


return


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function f = wfidel(u, v)
% WFIDEL worst case pure state fidelity between two unitaries.
%
%    April 6, 2006. M. Silva : m silva at iqc dot ca
%                   Join work with D. Kribs
%
%    For two unitaries U and V, WFIDEL(U,V) is the worst case
%    fidelity between U*a and V*a, where a is a complex vector
%    with norm 1. That is, WFIDEL(U,V) is equal to the minimum
%    of a'*U'*V*a over all complex vectors a with unit norm.
%
%    This implementation is, as of 6/4/6, untested, but the
%    theory is sound.
%
%    As of 15/6/6, it has been extensively tested for unitary
%      and valid input, and it seems solid.
%
%    $Id: wfidel.m,v 1.10 2006/06/19 20:03:15 msilva Exp $
%
% Requires
%    isnear.m
%    (from
%    http://www.mathworks.com/matlabcentral/fileexchange/loadFile.do?objectId=7098&objectType=FILE )
%

id = '$Id: wfidel.m,v 1.10 2006/06/19 20:03:15 msilva Exp $';

% some very basic error checking
su = size(u);
sv = size(v);
if su ~= sv
    error('MSILVA:BadMatrixSizes','Input matrices do not have the same size');
end;
% for debuging, maybe make sure that U and V are unitary?

% local variables
e  = eig(u'*v); % eigenvalues of u'*v
es = sortrows([e,angle(e)+pi],[2]);
d  = es(:,1);
f  = 0;

es(:,2)=es(:,2)-pi;
n = length(d);

% find the minimum distance from the convex hull
% of the distinct eigenvalues to the origin. if
% the origin is contained in the convex hull,
% that distance is defined as 0.
if n==1
    f = 1;
else
    if n==2
        f = orig_to_line(d(1),d(2));
    else
        dn = dist_to_neighbour(angle(d));
        dn(1) = 2*pi-sum(dn(2:n)); % the sum of the angular separations is 2*pi
        % and the boundary cases are funny
        % so it is best to calc it this way
        dn = find(dn > pi);
        if length(dn>0)==1
            f = orig_to_line(d(dn),d(mod(dn-2,n)+1));
        end;
    end;
end;

%
% Helper functions
%
function d = orig_to_line( c1, c2 )
% straight out of http://mathworld.wolfram.com/Point-LineDistance2-Dimensional.html
x1 = real(c1);
y1 = imag(c1);
x2 = real(c2);
y2 = imag(c2);
d = abs((x2-x1)*y1-x1*(y2-y1))/sqrt((x2-x1)^2+(y2-y1)^2);

function s = dist_to_neighbour( ar )
d = length(ar);
s = abs(ar - ar(mod([1:d]-2,d)+1));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function tf=isnear(a,b,tol)
%ISNEAR True Where Nearly Equal.
% ISNEAR(A,B) returns a logical array the same size as A and B that is True
% where A and B are almost equal to each other and False where they are not.
% A and B must be the same size or one can be a scalar.
% ISNEAR(A,B,TOL) uses the tolerance TOL to determine nearness. In this
% case, TOL can be a scalar or an array the same size as A and B.
%
% When TOL is not provided, TOL = SQRT(eps).
%
% Use this function instead of A==B when A and B contain noninteger values.

% D.C. Hanselman, University of Maine, Orono, ME 04469
% Mastering MATLAB 7
% 2005-03-09

%--------------------------------------------------------------------------
if nargin==2
    tol=sqrt(eps);
end
if ~isnumeric(a) | isempty(a) | ~isnumeric(b) | isempty(b) |...
        ~isnumeric(tol) | isempty(tol)
    error('Inputs Must be Numeric.')
end
if any(size(a)~=size(b)) & numel(a)>1 & numel(b)>1
    error('A and B Must be the Same Size or Either can be a Scalar.')
end
tf=abs((a-b))<=abs(tol);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


return

