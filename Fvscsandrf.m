function goodness =  Fvscsandrf(pulsestruc,rfrange,Hammats,Hamrange)

%This function calculates and plots the r.f. and cs profile of a pulse
% function goodness =  Fvscsandrf(pulsestruc,rfrange,Hammats,Hamrange)
%
% Hamrange is optional axis for Hamiltonian distribution

%If Hamrange is not input then assign it. 
if(nargin < 4)
    Hamrange = 1:1:length(Hammats);
elseif(length(Hamrange) ~= length(Hammats))
    error('Length of Hamrange does not match Hammats');
end

%Load the parameters
params = pulsestruc.params;

pulse = pulsestruc.pulse;

%Read the nucleus file and remove ignored spins
spins = read_nucleus_file(params.nucleifile);
if(isfield(params,'spins_ignore'))
    spins.ignore = params.spins_ignore;
    if(~isempty(spins.ignore))
        spins = spins_ignore(spins);
    end
end

%Update the spins.freqs to the pulsing frequencies
spins.freqs = spins.freqs - params.pulsefreq;

%Load the control matts
RFmatts = params.RFmatts;

%Create the natural Hamiltonian in the pulsing reference frame
HNAT = full(CreateHamiltonian);

if(params.Zfreedomflag)
    
    %Set up the corrmat: each column contains the diagonal of ZIII,IZII....
    corrmat = zeros(2^nbspins,nbspins);
    for ct = 1:1:nbspins
        reps = 2^(ct-1);
        corrmat(:,end-ct+1) =  repmat([zeros(reps,1);ones(reps,1)],2^nbspins/reps/2,1);
    end
    corrmat  = -2*corrmat +1;
    
    Zpre = 0; Zpost = 0;
    for ct = 1:1:size(corrmat,2)
        Zpre = Zpre + pulsestruc.zangles(1,ct)*corrmat(:,ct);
        Zpost = Zpost + pulsestruc.zangles(2,ct)*corrmat(:,ct);
    end
    
end  %params.Zfreedomflag if

goodness.average = zeros(length(Hammats),length(rfrange));
goodness.worstcase = zeros(length(Hammats),length(rfrange));

%Double loop
fprintf('\n    ');
for Hamct = 1:1:length(Hammats)
    
    Hamshift = Hammats{Hamct};
    
    Uwant = params.Uwant;
    rhoin = params.rhoin;
    rhogoal = params.rhogoal;
    
    %Modify the desired unitary or the goal states by the soft pulse
    %buffer free evolution
    if(params.softpulsebuffer ~= 0)
        if(params.searchtype == 1)
            Uwant = expm(1i*params.softpulsebuffer*(HNAT+Hamshift))*Uwant*expm(1i*params.softpulsebuffer*(HNAT+Hamshift));
        elseif(params.searchtype == 2)
            rhoin = expm(-1i*params.softpulsebuffer*(HNAT+Hamshift))*rhoin* ...
                expm(+1i*params.softpulsebuffer*(HNAT+Hamshift));
            rhogoal = expm(1i*params.softpulsebuffer*(HNAT+Hamshift))*rhogoal* ...
                expm(-1i*params.softpulsebuffer*(HNAT+Hamshift));
        end
    end 
    
    if(params.Zfreedomflag)
        if(params.searchtype == 1)
            Uwant = diag(exp(1i*pi*Zpost))*Uwant*diag(exp(1i*pi*Zpre));
        else
            rhoin = diag(exp(-1i*pi*Zpre))*rhoin*diag(exp(1i*pi*Zpre));
            rhogoal = diag(exp(1i*pi*Zpost))*rhogoal*diag(exp(-1i*pi*Zpost));
        end %params.searchtype if
    end
    
    for rfct = 1:1:length(rfrange)
        
         %Calculate the total propgator
            Usim = pulsefinder_calcprop(pulse,HNAT+Hamshift,RFmatts,rfrange(rfct),1);
            
            %Calculate the trace squared fidelity of the propagator
            if(params.searchtype == 1)
                goodness.average(Hamct,rfct) = (abs(trace(Usim'*Uwant)))^2/(4^spins.nb);
                goodness.worstcase(Hamct,rfct) = wfidel(Usim,Uwant);
            elseif(params.searchtype == 2)
                goodness.average(Hamct,rfct) = (abs(trace(Usim*rhoin*Usim'*rhogoal')))^2/(4^spins.nb);
            end
            fprintf('\b\b\b\b%3.0f%%',100*((Hamct-1)*length(rfrange) + rfct)/(length(Hammats)*length(rfrange)));
    end
end
fprintf('\n');

%Now plot the results
figure

%If we only have rf distribution
if(length(Hammats) == 1)
    plot(rfrange,goodness.average);
    hold on
    plot(rfrange,goodness.worstcase,'r');
    legend('Average Fidelity','Worst-Case Fidelity');
    title('Fidelity versus RF Power Mulitplier');
    xlabel('RF Power multiplier');
    ylabel('Fidelity');
elseif(length(rfrange) == 1)
    plot(Hamrange,goodness.average);
    hold on
    plot(Hamrange,goodness.worstcase,'r');
    legend('Average Fidelity','Worst-Case Fidelity');
    title('Fidelity versus Hamiltonian Distribution')
    xlabel('Hamiltonian Distribution');
    ylabel('Fidelity');
else
   [C,H] = contour(rfrange,Hamrange,goodness.average,[0:0.1:0.9 0.91:0.01:1]);
   clabel(C,H);
end
   


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

%id = '$Id: wfidel.m,v 1.10 2006/06/19 20:03:15 msilva Exp $';

% some very basic error checking
su = size(u);
sv = size(v);
if su ~= sv
    error('MSILVA:BadMatrixSizes','Input matrices do not have the same size');
end;
% for debuging, maybe make sure that U and V are unitary?

% local variables
e  = eig(u'*v); % eigenvalues of u'*v
es = sortrows([e,angle(e)+pi],2);
d  = es(:,1);
f  = 0;

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
s = abs(ar - ar(mod((1:d)-2,d)+1));



