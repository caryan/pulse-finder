function [prop,HTOT] = pulsefinder_calcprop(pulse,HNAT,RFmatts,rfmult,linesearchflag)


plength = size(pulse,1);
nbspins = log2(size(HNAT,1));

if(linesearchflag && nargout == 1)
    prop = eye(2^nbspins);
    for ct = 1:1:plength
        HTOT = HNAT;
        for RFmattct = 1:1:size(RFmatts,3)
            HTOT = HTOT + rfmult*pulse(ct,RFmattct+1)*RFmatts(:,:,RFmattct);
        end
        prop = expm(-i*pulse(ct,1)*HTOT)*prop;
    end

elseif(nargout == 2)
        %Initialize the propagator and Hamiltonian storage
        prop = zeros(2^nbspins,2^nbspins,plength);
        HTOT = repmat(HNAT,[1 1 plength]);
        for ct = 1:plength
            for RFmattct = 1:1:size(RFmatts,3)
                HTOT(:,:,ct) = HTOT(:,:,ct) + rfmult*pulse(ct,RFmattct+1)*RFmatts(:,:,RFmattct);
            end
            prop(:,:,ct) = expm(-i*pulse(ct,1)*HTOT(:,:,ct));
        end
else
        prop = zeros(2^nbspins,2^nbspins,plength);
        for ct = 1:plength
            HTOT = HNAT;
            for RFmattct = 1:1:size(RFmatts,3)
                HTOT = HTOT + rfmult*pulse(ct,RFmattct+1)*RFmatts(:,:,RFmattct);
            end
            prop(:,:,ct) = expm(-i*pulse(ct,1)*HTOT);
        end
end %linesearchflag nargout if
