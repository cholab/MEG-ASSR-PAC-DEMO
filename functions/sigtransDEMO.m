function [Hamp,Hphase] = sigtransDEMO(ampdata,phasedata)
Hamp = abs(hilbert(permute(ampdata,[2,1,3]))).^2;
Hphase = angle(hilbert(permute(phasedata,[2,1,3])));
Hamp = permute(Hamp,[2,1,3]);
Hphase = permute(Hphase,[2,1,3]);
end