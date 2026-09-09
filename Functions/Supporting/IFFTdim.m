%IFFT along specified dimensions
function img = IFFTdim(ksp,dims)
    for i = dims
        ksp = sqrt(size(ksp,i)) * fftshift(ifft(ifftshift(ksp,i),[],i),i);
    end
    img = ksp;
end