%FFT along specified dimensions
function ksp = FFTdim(img, dims)
    for i = dims
        img = (1/sqrt(size(img,i))) * fftshift(fft(ifftshift(img,i),[],i),i);
    end
    ksp = img;
end
