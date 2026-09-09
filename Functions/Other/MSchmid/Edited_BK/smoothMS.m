function smoothedData = smoothMS(data, deg, m)
  if (nargin ~= 3)
    error("Usage: smoothMS(dataRowVector, degree, m)");
  end
  %if (columns(data) < 2 || rows(data) ~= 1) %OCTAVE
  if (size(data,2) < 2 || size(data,1) ~= 1)
    error("Less than two data points or not a row vector");
  end
  kernel = kernelMS(deg, m);
  fitWeights = edgeWeights(deg, m);
  extData = extendData(data, m, fitWeights);
  smoothedExtData = conv(extData, kernel, "same");
  smoothedData = smoothedExtData(m+1 : end-m);
end
