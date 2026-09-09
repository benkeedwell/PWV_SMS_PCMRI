function smoothedData = smoothMS1(data, deg, m)
  if (nargin ~= 3)
    error("Usage: smoothMS1(dataRowVector, degree, m)");
  end
  %if (columns(data) < 2 || rows(data) ~= 1) %OCTAVE
  if (size(data,2) < 2 || size(data,1) ~= 1)
    error("Less than two data points or not a row vector");
  end
  kernel = kernelMS1(deg, m);
  fitWeights = edgeWeights1(deg, m);
  extData = extendData(data, m, fitWeights);
  smoothedExtData = conv(extData, kernel, "same");
  smoothedData = smoothedExtData(m+1 : end-m);
end
