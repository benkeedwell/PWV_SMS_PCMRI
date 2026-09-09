function [offset, slope] = fitWeighted(xData, yData, weights)
  sumWeights = sum(weights);
  sumX  = sum(xData.*weights);
  sumY  = sum(yData.*weights);
  sumX2 = sum(xData.*xData.*weights);
  sumXY = sum(xData.*yData.*weights);
  %sumX2
  varX2 = sumX2*sumWeights - sumX*sumX;
  if (varX2 == 0)
    slope = 0;
  else
    slope = (sumXY*sumWeights - sumX*sumY)/varX2;
  end
  offset = (sumY-slope*sumX)/sumWeights;
end