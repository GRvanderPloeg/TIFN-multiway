function Xclr=transformRCLR(X)
% Centered-log ratio transformation of compositional data
% Zero-handling implemented as suggested in Quinn et al 2019: 10.1093/gigascience/giz107
% Currently limited to two-way datasets
%
% Input:
% X     : input dataset
%
% Output:
% Xclr  : CLR transformed data

[I,J] = size(X);

% Set all zero values to NaN
idx = (X == 0);
X(idx) = NaN;

Xlong = X;
Xclr = Xlong;

for i=1:I
    dataVector = Xlong(i,:);
    Xclr(i,:) = log(dataVector/geomean(dataVector, "omitnan"));
end
