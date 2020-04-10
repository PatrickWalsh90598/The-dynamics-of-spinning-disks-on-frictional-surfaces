%% Lyapunov Exponent Plotter
%  Measuring the Lyapunov exponent involves taking the slope of the line of
%  best fit of a suitable set of linear datapoints near the beginning of
%  the plot. However, the system is extremely sensitive to the initial 
%  conditions and so the exact location of this set of datapoints varied 
%  erratically as theta was varied. This made the range of datapoints used 
%  for the LE calculation difficult to automate. The LE was therefore 
%  calculated manually for each value of Theta.
Theta = [28:89]';
LyapunovExponents = ...
    [72.6450
    70.6024
    68.5794
    72.0817
    77.2120
    85.5336
    91.7740
    98.0402
    105.2798
    107.0649
    112.8604
    122.9429
    129.5178
    136.8976
    148.8340
    152.5570
    167.0119
    173.2172
    163.9791
    153.0496
    143.2071
    136.7544
    126.0812
    122.2405
    120.3661
    113.6397
    109.9399
    107.6846
    104.8207
    97.6202
    94.6429
    93.2903
    92.2087
    89.2437
    85.3296
    81.9046
    81.5086
    78.3085
    77.7820
    77.8301
    78.7930
    79.6218
    81.5838
    83.7317
    86.7674
    88.5410
    95.8728
    99.0506
    107.5283
    109.9647
    113.7994
    127.4730
    147.8338
    164.9343
    186.3686
    205.6586
    271.4705
    275.8553
    379.4066
    631.5421
    890.1528
    2.3242e+03];

% The LE values are themselves highly erratic and are very sensitive to the
% exact range chosen to calculate them. Care was taken to ensure that the
% method for selecting the range was kept as consistent as possible for all
% angles. For this reason the exact values calculated for the Lyapunov 
% exponents may not be accurate and should not be taken as exact results
% but rather simply as an indication of the pattern of the relative 
% magnitude of the Lyapunov exponent with Theta. The values of the Lyapunov
% exponents are therefore all normalised against the minimum value.
NormVal = min(LyapunovExponents);
figure(1)
plot(Theta,smooth(LyapunovExponents/NormVal));
xlabel('\theta');
ylabel('\lambda');
savefig(figure(1),'Normalised Lyapunov Exponent vs Theta(Entire).fig')

