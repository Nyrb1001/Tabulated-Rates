## Rates
The rates are sorted into folders by proton number, followed by files for each isotope. The folder structure is given by z{ZZZ}/z{ZZZ}n{NNN}.dat where Z and N are integers for the respective number of protons/neutrons. 

Each file is a wide-format grid. Each column represents a different temperature, and each row is for a different chemical potnential.

There are 29 temperature values (columns) in units of GK that are set as:  ``` T =      [0.0001, 0.0005, 0.001, 0.005, 0.01, 0.05, 0.1, 0.2, 0.25, 0.3, 0.4, &
    0.5, 0.6, 0.7, 0.8, 0.9, 1.0, 1.5, 2.0, 2.5, 3.0, 3.5, 4.0, 5.0, 6.0, &
    7.0, 8.0, 9.0, 10.0]  ``` 
    
There are also 89 Chemical potential values (rows) in units of MeV. The first 10 points range from (and including) 0 to 0.1 MeV with a spacing of 0.01 MeV, followed by the remaining 79 points with a spacing of 0.05 MeV until a maximum of 4.0 MeV. 

Each data entry is in %.6e formatting, with single spaces between columns. The units of the reaction rates are [$cm^3 mol^{-1} s^{-1}$]

## Routine

The control structure of the routine goes as follows:
- find the file corresponding to the n,z isotope identifier
- Retrieve a 4x4 grid of points surrounding the Temperature and chemical potential values you wish to interpolate
- the inner most 2x2 matrix is mainly used for the actual interpolation (bilinear interpolation), the outer values (ie indicies 1 and 4) are used for calculating the second derivative.
- Once the smaller matrix has been found and inteprolated, the error is estimated via the second derivative. For linear interpolation, the error is bounded by the second derivative within the two interpolant points.
- Because this is tabulated data, the second derivatives are calculated at each end point via finite differences, and the maximum is selected to estimate the interpolation error.
- return the rate and (upper bound) error.

I wasn't sure how you would interface/index your arrays so the read structure of this is pretty suboptimal. It could easily be changed by splitting the get-data() routine into two, with one routine storing the entire rate matrix, and another function call that references that stored matrix. Furthermore, the derivatives could be precalculated to save some compute. But, if you want something that works out of the box this should be fine.

A few notes:
- The error is significantly larger than what I'm getting when performing something such as an RMSE - these are the absolute worst case error values. In actuality, the interpolation is very close to the true value, especially once the chemical potential is > 0.1 MeV.
- I assume "clamped" end points on the second derivatives close to the tabulated limits - ie df/dx = 0. I *could* do a forward/backwards difference near the end but I don't think its going to make a huge difference.
- If important isotopes are missing, please let me know and I can add them asap. The turn around of this should be much faster now that I'm more confident in my calculations.
