# Wavelets_jl

Wavelets_jl is a set of functions that allow to perform the forward and inverse discrete Wavelet transform to 1-D signals. the code includes functions for wavelett based compression and denoising.

## Installation

download in the same folder the .jl scripts and the .csv files. From the main.jl file you can call all the functions, wavelets and test signals.

## Usage

You may use the functions from the wavelet general file by including the file in your script. However you need to read the csv_loading.jl file to load the coefficients of each wavelet transform.
The wavelets currently supported are:
Haar, daubechies2, daubechies4, daubechies5, daubechies6, daubechies8, daubechies 10, daubechies12, daubechies14 daubechies16, daubechies 18, daubechies20, coiflet5 and symlet5

The test signals that can be used are the Blocks, Bumps, Doppler and Heavy-sine functions.

## Contributing
Pull requests are welcome. For major changes, please open an issue first to discuss what you would like to change.

Please make sure to update tests as appropriate.

## License
This file is under GNU GPL license.
