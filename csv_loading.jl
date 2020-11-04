using DataFrames,CSV
### To invert the coefficients and change the sign
function wav_coeff(value_filters,n,pos)
    ans = zeros(n,2)
    ans[:,1] = value_filters[1:n,pos]
    cont = 1
    for i in n:-1:1
        if i % 2 != 0
            ans[cont,2] =(-1) * (value_filters[i,pos])
        else
            ans[cont,2] = (value_filters[i,pos])
        end
        cont +=1
    end
    return ans
end

## Loading csv with values
value_filters = CSV.read("filter_coeff.csv", types=[Float64,Float64,Float64,Float64,Float64,Float64,Float64,Float64,Float64,Float64,Float64,Float64,Float64,Float64,Float64,Float64,Float64,Float64,Float64,Float64,Float64,Float64,Float64,Float64,Float64,Float64,Float64,Float64,Float64,Float64,Float64,Float64])
###Coefficients filters
haar = wav_coeff(value_filters,2,1)
daub2 = wav_coeff(value_filters,4,9)
daub4 = wav_coeff(value_filters,8,11)
daub5 = wav_coeff(value_filters,10,3)
daub6 = wav_coeff(value_filters,12,13)
daub8 = wav_coeff(value_filters,16,15)
daub10 = wav_coeff(value_filters,20,17)
daub12 = wav_coeff(value_filters,24,19)
daub14 = wav_coeff(value_filters,28,21)
daub16 = wav_coeff(value_filters,32,23)
daub18 = wav_coeff(value_filters,36,25)
daub20 = wav_coeff(value_filters,40,27)

daub3 = wav_coeff(value_filters,6,29)
coiff3 = wav_coeff(value_filters,18,31)
sym5 = wav_coeff(value_filters,10,5)
coiff5 = wav_coeff(value_filters,30,7)


#daub4_2_2 = daub4_2[end:-1:1,:]
## Loading csv with the signals
signals = CSV.read("signals.csv"; types=[Float64,Float64,Float64,Float64,Float64,Float64,Float64,Float64,Float64,Float64])
###signals
blocks = zeros(1024,1)
blocks[1:1024,1] = signals[1:1024,1]
bumps = zeros(1024,1)
bumps[1:1024,1] = signals[1:1024,2]
doppler = zeros(1024,1)
doppler[1:1024,1] = signals[1:1024,3]
heaviSine= zeros(1024,1)
heaviSine[1:1024,1] = signals[1:1024,4]
cuarentapi = zeros(1024,1)
cuarentapi = signals[1:1024,9]
preg5 = zeros(1024,1)
preg5 = signals[1:1024,10]
## Loading the signals with additive white gaussian noise
blocks_w = zeros(1024,1)
blocks_w[1:1024,1] = signals[1:1024,5]
bumps_w = zeros(1024,1)
bumps_w[1:1024,1] = signals[1:1024,6]
doppler_w = zeros(1024,1)
doppler_w[1:1024,1] = signals[1:1024,7]
heaviSine_w= zeros(1024,1)
heaviSine_w[1:1024,1] = signals[1:1024,8]

## Load thresholds
thresh = CSV.read("thresholds.csv")

blocks_universal = thresh[1,1]
blocks_rigsure = thresh[2,1]
blocks_minimax = thresh[3,1]

bump_universal = thresh[4,1]
bump_sure = thresh[5,1]
bump_minimax = thresh[6,1]

heavy_universal = thresh[7,1]
heavy_sure = thresh[8,1]
heavy_minimax = thresh[9,1]

doppler_universal = thresh[10,1]
doppler_rigsure = thresh[11,1]
doppler_minimax = thresh[12,1]


## Thresholds sure
thresh_sure = CSV.read("Thresholds_sure.csv")

blocks_coiff = thresh_sure[1:5,1]
blocks_daub = thresh_sure[1:5,2]
blocks_haar = thresh_sure[1:5,3]
blocks_sym = thresh_sure[1:5,4]

bump_coiff = thresh_sure[1:5,5]
bump_daub = thresh_sure[1:5,6]
bump_haar = thresh_sure[1:5,7]
bump_sym = thresh_sure[1:5,8]

heavy_coiff = thresh_sure[1:5,13]
heavy_daub = thresh_sure[1:5,14]
heavy_haar = thresh_sure[1:5,15]
heavy_sym = thresh_sure[1:5,16]

doppler_coiff = thresh_sure[1:5,9]
doppler_daub = thresh_sure[1:5,10]
doppler_haar = thresh_sure[1:5,11]
doppler_sym = thresh_sure[1:5,12]
