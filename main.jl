#### Include all the .jl files with the necessary implemented functions
include("csv_loading.jl")
include("wavelet_general.jl")
# Adding white noise to signals
#=
blocks_wg = awgn(blocks,3)
bumps_wg = awgn(bumps,3)
doppler_wg = awgn(doppler,3)
heavySine_wg = awgn(heaviSine,3)
=#
#######################################################################################################
#################################### JUST FOR TEST THE VECTOR WRAPP AROUND ############################
no_samples = 1024
noLevels = 10
sampled_signal = heaviSine_w
#plot(1:1024,doppler,label="",title="Doppler function ",xaxis = "samples")
wavelet = daub5
a,d = wavelet_analysis(no_samples,sampled_signal[:,1],wavelet)
currentLevel = 10 ##level  from where reconstruction starts
a2= @time waveletSynthesis(no_samples,a[:,currentLevel+1],d,currentLevel,wavelet)
### level to which the denoising algorithm will iteratet to
CSV.write("approximations.csv",  DataFrame(a), writeheader=false)
CSV.write("details.csv",  DataFrame(d), writeheader=false)
CSV.write("reconstruction.csv",  DataFrame(a2), writeheader=false)

### Ortogonalyty check

no_samples = 1024
noLevels = 10
sampled_signal = doppler
wavelet = daub5

anst = daub5[:,1] .* daub5[:,2]
anst = sum(anst)
plot(1:10,daub5[:,1])

## Returns the mse and the optimal level for denoising
no_samples = 1024
noLevels = 5
threshold = heavy_sym
signal_noisy = heaviSine_w
original_signal = heaviSine
wavelet = sym5
## Find the mse value and the optimal level to perform denoisisng
value,level = @time optimal_level_denoising(wavelet,original_signal,signal_noisy,no_samples,threshold,noLevels)
#value,level = @time optimal_level_denoising_sure(wavelet,original_signal,signal_noisy,no_samples,threshold,noLevels)
### Obtain the figure of the denoised signal
a,d = wavelet_analysis(no_samples,signal_noisy[:,1],wavelet)
#clean_signal = denoising(a,d,threshold,level,no_samples,signal_noisy,noLevels,wavelet)
clean_signal = denoising_sure(a,d,threshold,level,no_samples,signal_noisy,noLevels,wavelet)
### Plot the signal
p0 =plot(1:no_samples,original_signal,title="Original Signal",label="",xlabel = "samples", titlefont="black",fontfamily=font(42, "Arial"))
p1 = plot(1:no_samples,signal_noisy,title="Signal with Additive Gaussian Noise",label="",xlabel = "samples")
p2 = plot(1:no_samples,clean_signal[:,1],title = "Denoised Signal",label="",xlabel = "samples")
plot_f=plot(p0,p1,p2,layout = (3,1))
### Save figure as png
mse_string = string(value)
level_string = string(level)
wavelet_name = "sym_"
filter_string = "sure_"
signal_name = "heavy_"
name = wavelet_name * signal_name *filter_string* level_string* "_" * mse_string * ".png"
savefig(plot_f,name)
#CSV.write("Clean_signal.csv",  DataFrame(aClean), writeheader=false)

#=
## Explanation of the proposed algorithm for reconstruction
### DAUB2
z = (a_test[1,11]*daub2[:,1]) + (d_test[1,11]*daub2[:,2])
z1 = z[1] + z[3]
z2 = z[2] + z[4]

z_1 = (a_test[1,10]*daub2[:,1]) + (d_test[1,10]*daub2[:,2])
z_2 = (a_test[2,10]*daub2[:,1]) + (d_test[2,10]*daub2[:,2])

z1 = z_1[1] + z_2[3]
z2 = z_1[2] + z_2[4]
z3 = z_1[3] + z_2[1]
z4 = z_1[4] + z_2[2]

z_1 = (a_test[1,9]*daub2[:,1]) + (d_test[1,9]*daub2[:,2])
z_2 = (a_test[2,9]*daub2[:,1]) + (d_test[2,9]*daub2[:,2])
z_3 = (a_test[3,9]*daub2[:,1]) + (d_test[3,9]*daub2[:,2])
z_4 =(a_test[4,9]*daub2[:,1]) + (d_test[4,9]*daub2[:,2])

# Ultimo nivel DAUB4
z = (a_test[1,11]*daub4[:,1]) + (d_test[1,11]*daub4[:,2])
z1 = z[1]+z[3]+z[5]+z[7]
z2 = z[2]+z[4]+z[6]+z[8]

# Penultimo nivel
z_1 = (a_test[1,10]*daub4[:,1]) + (d_test[1,10] * daub4[:,2])
z_2 = (a_test[2,10]*daub4[:,1]) + ( d_test[2,10] * daub4[:,2])

z1 = z_1[1] + z_1[5] + z_2[3] + z_2[7]
z2 = z_1[2] + z_1[6] + z_2[4] + z_2[8]
z3 = z_1[3] + z_1[7] + z_2[1] + z_2[5]
z4 = z_1[4] + z_1[8] + z_2[6] + z_2[2]


z_1 = (a_test[1,9]*daub4[:,1]) + (d_test[1,9]*daub4[:,2])
z_2 = (a_test[2,9]*daub4[:,1]) + (d_test[2,9]*daub4[:,2])
z_3 = (a_test[3,9]*daub4[:,1]) + (d_test[3,9]*daub4[:,2])
z_4 =(a_test[4,9]*daub4[:,1]) + (d_test[4,9]*daub4[:,2])

z1 = z_1[1] + z_2[7] +z_3[5] + z_4[3]
z2 = z_1[2] + z_2[8] +z_3[6] + z_4[4]
z3 = z_1[3] + z_2[1] +z_3[7] + z_4[5]
z4 = z_1[4] + z_2[2] +z_3[8] + z_4[6]
z5 = z_1[4] + z_2[2] +z_3[8] + z_4[6]
z6 = z_1[4] + z_2[2] +z_3[8] + z_4[6]
z7 = z_1[4] + z_2[2] +z_3[8] + z_4[6]
z8 = z_1[4] + z_2[2] +z_3[8] + z_4[6]
# Ultimo nivel daub5
z = (a_test[1,11]*daub5[:,1]) + (d_test[1,11]*daub5[:,2])
ho = z[1]+z[3]+z[5]+z[7]+z[9]
hu = z[2]+z[4]+z[6]+z[8]+z[10]

# Penultimo nivel
z1 = (a[1,10]*daub5[:,1]) + (d[1,10] * daub5[:,2])
z2 = (a[2,10]*daub5[:,1]) + (d[2,10] * daub5[:,2])

z_ans = zeros(4)
z_ans[1] = z1[1] + z1[5] + z1[9] + z2[3] + z2[7]
z_ans[2] = z1[2] + z2[4] + z1[6] + z2[8]+ z1[10]
z_ans[3] = z1[3] + z2[5] + z1[7] + z2[9]+ z2[1]
z_ans[4] = z1[4] + z2[6] + z1[8] + z2[10]+ z2[2]

#antepenultimo nivel
z1 = (a_test[1,9]*daub5[:,1]) + (d_test[1,9] * daub5[:,2])
z2 = (a_test[2,9]*daub5[:,1]) + ( d_test[2,9] * daub5[:,2])
z3 = (a_test[3,9]*daub5[:,1]) + ( d_test[3,9] * daub5[:,2])
z4 = (a_test[4,9]*daub5[:,1]) + ( d_test[4,9] * daub5[:,2])

zans = zeros(8)
zans[1] = z1[1]+z1[9]+z2[7]+z3[5]+z4[3]
zans[2] = z1[2]+z1[10]+z2[8]+z3[6]+z4[4]
zans[3] = z1[3]+z2[1]+z2[9]+z3[7]+z4[5]
=#
