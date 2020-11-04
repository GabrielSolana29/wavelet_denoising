using Plots
function samplingData(size)
     x = 0
     y = zeros(size)
     for i =1:size
         x = i
         y[i] = 2*(x^2) - 3*x + 1
     end
     display(plot(1:size,y[1:size]))
     #gui(plot(1:size,y[1:size]))
     return y
end
function samplingData2(size)
     x = 0
     y = zeros(size)
     for i =0:.001:size
         x = i
         y[i] = (20*x^4)*((1-x)^6) * cos(48*pi*x)
     end
     display(plot(1:size,y[1:size]))
     #gui(plot(1:size,y[1:size]))
     return y
end
function noLevels_possible(noLevels_prov,x_coeff)
    s = 0
    for i = 1:noLevels_prov
        s = 2^i
        if s > x_coeff
            return Int32(noLevels_prov - i)
        end
    end
    return Int32(noLevels_prov)
end

### Adding white gaussian noise to the signal x
function awgn(X,SNR)
    #Assumes X to be a matrix and SNR a signal-to-noise ratio specified in decibel (dB)
    #Implented by author, inspired by https://www.gaussianwaves.com/2015/06/how-to-generate-awgn-noise-in-matlaboctave-without-using-in-built-awgn-function/
    N=length(X) #Number of elements in X
    signalPower = sum(X[:].^2)/N
    linearSNR = 10^(SNR/10)
    a,b=size(X)
    noiseMat = randn((a,b)).*√(signalPower/linearSNR) #Random gaussian noise, scaled according to the signal power of the entire matrix (!) and specified SNR

    return solution = X + noiseMat
end

function optimal_level_denoising(wavelet,original_signal,noise_signal,no_samples,threshold,no_levels)
    a_n,d_n = wavelet_analysis(no_samples,noise_signal[:,1],wavelet)
    mse_vec = zeros(no_levels,1)
    for i = 1:no_levels
        clean_signal = denoising(a_n,d_n,threshold,i,no_samples,noise_signal,no_levels,wavelet)
        mse_vec[i] = meanSquareError(original_signal[:,1],clean_signal[:,1])
    end
    position = argmin(mse_vec)
    value = mse_vec[position[1],1]
    return value, position[1]
end

function optimal_level_denoising_sure(wavelet,original_signal,noise_signal,no_samples,threshold,no_levels)
    a_n,d_n = wavelet_analysis(no_samples,noise_signal[:,1],wavelet)
    mse_vec = zeros(no_levels,1)
    for i = 1:no_levels
        clean_signal = denoising_sure(a_n,d_n,threshold,i,no_samples,noise_signal,no_levels,wavelet)
        mse_vec[i] = meanSquareError(original_signal[:,1],clean_signal[:,1])
    end
    position = argmin(mse_vec)
    value = mse_vec[position[1],1]
    return value, position[1]
end

function denoising_sure(ap,de,threshold,level,no_samples,signal,noLevels,wavelet)
    a_a = deepcopy(ap)
    d_d = deepcopy(de)
    for j = 1:level
        size_level = Int(no_samples / 2^(j))
        for i = 1:size_level
                if sqrt(de[i,j+1]^2) <= threshold[j]
                    d_d[i,j+1] = 0
                end
        end
    end
    ans = @time waveletSynthesis(no_samples,a_a[:,level+1],d_d,level,wavelet)
    return ans
end


### Obtain the energy of the original signal
function energy_magnitude_aproximation_original(sampled_signal)
    x = size(sampled_signal)
    sum_a = 0
    for i = 1:x[1]
        sum_a = sum_a + (sampled_signal[i])^2
    end
    return sum_a
end

function energy_magnitude_details(a,d,currentLevel,noLevels)
    x,y = size(d)
    sum = 0
    #for l = 2:currentLevel+1
        for i = 1:x
            sum = sum + (d[i,currentLevel+1])^2
        end
    #end
    return sum
end

function energy_magnitude_aproximation(a,d,currentLevel,noLevels)
    x,y = size(a)
    sum = 0
    for i = 1:x
        sum = sum + (a[i,currentLevel+1])^2
    end
    return sum
end

function energy_levels(a,d,noLevels,complete)
    energy_vect = zeros(noLevels,2)
    #complete = energy_magnitude_aproximation(a,d,0,noLevels)
    for i =1:noLevels
        energy_vect[i,1] = (energy_magnitude_aproximation(a,d,i,noLevels) * 100)/complete
        energy_vect[i,2] = (energy_magnitude_details(a,d,i,noLevels) * 100)/complete
    end
    display(plot(1:noLevels,energy_vect[:,1],title = "Aproximation energy",xaxis = "Level of decomposition", yaxis = "Energy Percentage",label="aproximation",legend =:bottomleft)) #legend=false
    display(plot!(1:noLevels,energy_vect[:,2],title = "Details energy",xaxis = "Level of decomposition", yaxis = "Energy Percentage",label="details"))
    return energy_vect
end
### Receives the original signal and the denoised signal as vectors
function meanSquareError(s1,s2)
    x =  size(s1)
    ans = 0
    for i =1:x[1]
        ans = ans + ((s1[i,1]-s2[i,1])*(s1[i,1]-s2[i,1]))
    end
    ans_f = ans/x[1]
    return ans_f
end

function rootMeanSquareError(s1,s2)
    x =  size(s1)
    ans = 0
    for i =1:x[1]
        ans = ans + (s1[i,1]-s2[i,1])^2
    end
    ans = sqrt(ans/x[1])
    return ans
end

function denoising(ap,de,threshold,level,no_samples,signal,noLevels,wavelet)
    a_a = deepcopy(ap)
    d_d = deepcopy(de)
    for j = 1:level
        size_level = Int(no_samples / 2^(j))
        for i = 1:size_level
                if sqrt(de[i,j+1]^2) <= threshold
                    d_d[i,j+1] = 0
                end
        end
    end
    ans = @time waveletSynthesis(no_samples,a_a[:,level+1],d_d,level,wavelet)
    return ans
end

function sup_val(level,support)
    x = sqrt(1/support)
    return x
end
function reconstruct_aproximations(comp,no_samples,level)
    a_recon = zeros(no_samples)
    a_1 = zeros(no_samples)
    support =  Int(2^level)
    ## Obtain the value for the support if the approximated signals is reconstructed with the coeff
    val_support= sup_val(level,support)
    sp = zeros(support)
    sp .= val_support
    new = sp * comp[1]
    for j = 2:Int(no_samples/(2^level))
        hlp = sp * comp[j]
        new = vcat(new,hlp)
    end
    return new
end

## Analysis wavelet
function approximation_coeff_general(no_samples,currentLevel,coeff,a)
    size_coeff = Int(no_samples / 2^(currentLevel))
    size_prev_level = Int(no_samples / 2^(currentLevel-1))
    ### Create the matrix filter
    x,y = size(coeff)
    filt_mat = zeros(size_coeff,size_prev_level)

    cont = 0
    for i = 1:size_coeff
        if i + cont + x - 1 <= size_prev_level
            filt_mat[i,(i+cont):(i+cont+x-1)] = coeff[:,1]'
        else
            ## Wrapp around
            #rem = size_prev_level-(i+cont)
            rem = size_prev_level-(i+cont-1)
            #filt_mat[i,(i+cont+1):size_prev_level] = coeff[1:rem,1]'
            filt_mat[i,(i+cont):size_prev_level] = coeff[1:rem,1]'
            filt_mat[i,1:x-rem] = coeff[rem+1:x,1]'
        end
        cont +=1
    end

    ## Extract the value of the aproximation coefficients from previous level
    approximation = a[1:size_prev_level,currentLevel]
    filter_coeff = filt_mat * approximation
    ## Concatenate the matrix to be the same size as the array a
    filter_mat_complete = vcat(filter_coeff,zeros(no_samples-size_coeff))
    return filter_mat_complete
end


function wavelet_coeff_general(no_samples,level,coeff,a)
    size_coeff = Int(no_samples / 2^(level))
    size_prev_level = Int(no_samples / 2^(level-1))
    ### Create the matrix filter
    x,y = size(coeff)
    filt_mat = zeros(size_coeff,size_prev_level)
    cont = 0
    for i = 1:size_coeff
        if i + cont + x - 1 <= size_prev_level
            filt_mat[i,(i+cont):(i+cont+x-1)] = coeff[:,2]'
        else
            ## Wrapp around
            #rem = size_prev_level-(i+cont)
            rem = size_prev_level-(i+cont-1)
            #filt_mat[i,(i+cont+1):size_prev_level] = coeff[1:rem,1]'
            filt_mat[i,(i+cont):size_prev_level] = coeff[1:rem,2]'
            filt_mat[i,1:x-rem] = coeff[rem+1:x,2]'
        end
        cont +=1
    end
    ## Extract the value of the aproximation coefficients from previous level
    approximation = a[1:size_prev_level,level]
    ## Multiply filter with previous approximation
    filter_coeff = filt_mat*approximation
    ## Concatenate the matrix to be the same size as the array a
    filter_mat_complete = vcat(filter_coeff,zeros(no_samples-size_coeff))
    return filter_mat_complete
end

function wavelet_general_Analysis(no_samples,sampled_signal,coeff)
    x_coeff,y_coeff=size(coeff)
    noLevels_prov = Int32(log(no_samples)/log(2))
    ## Zeros matrices with approximation coefficients and wavelets coefficients
    a = zeros(no_samples,noLevels_prov+1)
    d = zeros(no_samples,noLevels_prov+1)
    # First column of the matrix is the sampled signal or level 0 of decomposition (julia starts the indexes in 1 instead of 0)
    a[:,1] = sampled_signal'
    ## Calculate the number of levels depending on the size of the coefficients from the wavelet
    noLevels = noLevels_possible(noLevels_prov+1,x_coeff)
    #obtain the coefficients at each level
    for i = 1:noLevels+1
        if i != noLevels+1
            ## Calculate the coefficients
            a[:,i+1] = approximation_coeff_general(no_samples,i,coeff,a)
            d[:,i+1] = wavelet_coeff_general(no_samples,i,coeff,a)
        end
    end
    return a,d
end

## Synthesis wavelets
function synthesis_approximation_coeff_general(no_samples,level,coeff,a2,d)
    size_coeff = Int(no_samples / 2^(level))
    size_prev_level = Int(no_samples / 2^(level-1))

    ### Create the matrix filter aproximation
    x,y = size(coeff)
    filt_mat = zeros(size_coeff,size_prev_level)

    cont = 0
    for i = 1:size_coeff
        if i + cont + x - 1 <= size_prev_level
            filt_mat[i,(i+cont):(i+cont+x-1)] = coeff[:,1]'
        else
            ## Wrapp around
            #rem = size_prev_level-(i+cont)
            rem = size_prev_level-(i+cont-1)
            #filt_mat[i,(i+cont+1):size_prev_level] = coeff[1:rem,1]'
            filt_mat[i,(i+cont):size_prev_level] = coeff[1:rem,1]'
            filt_mat[i,1:x-rem] = coeff[rem+1:x,1]'
        end
        cont +=1
    end
    filt_mat_approx = filt_mat
    ### Create the coeff filter wavelet
    x,y = size(coeff)
    filt_mat = zeros(size_coeff,size_prev_level)

    cont = 0
    for i = 1:size_coeff
        if i + cont + x - 1 <= size_prev_level
            filt_mat[i,(i+cont):(i+cont+x-1)] = coeff[:,2]'
        else
            ## Wrapp around
            #rem = size_prev_level-(i+cont)
            rem = size_prev_level-(i+cont-1)
            #filt_mat[i,(i+cont+1):size_prev_level] = coeff[1:rem,1]'
            filt_mat[i,(i+cont):size_prev_level] = coeff[1:rem,2]'
            filt_mat[i,1:x-rem] = coeff[rem+1:x,2]'
        end
        cont +=1
    end

    filt_mat_detail = filt_mat
    #if level == 6
    #    CSV.write("testFiter.csv",  DataFrame(filt_mat_detail), writeheader=false)
    #end
    ## Extract the value of the aproximation coefficients from next level
    approximation = a2[1:size_coeff,1]
    ## Multiply filter with previous approximation (the matrix is transposed)
    filter_coeff = (filt_mat_approx' * approximation) + (filt_mat_detail' *  d[1:size_coeff,1])
    #return filter_coeff,approximation,d
    x = size(filter_coeff)
    ## Concatenate the matrix to be the same size as the array a
    filter_mat_complete = vcat(filter_coeff,zeros(no_samples-x[1]))

    return filter_mat_complete
end



function waveletSynthesis_general(noLevels,no_samples,a,d,currentLevel,coeff)
    ## Zeros matrices with approximation coefficients and wavelets coefficients
    a2 = zeros(no_samples,noLevels+1)
    a2[:,currentLevel+1] = a[:,1]
    x_coeff,y_coeff=size(coeff)
    noLevels_prov = Int32(log(no_samples)/log(2))
    ## Calculate the number of levels depending on the size of the coefficients from the wavelet
    noLevels = noLevels_possible(noLevels_prov+1,x_coeff)

    ## Obtain coefficients for each level
    for i in currentLevel+1:-1:1
        if i <= noLevels
            if i != currentLevel+1
                ## Calculate the coefficients
                #### ESTA FUNCIÓN FUNCIONA SOLO CUANDO NO ES MAYOR QUE EL NUMERO DE NIVELES POSIBLES
                a2[:,i] =  synthesis_approximation_coeff_general(no_samples,i,coeff,a2[:,i+1],d[:,i+1])
            end
        else
            ### PROGRAMAR UNA FUNCÍON CON LA MODIFICACIÓN
        end
    end
    return a2
end

#######################################################################################################
#################################### JUST FOR TEST THE VECTOR WRAPP AROUND ############################

function wavelet_analysis(no_samples,sampled_signal,coeff)
    x_coeff,y_coeff=size(coeff)
    noLevels_prov = Int32(log(no_samples)/log(2))
    ## Zeros matrices with approximation coefficients and wavelets coefficients
    a = zeros(no_samples,noLevels_prov+1)
    d = zeros(no_samples,noLevels_prov+1)
    # First column of the matrix is the sampled signal or level 0 of decomposition (julia starts the indexes in 1 instead of 0)
    a[:,1] = sampled_signal'
    ## Calculate the number of levels depending on the size of the coefficients from the wavelet
    noLevels = noLevels_possible(noLevels_prov+1,x_coeff)
    #obtain the coefficients at each level
    for i = 1:noLevels_prov
        if i < noLevels+1
            if i != noLevels+1
                ## Calculate the coefficients
                a[:,i+1] = approximation_coeff_general(no_samples,i,coeff,a)
                d[:,i+1] = wavelet_coeff_general(no_samples,i,coeff,a)
            end
        elseif i == noLevels_prov # Last level you are going to write the last coef
            a[:,i+1] = scaling_coeff_last(no_samples,i,coeff,a)
            d[:,i+1] = wavelet_coeff_last(no_samples,i,coeff,a)
        else
            #a[:,i+1] = fastImple_last(no_samples,i,coeff,a,noLevels)
            a[:,i+1] = scaling_coeff_lasts(no_samples,i,coeff,a)
            d[:,i+1] = wavelet_coeff_lasts(no_samples,i,coeff,a)
        end
    end
    return a,d
end

function scaling_coeff_lasts(no_samples,level,coeff,a)
    size_coeff = Int(no_samples / 2^(level))
    size_prev_level = Int(no_samples / 2^(level-1))
    x,y = size(coeff)
    s = zeros(size_coeff,1)
    cont2 = 0
    for c = 1:size_coeff
        suma = 0
        cont = cont2 + 1
        for i = 1:x
            if cont > size_prev_level
                cont = 1
            end
            suma = suma + (coeff[i,1] * a[cont,level])
            cont = cont + 1
        end
        cont2 = cont2 + 2
        s[c] = suma
    end
    s = vcat(s[:,1],zeros(no_samples-size_coeff))
    return s
end

# Ya esta bien
function scaling_coeff_last(no_samples,level,coeff,a)
    size_coeff = Int(no_samples / 2^(level))
    size_prev_level = Int(no_samples / 2^(level-1))
    x,y = size(coeff)
    s = zeros(1)
    cont = 1
    for i = 1:x
        if cont > size_prev_level
            cont = 1
        end
        s[1] = s[1] + (coeff[i,1] * a[cont,level])
        cont = cont + 1
    end
    s = vcat(s,zeros(no_samples-1))
    return s
end

function wavelet_coeff_lasts(no_samples,level,coeff,a)
    size_coeff = Int(no_samples / 2^(level))
    size_prev_level = Int(no_samples / 2^(level-1))
    x,y = size(coeff)
    s = zeros(size_coeff,1)
    cont2 = 0
    for c = 1:size_coeff
        suma = 0
        cont = cont2 + 1
        for i = 1:x
            if cont > size_prev_level
                cont = 1
            end
            suma = suma + (coeff[i,2] * a[cont,level])
            cont = cont + 1
        end
        cont2 = cont2 + 2
        s[c] = suma
    end
    s = vcat(s[:,1],zeros(no_samples-size_coeff))
    return s
end

# Ya esta bien
function wavelet_coeff_last(no_samples,level,coeff,a)
    size_coeff = Int(no_samples / 2^(level))
    size_prev_level = Int(no_samples / 2^(level-1))
    x,y = size(coeff)
    s = zeros(1)
    cont = 1
    for i = 1:x
        if cont > size_prev_level
            cont = 1
        end
        s[1] = s[1] + (coeff[i,2] * a[cont,level])
        cont = cont + 1
    end
    s = vcat(s,zeros(no_samples-1))
    return s
end

### Synthesis process
function waveletSynthesis(no_samples,a,d,currentLevel,coeff)
    x_coeff,y_coeff=size(coeff)
    ## Calculate the number of levels depending on the size of the coefficients from the wavelet
    noLevels_prov = Int32(log(no_samples)/log(2))
    noLevels = noLevels_possible(noLevels_prov+1,x_coeff)
    ## Zeros matrices with approximation coefficients and wavelets coefficients
    a2 = zeros(no_samples,noLevels_prov+1)
    a2[:,currentLevel+1] = a[:,1]
    ## Obtain coefficients for each level
    for i in currentLevel+1:-1:1
        if i < noLevels+1
            if i != currentLevel+1
                ## Calculate the coefficients
                #### ESTA FUNCIÓN FUNCIONA SOLO CUANDO NO ES MAYOR QUE EL NUMERO DE NIVELES POSIBLES
                a2[:,i] = synthesis_approximation_coeff_general(no_samples,i,coeff,a2[:,i+1],d[:,i+1])
            end
        elseif i == currentLevel+1
            #a2[:,i] = synthesis_filt_test_ultimo_nivel(no_samples,i,coeff,a2[:,i],d[:,i+1])
            a2[:,i-1] = synthesis_filt(no_samples,i-1,coeff,a2[:,i],d[:,i])
        else
            a2[:,i-1] = synthesis_filt(no_samples,i-1,coeff,a2[:,i],d[:,i])
        end
    end
    return a2
end

function synthesis_filt(no_samples,level,coeff,a,d)
    #Con esto obtengo las matrices y necesito transponerlas para que funcione y luego se sume
    if no_samples / 2^(level) == .5
        size_coeff = 1
    else
        size_coeff = Int(no_samples / 2^(level))
    end
    size_prev_level = Int(no_samples / 2^(level-1))
    x,y = size(coeff)
    s = zeros(size_prev_level,1)
    s_prov = zeros(x,size_coeff)

    for i = 1:size_coeff
        s_prov[:,i] = (a[i,1]*coeff[:,1]) + (d[i,1]*coeff[:,2])
    end
    j2 = 0
    for k = 1:size_prev_level
        cont = 0
        for i = 1:size_coeff
                if cont == 0
                    j2 = 1
                else
                    j2 = cont+1
                end
                for j = 1:x
                    if k == j2
                        s[k,1] += s_prov[j,i]
                    end
                    if j2 >= size_prev_level
                        j2 = 1
                    else
                        j2 +=1
                    end
                end
                cont = cont + 2
            end
    end
    s = vcat(s[:,1],zeros(no_samples-size_prev_level))
    return s
end



## Compression

### Function for performing compression
function compressionAprox_coef(a,d,level,no_samples,energy_threshold)
    ## Obtain the enrgy of the original signal
    ax,ay = size(a)
    maxEnergy = 0
    for i = 1:ax
         #global maxEnergy = maxEnergy + (a[i,1])^2
         maxEnergy = maxEnergy + (a[i,1])^2
    end
    # Find the energy of each coefficient
    no_samples_level = Int(no_samples / (2^level))
    newVec = zeros(no_samples_level)
    for i = 1:no_samples_level
        #global level
        #global newVec[i] = (a[i,level+1])^2
         newVec[i] = (a[i,level+1])^2
    end
    ## this method returns the order by the smallest to the largest of the vector
    index_vec = sortperm(newVec)
    ## Eliminate half of the coefficients
    #respaldo = copy(newVec)
    half_coeff =  copy(newVec)
    energy_threshold_vec = copy(newVec)

    ## Find the compression by eleiminating half of the smallest coefficients
    for i = 1:no_samples_level
        if index_vec[i] <= no_samples_level/2
            half_coeff[i] = 0
        end
    end
    ## find compression by setting a threshold
    sum_energy = 0
    cont = no_samples_level
    cont_2= 0
    sum_energy_percentage = 0
    while sum_energy_percentage < energy_threshold

        if cont == 0
            print("\n Threshold couldn't be reached!!")
            return 0
        end
        ## To avoid infinite looping
        ## iterates in reverse because it goes from the larges to the smallest coefficient
        sum_energy = sum_energy + energy_threshold_vec[index_vec[cont]]
        sum_energy_percentage = ((sum_energy * 100) / maxEnergy)
        cont_2 = cont_2 + 1

        cont = cont - 1
    end

    for i = 1:no_samples_level
        if i > cont_2
            energy_threshold_vec[index_vec[no_samples_level + 1 - i]] = 0
        end
    end

    energy_threshold_percentage = (sum(energy_threshold_vec)*100)/maxEnergy
    energy_threshold_vec = sqrt.(energy_threshold_vec)
    half_coeff_percentage = (sum(half_coeff) * 100) / maxEnergy
    half_coeff = sqrt.(half_coeff)
    print("\n The number of samples required to reach the threshold are: ", cont_2)
    print("\n The energy percentage obtained by removing half the coefficients is: ", half_coeff_percentage)
    print("\n The energy percentage obtained by using the energy threshold is: ", sum_energy_percentage)

    return half_coeff_percentage,half_coeff,energy_threshold_percentage,energy_threshold_vec, cont_2
end
