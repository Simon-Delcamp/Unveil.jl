###################################################################
# Fonctions used for analyses purpose. Computing rms of maps, snr, produce gaussian.
# Called these functions by calling (or writting on a script ):
#       >include("../src/Data_analysis.jl")
#       >using .Data_analysis
#       >output = Data_analysis.NameOfTheFunction(input)
###################################################################


module Data_analysis


include("Functionforpca.jl") #Calculations of PCA
include("Functionforcvi.jl") #Calculations of CVI
include("Data_preparation.jl") #Read and write fits

using .Functionforcvi
using .Functionforpca
using .Data_preparation

using MultivariateStats, Plots, Statistics, Format, Profile, Distributions, StatsBase
using StatsPlots, Measures, Interact

export calcmetric
export calcmetricOW
export find_indices
export find_threshold
export fourmoments
export generate_gaussian
export integration_bysection
export rms_analytic
export rms_analytic_field
export snr
export snr_allfield




"""
    calcmetric(mom1,mom2,mom3,mom4,dv)


Calculate the following metric on each values of mom1,2,3 and 4: sqrt((mom1/dv)**2+(mom2/dv)**2+(mom3)**2+(mom4-3)**2)
"""
function calcmetric(mom1,mom2,mom3,mom4,dv)
    metric = similar(mom1)
    for ix=1:size(metric)[1]
        metric[ix] = sqrt((mom1[ix]/dv)^2 + (mom2[ix]/dv)^2+(mom3[ix])^2+(mom4[ix]-3)^2)
    end
    return(metric)

end#calcmetric



"""
    calcmetricOW(mom1,mom2,mom3,mom4,dv,SIGMAT)


Calculate the following metric on each values of mom1,2,3 and 4: sqrt((mom1)**2+((mom2-SIGMAT))**2+(mom3)**2+(mom4-3)**2)
"""
function calcmetricOW(mom1,mom2,mom3,mom4,dv,SIGMAT)
    metric = similar(mom1)
    for ix=1:size(metric)[1]
        metric[ix] = sqrt((mom1[ix])^2 + ((mom2[ix]-SIGMAT))^2+(mom3[ix])^2+(mom4[ix]-3)^2)
    end
    return(metric)

end#calcmetricOW


"""
    find_indices(specdiff,arra ; dims = size(arra))

Return an array of the same dimension than arra where the indices in specdiff are equal to 1000, and others indices are equal to 0. Dims indicates the dimensions of the returned array (by default same as in input). Use this if some dimensions should not be used.
"""
function find_indices(specdiff,arra ; dims = size(arra))
    specdiff_map = zeros(Float64,dims)
    #specdiff_map = similar(arra,Float64,dims)
    count           = 1 # for iteration on the reconstructed data without missing values
    for ix in eachindex(specdiff_map)
        count>size(specdiff[])[1] && break
        ix==specdiff[][count] && (specdiff_map[ix]=1000.0)
        ix==specdiff[][count] && (count+=1)
    end
    return(specdiff_map)
end



"""
    find_threshold(array1,array2,thresh)

Return an array with each value being the index where abs(array1[indice]-array2[indice]) is greater than a threshold.
"""
function find_threshold(array1,array2,thresh)
    differences    = abs.(array1.-array2)
    differences    = Data_preparation.replace_missingtoblank(differences,0)
    spectratokeep  = []
    push!(spectratokeep,findall(x->x>=thresh,differences))
    return(spectratokeep)
end


""" 
    fourmoments(cube;dim=2)

Compute the first four moments order of a given cube of dimension 'dim'. If 2D PV cube, will compute the moments on each row. Return the four moments. They are Vector of dimension=size(cube)[1]

"""
function fourmoments(cube;dim=2,bin=50)
    if dim==2
        mom1 = Array{Float64}(undef,size(cube)[2])
        mom2 = similar(mom1)
        mom3 = similar(mom1)
        mom4 = similar(mom1)
        for dx=1:size(cube)[2]

            hist  = StatsBase.fit(Histogram,cube[:,dx],nbins=bin)
            histnall = StatsBase.normalize(hist,mode=:pdf)
            histn = histnall.weights
            edges = Data_preparation.flatiterator(histnall.edges)[1,:]
            mudiftemp = Array{Float64}(undef,size(histn)[1])
            sigdiftemp = similar(mudiftemp)
            skediftemp = similar(mudiftemp)
            kurdiftemp = similar(mudiftemp)
            for nx=1:size(histn)[1]
                mudiftemp[nx] = histn[nx]*edges[nx]
            end
            mudif= sum(mudiftemp)
            mudiftempp = sum(histn) 
            mom1[dx] = mudif/mudiftempp

            for nx=1:size(histn)[1]
                sigdiftemp[nx] = histn[nx]*(edges[nx]-mom1[dx])^2
                skediftemp[nx] = histn[nx]*(edges[nx]-mom1[dx])^3
                kurdiftemp[nx] = histn[nx]*(edges[nx]-mom1[dx])^4
            end
            sigdif = sum(sigdiftemp)
            mom2[dx] = sqrt(sigdif/mudiftempp)

            skedif = sum(skediftemp)
            mom3[dx] = skedif/mudiftempp/mom2[dx].^3

            kurdif = sum(kurdiftemp)
            mom4[dx] = kurdif/mudiftempp/mom2[dx].^4
        end #for dx
    elseif dim==1
        #mom1 = Array{Float64}(undef,size(cube)[1])
        #mom2 = similar(mom1)
        #mom3 = similar(mom1)
        #mom4 = similar(mom1)
        hist  = StatsBase.fit(Histogram,cube,nbins=bin)
        histnall = StatsBase.normalize(hist,mode=:pdf)
        histn = histnall.weights
        edges = Data_preparation.flatiterator(histnall.edges)[1,:]
        mudiftemp = Array{Float64}(undef,size(histn)[1])
        sigdiftemp = similar(mudiftemp)
        skediftemp = similar(mudiftemp)
        kurdiftemp = similar(mudiftemp)
        for nx=1:size(histn)[1]
            mudiftemp[nx] = histn[nx]*edges[nx]
        end
        mudif= sum(mudiftemp)
        mudiftempp = sum(histn) 
        mom1 = mudif/mudiftempp

        for nx=1:size(histn)[1]
            sigdiftemp[nx] = histn[nx]*(edges[nx]-mom1)^2
            skediftemp[nx] = histn[nx]*(edges[nx]-mom1)^3
            kurdiftemp[nx] = histn[nx]*(edges[nx]-mom1)^4
        end
        sigdif = sum(sigdiftemp)
        mom2 = sqrt(sigdif/mudiftempp)

        skedif = sum(skediftemp)
        mom3 = skedif/mudiftempp/mom2.^3

        kurdif = sum(kurdiftemp)
        mom4 = kurdif/mudiftempp/mom2.^4
    end #if dim
    return(mom1,mom2,mom3,mom4)
end



"""
    generate_gaussian(xarray,mu,sig)

Produce a gaussian and a normalized gaussian, following mu and sigma.
"""
function generate_gaussian(xarray,mu,sig)
    # nor : normal distribution
    nor = Normal(mu,sig)
    dif = ((xarray[:]).-(mean(nor)))
    gau = exp.(-(xarray.-mu).^2/2.0/sig^2)
    gau_normalized = exp.(-0.5*(dif/std(nor)).^2)/(std(nor)*sqrt(2pi))
    return(gau,gau_normalized)
end



"""
    integration_bysection(data2D,integration_range,datadim,increment)

Compute an integration of spectra with a given increment, for each spectra in a map. Return an array with the first dimension giving a different spectra, and the second dimension giving the integration sections.
"""
function integration_bysection(data2D,integration_range,datadim,velocity_increment)

    numberofsection     = trunc(abs(datadim[2]*velocity_increment)/integration_range) |> Int
    sectionsize_indices = trunc(datadim[2]/numberofsection)|> Int
    integrated          = zeros(Float64,datadim[1],numberofsection)
    for ix=1:numberofsection
        if ((ix+1)*sectionsize_indices)>datadim[2]
                integrated[:,ix] = abs.(sum(data2D[:,ix*sectionsize_indices:end],dims=2).*abs(velocity_increment))
        else
                integrated[:,ix] = abs.(sum(data2D[:,ix*sectionsize_indices:(ix+1)*sectionsize_indices],dims=2).*abs(velocity_increment))
        end
    end
    return(integrated)
end







"""
    rms_analytic(yarr,rms,xarr)

Return the uncertainty of the first moment order in velocity, and some part of the analytical expression. Yarr represents the emission part of the data, xarr the velocity vector.
"""
function rms_analytic(yarr,rms,xarr)
    incr  = abs(xarr[2]-xarr[1])
    nchan = size(xarr)[1]
    # Denominator
    den   = sum(yarr)*incr
    # Numerator
    num1  = nchan*incr^2/den^2
    num2  = sum(xarr.^2)/(sum(xarr.*yarr))^2
    # 1st order moment
    mu    = abs(mean(xarr,aweights(yarr)))
    # Result
    res1  = (rms*mu)^2*num1
    res2  = (rms*mu)^2*num2
    res   = sqrt(res1+res2)
    return(res,res1,res2,den)
end



"""
    rms_analytic_field(array1,xarr,noise_canals)

Return the uncertainty of the first moment order in velocity of a field, and some part of the analytical expression. Xarr represents the velocity vector. Array1 is a 2D data cube, with each row a new pixel and the columns forms the spectra.
"""
function rms_analytic_field(cube,xarr,noise_canals)
    rms = mean((std(cube[:,noise_canals[1]])+std(cube[:,noise_canals[2]]))/2)
    uncertainty = zeros(Float64,size(cube)[1])
    for ix=1:size(cube)[1]
        uncertainty[ix] = rms_analytic(cube[ix,:],rms,xarr)[1]
    end
    return(uncertainty)
end



"""
    snr(yarr,noise_canals)

    Return the SNR of a signal (spectra for example). yarr have to be a vector.
"""
function snr(yarr,noise_canals)
    snr = maximum(yarr[:])/std(yarr[noise_canals])
    return(snr)
end



"""
    snr_allfield(arr,noise_canals)

Return the mean of the SNR calculated on several spectra. Arr can be 2D (PV) or 3D (PPV) but will be faster if 2D. Better to avoid missing values in your data.
"""
function snr_allfield(arr,noise_canals)
    typeof(size(arr))==Tuple{Int,Int,Int} && (arr=reshape(arr,size(arr)[1]*size(arr)[2],size(arr)[3]))
    snr_field = zeros(Float64,size(arr)[1])
    for ix=1:size(arr)[1]
           snr_field[ix]=(Data_analysis.snr(arr[ix,:],noise_canals[1])+Data_analysis.snr(arr[ix,:],noise_canals[2]))/2
    end
    return(mean(snr_field))
end

end
