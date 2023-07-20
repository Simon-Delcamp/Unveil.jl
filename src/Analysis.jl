module Analysis

include("Dataprep.jl") #Read and write fits
using .Dataprep

import MultivariateStats, StatsBase, Statistics
using  FFTW, AbstractFFTs

export fourmoments
export metricOW
export metricPCA
export rmscube


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

            hist  = StatsBase.fit(StatsBase.Histogram,cube[:,dx],nbins=bin)
            histnall = StatsBase.normalize(hist,mode=:pdf)
            histn = histnall.weights
            edges = Dataprep.flatiterator(histnall.edges)[1,:]
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
        hist  = StatsBase.fit(StatsBase.Histogram,cube,nbins=bin)
        histnall = StatsBase.normalize(hist,mode=:pdf)
        histn = histnall.weights
        edges = Dataprep.flatiterator(histnall.edges)[1,:]
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
    metricOW(mom1,mom2,mom3,mom4,dv,SIGMAT)


Calculate the following metric on each values of mom1,2,3 and 4: sqrt((mom1)**2+((mom2-SIGMAT))**2+(mom3)**2+(mom4-3)**2)
"""
function metricOW(mom1,mom2,mom3,mom4,dv,SIGMAT)
    metric = similar(mom1)
    for ix=1:size(metric)[1]
        metric[ix] = sqrt((mom1[ix])^2 + ((mom2[ix]-SIGMAT))^2+(mom3[ix])^2+(mom4[ix]-3)^2)
    end
    return(metric)

end#metricOW


"""
    metricPCA(mom1,mom2,mom3,mom4,dv)


Calculate the following metric on each values of mom1,2,3 and 4: sqrt((mom1/dv)**2+(mom2/dv)**2+(mom3)**2+(mom4-3)**2)
"""
function metricPCA(mom1,mom2,mom3,mom4,dv)
    metric = similar(mom1)
    for ix=1:size(metric)[1]
        metric[ix] = sqrt((mom1[ix]/dv)^2 + (mom2[ix]/dv)^2+(mom3[ix])^2+(mom4[ix]-3)^2)
    end
    return(metric)

end#metric





"""
    power_spectra(arr,imsize ; fitted=true)

Compute the power spectrum of arr, and return it with the Fourier domain vector.
"""
function power_spectra(arr,imsize ;BLANK=-1000)
    Np   = trunc(Int,imsize/2)
    karr =  rfftfreq(imsize)

    arr = Dataprep.replace_nantomissing(arr)   
    arr  = Dataprep.replace_blanktomissing(arr,BLANK)
    arr  = Dataprep.replace_missingtoblank(arr,0)
    arr = convert(Matrix{Float64},arr)
    arr  = abs.(fft(arr)).^2

    sumi = Array{Float64}(undef,size(arr)[2])
    mea  = similar(sumi)
    many = similar(sumi)
    many .= 0
    for px = 1:size(arr)[2]
        sumi[px] = sum(skipmissing(arr[:,px]))
        for py=1:size(arr)[1]
            if ismissing(arr[py,px])==false
                many[px] += 1
            end 
        end
    end
    mea .= sumi./many
    return(mea, karr)


   # return(coeff)
end




"""
    rms_cube(cube,can)

Compute the rms on velocity canal given as input of a cube. Return a 2D map with given rms on each pixel and the averaged rms accross the map.
"""
function rms_cube(cube,can)
    map = Array{Float64}(undef,size(cube)[1],size(cube)[2])
    for ix=1:size(cube)[2]

        for jx=1:size(cube)[1]
            map[jx,ix] = StatsBase.moment((cube[jx,ix,can]),2)
        end
    end
    rmsavr = StatsBase.mean(skipmissing(map))
    return(map,rmsavr)


end






end #module