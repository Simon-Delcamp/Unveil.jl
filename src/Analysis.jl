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
            mom1[dx] = StatsBase.moment(cube[:,dx],1,0)
            mom2[dx] = sqrt(StatsBase.moment(cube[:,dx,],2))


            mom3[dx] = StatsBase.moment(((cube[:,dx].-mom1[dx])./mom2[dx]).^3,1,0)
            mom4[dx] = StatsBase.moment(((cube[:,dx].-mom1[dx])./mom2[dx]).^4,1,0)

        end #for dx
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
    metricPCA(mom1,mom2,mom3,mom4)


Calculate the following metric on each values of mom1,2,3 and 4: sqrt((mom1/dv)**2+(mom2/dv)**2+(mom3)**2+(mom4-3)**2)
"""
function metricPCA(mom1,mom2,mom3,mom4,PCOPT,dv)
    metric = similar(mom1)

    #norm(i) = ((i-moment(i[PCOPT:end]),1,0)/moment(i[PCOPT:end],2))^2
    a = (mom1.-StatsBase.moment(mom1[PCOPT:end],1,0))./StatsBase.moment(mom1[PCOPT:end],2)
    b = (mom2.-StatsBase.moment(mom2[PCOPT:end],1,0))./StatsBase.moment(mom2[PCOPT:end],2)
    c = (mom3.-StatsBase.moment(mom3[PCOPT:end],1,0))./StatsBase.moment(mom3[PCOPT:end],2)
    d = (mom4.-StatsBase.moment(mom4[PCOPT:end],1,0))./StatsBase.moment(mom4[PCOPT:end],2)
    #metric .= sqrt.(a.^2+b.^2+c.^2+d.^2)
    a = (mom1./0.05)#.-StatsBase.moment(mom1[PCOPT:end],1,0))./StatsBase.moment(mom1[PCOPT:end],2)
    b = (mom2)#.-StatsBase.moment(mom2[PCOPT:end],1,0))./StatsBase.moment(mom2[PCOPT:end],2)
    c = (mom3)#.-StatsBase.moment(mom3[PCOPT:end],1,0))./StatsBase.moment(mom3[PCOPT:end],2)
    d = (mom4.-3)#.-StatsBase.moment(mom4[PCOPT:end],1,0))./StatsBase.moment(mom4[PCOPT:end],2)

    metric .= sqrt.(a.^2+b.^2+c.^2+d.^2)
    #metric .= sqrt(norm.(mom1)+norm.(mom2)+norm.(mom3)+norm.(mom4))
    #for ix=1:size(metric)[1]
    #    metric[ix] = sqrt((mom1[ix]/dv)^2 + (mom2[ix]/dv)^2+(mom3[ix])^2+(mom4[ix]-3)^2)
       # metric[ix] = sqrt((mom1[ix]/dv)^2 + (mom2[ix]/dv)^2+(mom3[ix]/dv^3)^2+((mom4[ix]-3)/dv^4)^2)
       # metric[ix] = sqrt(norm(1)+norm(2)+norm(3)+norm(4))
    #end
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

Compute the dispersion (e.g. sqrt(moment order 2)) on velocity canal given as input of a cube. Return a 2D map with given dispersion on each pixel and the averaged dispersion accross the map.
"""
function rms_cube(cube::Array{Float64,3},can)
    map = Array{Float64}(undef,size(cube)[1],size(cube)[2])
    for ix=1:size(cube)[2]

        for jx=1:size(cube)[1]
            map[jx,ix] = sqrt(StatsBase.moment(cube[jx,ix,can],2))
        end
    end
    rmsavr = StatsBase.mean(skipmissing(map))
    return(map,rmsavr)

end





"""
    rms_cube(cube,can)

Compute the dispersion (e.g. sqrt(second moment order)) on velocity canal given as input of a cube. Return a 2D map with given dispersion on each pixel and the averaged variance accross the map.
"""
function rms_cube(cube::Array{Float64,2},can)
    map = Array{Float64}(undef,size(cube)[1])
    for jx=1:size(cube)[1]
        map[jx] = sqrt(StatsBase.moment(cube[jx,can],2))
    end
    rmsavr = StatsBase.mean(skipmissing(map))
    return(map,rmsavr)

end




end #module