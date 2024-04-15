module SWO

include("Dataprep.jl")
include("Analysis.jl")

using .Dataprep
using .Analysis

using StatsBase

export convswo
export swo
export swo1D

using Plots


"""
   swo(cubesource,DATADIMENSION_NOMISSING,VELOCITYVECTOR,NOISECAN)

Method SWO (Spectral Window Optimisation). For each spectra of the *cubesource*, look for the maximum SNR given by :
```math
SNR(v_\\text{i})=\\frac{\\sum_{\\text{i}=1}^\\text{m}T_\\text{i}}{\\sqrt{m\\sigma_{1:m}}}
```
In this equation, ``T_\\text{i}`` is the intensity at velocity channel ``v_\\text{i}``, and ``\\sigma_{1:m}`` the dispersion computed between velocity channels 1 to m. The SNR is computed in increasing velocity channels, then in decreasing velocity channels. Each maxima of these two computations will give one window limit containing the emission.
    
INPUT : *cubesource*, the cube whom you want to compute the SWO. It should be in 2D (PxP,V), as *DATADIMENSION_NOMISSING*, which gives the total dimensions of the cube. *VELOCITYVECTOR* is the velocity vector, computed and given as output of the function **Dataprep.read\\_fits\\_ppv**. *NOISECAN* are some of the velocity channels of the cube where there is only noise. Used to compute the dispersion of the noise of *cubesource*.
    
OUTPUT : [1] Cubesource with SWO mask applied (inside mask windows are equal to cubesource, outside is 0) 
         [2] Mask of SWO method computed on cubesource
"""
function swo(cubesource,DATADIMENSION_NOMISSING,VELOCITYVECTOR,NOISECAN)
    cubeout = similar(cubesource)
    cubeout .= 0
    mask = Array{Float64}(undef,DATADIMENSION_NOMISSING[1],DATADIMENSION_NOMISSING[2])
    mask .= 0
    step = 1
    for sx=1:DATADIMENSION_NOMISSING[1]
        sigma = moment(VELOCITYVECTOR,2,aweights(cubesource[sx,:]),0)
        sigmarms = std(cubesource[sx,NOISECAN])
        wininf = Array{Float64}(undef,DATADIMENSION_NOMISSING[2])
        it = 0
        for vx=2:DATADIMENSION_NOMISSING[2]
            it += 1
            rms = std(cubesource[sx,1:vx])
            #wininf[it] = sum(cubesource[sx,1:vx])/sqrt(it)/sqrt(rms)
            wininf[it] = sum(cubesource[sx,1:vx])/sqrt(it)/sqrt(sigmarms)
        end
        caninf = findall(x->x==maximum(wininf[1:end-1]),wininf[1:end-1])[1]

        winsup = Array{Float64}(undef,DATADIMENSION_NOMISSING[2])
        it = 0
        for vx=1:DATADIMENSION_NOMISSING[2]-1
            it += 1
            rms = std(cubesource[sx,DATADIMENSION_NOMISSING[2]-vx:DATADIMENSION_NOMISSING[2]])
            #winsup[it] = sum(cubesource[sx,DATADIMENSION_NOMISSING[2]-vx:DATADIMENSION_NOMISSING[2]])/sqrt(it)/sqrt(rms)
            winsup[it] = sum(cubesource[sx,DATADIMENSION_NOMISSING[2]-vx:DATADIMENSION_NOMISSING[2]])/sqrt(it)/sqrt(sigmarms)
        end
        cansup = findall(x->x==maximum(winsup[1:end-1]),winsup[1:end-1])[1]
        
        if ((DATADIMENSION_NOMISSING[2]-cansup)>caninf) || (abs((DATADIMENSION_NOMISSING[2]-cansup)-caninf)*abs(VELOCITYVECTOR[2]-VELOCITYVECTOR[1])<0.1*sigma) || (std(cubesource[sx,DATADIMENSION_NOMISSING[2]-cansup:caninf])<2*sigmarms) #|| (mean(cubesource[sx,DATADIMENSION_NOMISSING[2]-cansup:caninf])<2*sigmarms)
            nsmooth = floor(Int,DATADIMENSION_NOMISSING[2]*0.1)
            ninterv = floor(Int,DATADIMENSION_NOMISSING[2]/nsmooth)
            specsmooth = Array{Float64}(undef,nsmooth)
            for kx = 1:nsmooth
                specsmooth[kx] = mean(cubesource[sx,(kx-1)*ninterv+1:ninterv*kx])
            end

            wininf = Array{Float64}(undef,nsmooth)
            it = 0
            for vx=2:nsmooth
                it += 1
                rms = std(specsmooth[1:vx])
                #wininf[it] = sum(specsmooth[1:vx])/sqrt(it)/sqrt(rms)
                wininf[it] = sum(specsmooth[1:vx])/sqrt(it)/sqrt(sigmarms)
            end
            caninf = findall(x->x==maximum(wininf[1:end-1]),wininf[1:end-1])[1]*2

            winsup = Array{Float64}(undef,nsmooth)
            it = 0
            for vx=1:nsmooth-1
                it += 1
                rms = std(specsmooth[nsmooth-vx:nsmooth])
                #winsup[it] = sum(specsmooth[nsmooth-vx:nsmooth])/sqrt(it)/sqrt(rms)
                winsup[it] = sum(specsmooth[nsmooth-vx:nsmooth])/sqrt(it)/sqrt(sigmarms)
            end
            cansup = findall(x->x==maximum(winsup[1:end-1]),winsup[1:end-1])[1]*2
        end
        
        if (DATADIMENSION_NOMISSING[2]-cansup)>caninf 
            caninf = DATADIMENSION_NOMISSING[2]
            cansup = DATADIMENSION_NOMISSING[2]-1
        end

        if abs((DATADIMENSION_NOMISSING[2]-cansup)-caninf)*abs(VELOCITYVECTOR[2]-VELOCITYVECTOR[1])<0.1*sigma || (std(cubesource[sx,DATADIMENSION_NOMISSING[2]-cansup:caninf])<2*sigmarms) || (mean(cubesource[sx,DATADIMENSION_NOMISSING[2]-cansup:caninf])<2*sigmarms)
            caninf = DATADIMENSION_NOMISSING[2]
            cansup = DATADIMENSION_NOMISSING[2]-1
            
        end


        mask[sx,DATADIMENSION_NOMISSING[2]-cansup:caninf] .= 1

    end
    cubeout .= mask.*cubesource
    return(cubeout,mask)
end #bestsnr
   


function newswo(cubesource,DATADIMENSION_NOMISSING,VELOCITYVECTOR,NOISECAN)
    cubeout = similar(cubesource)
    cubeout .= 0
    
    mask = Array{Float64}(undef,DATADIMENSION_NOMISSING[1],DATADIMENSION_NOMISSING[2])
    mask .= 0

    posimap = Array{Float64}(undef,size(cubesource)[1],2)
    posimap .= 0
    step = 1
    for sx=1:DATADIMENSION_NOMISSING[1]
        sigma = moment(VELOCITYVECTOR,2,aweights(cubesource[sx,:]),0)
        sigmarms = std(cubesource[sx,NOISECAN])
        wininf = Array{Float64}(undef,DATADIMENSION_NOMISSING[2])
        it = 0
        for vx=2:DATADIMENSION_NOMISSING[2]
            it += 1
            rms = std(cubesource[sx,1:vx])
            #wininf[it] = sum(cubesource[sx,1:vx])/sqrt(it)/sqrt(rms)
            wininf[it] = sum(cubesource[sx,1:vx])/sqrt(it)/sqrt(sigmarms)
        end
        caninf = findall(x->x==maximum(wininf[1:end-1]),wininf[1:end-1])[1]

        winsup = Array{Float64}(undef,DATADIMENSION_NOMISSING[2])
        it = 0
        for vx=1:DATADIMENSION_NOMISSING[2]-1
            it += 1
            rms = std(cubesource[sx,DATADIMENSION_NOMISSING[2]-vx:DATADIMENSION_NOMISSING[2]])
            #winsup[it] = sum(cubesource[sx,DATADIMENSION_NOMISSING[2]-vx:DATADIMENSION_NOMISSING[2]])/sqrt(it)/sqrt(rms)
            winsup[it] = sum(cubesource[sx,DATADIMENSION_NOMISSING[2]-vx:DATADIMENSION_NOMISSING[2]])/sqrt(it)/sqrt(sigmarms)
        end
        cansup = findall(x->x==maximum(winsup[1:end-1]),winsup[1:end-1])[1]
        
        if ((DATADIMENSION_NOMISSING[2]-cansup)>caninf) || (abs((DATADIMENSION_NOMISSING[2]-cansup)-caninf)*abs(VELOCITYVECTOR[2]-VELOCITYVECTOR[1])<0.1*sigma) || (std(cubesource[sx,DATADIMENSION_NOMISSING[2]-cansup:caninf])<2*sigmarms) || (mean(cubesource[sx,DATADIMENSION_NOMISSING[2]-cansup:caninf])<2*sigmarms)
            nsmooth = floor(Int,DATADIMENSION_NOMISSING[2]*0.1)
            ninterv = floor(Int,DATADIMENSION_NOMISSING[2]/nsmooth)
            specsmooth = Array{Float64}(undef,nsmooth)
            for kx = 1:nsmooth
                specsmooth[kx] = mean(cubesource[sx,(kx-1)*ninterv+1:ninterv*kx])
            end

            wininf = Array{Float64}(undef,nsmooth)
            it = 0
            for vx=2:nsmooth
                it += 1
                rms = std(specsmooth[1:vx])
                #wininf[it] = sum(specsmooth[1:vx])/sqrt(it)/sqrt(rms)
                wininf[it] = sum(specsmooth[1:vx])/sqrt(it)/sqrt(sigmarms)
            end
            caninf = findall(x->x==maximum(wininf[1:end-1]),wininf[1:end-1])[1]*2

            winsup = Array{Float64}(undef,nsmooth)
            it = 0
            for vx=1:nsmooth-1
                it += 1
                rms = std(specsmooth[nsmooth-vx:nsmooth])
                #winsup[it] = sum(specsmooth[nsmooth-vx:nsmooth])/sqrt(it)/sqrt(rms)
                winsup[it] = sum(specsmooth[nsmooth-vx:nsmooth])/sqrt(it)/sqrt(sigmarms)
            end
            cansup = findall(x->x==maximum(winsup[1:end-1]),winsup[1:end-1])[1]*2
        end
        
        if (DATADIMENSION_NOMISSING[2]-cansup)>caninf 
            caninf = DATADIMENSION_NOMISSING[2]
            cansup = DATADIMENSION_NOMISSING[2]-1
        end

        if abs((DATADIMENSION_NOMISSING[2]-cansup)-caninf)*abs(VELOCITYVECTOR[2]-VELOCITYVECTOR[1])<0.1*sigma || (std(cubesource[sx,DATADIMENSION_NOMISSING[2]-cansup:caninf])<2*sigmarms) || (mean(cubesource[sx,DATADIMENSION_NOMISSING[2]-cansup:caninf])<2*sigmarms)
            caninf = DATADIMENSION_NOMISSING[2]
            cansup = DATADIMENSION_NOMISSING[2]-1
            
        end


        mask[sx,DATADIMENSION_NOMISSING[2]-cansup:caninf] .= 1
        posimap[sx,1] = DATADIMENSION_NOMISSING[2]-cansup
        posimap[sx,2] = caninf

    end
    cubeout .= mask.*cubesource
    return(cubeout,mask,posimap)
end #bestsnr


"""
   petysnr(cubesource,DATADIMENSION_NOMISSING,VELOCITYVECTOR,NOISECAN)

Method from Pety+2003, similar to the SWO method and used in the same way.

INPUT : *cubesource*, the cube whom you want to compute the SWO. It should be in 2D (PxP,V), as *DATADIMENSION_NOMISSING*, which gives the total dimensions of the cube. *VELOCITYVECTOR* is the velocity vector, computed and given as output of the function **Dataprep.read\\_fits\\_ppv**. *NOISECAN* are some of the velocity channels of the cube where there is only noise. Used to compute the dispersion of the noise of *cubesource*.
    
OUTPUT : [1] Cubesource with SWO mask applied (inside mask windows are equal to cubesource, outside is 0) 
         [2] Mask of SWO method computed on cubesource
"""
function petysnr(cubesource,DATADIMENSION_NOMISSING,VELOCITYVECTOR,NOISECAN)
    cubeout = similar(cubesource)
    cubeout .= 0
    mask = Array{Float64}(undef,DATADIMENSION_NOMISSING[1],DATADIMENSION_NOMISSING[2])
    mask .= 0
    step = 1
    for sx=1:DATADIMENSION_NOMISSING[1]
        sigma = moment(VELOCITYVECTOR,2,aweights(cubesource[sx,:]),0)
        sigmarms = std(cubesource[sx,NOISECAN])
        wininf = Array{Float64}(undef,DATADIMENSION_NOMISSING[2])
        it = 0
        for vx=2:DATADIMENSION_NOMISSING[2]
            it += 1
            rms = std(cubesource[sx,1:vx])
            wininf[it] = sum(cubesource[sx,1:vx])/sqrt(it)/sigmarms
        end
        caninf = findall(x->x==maximum(wininf[1:end-1]),wininf[1:end-1])[1]

        winsup = Array{Float64}(undef,DATADIMENSION_NOMISSING[2])
        it = 0
        for vx=1:DATADIMENSION_NOMISSING[2]-1
            it += 1
            rms = std(cubesource[sx,DATADIMENSION_NOMISSING[2]-vx:DATADIMENSION_NOMISSING[2]])
            winsup[it] = sum(cubesource[sx,DATADIMENSION_NOMISSING[2]-vx:DATADIMENSION_NOMISSING[2]])/sqrt(it)/sigmarms
        end
        cansup = findall(x->x==maximum(winsup[1:end-1]),winsup[1:end-1])[1]
        
        if ((DATADIMENSION_NOMISSING[2]-cansup)>caninf) || (abs((DATADIMENSION_NOMISSING[2]-cansup)-caninf)*abs(VELOCITYVECTOR[2]-VELOCITYVECTOR[1])<0.1*sigma) || (std(cubesource[sx,DATADIMENSION_NOMISSING[2]-cansup:caninf])<2*sigmarms) #|| (mean(cubesource[sx,DATADIMENSION_NOMISSING[2]-cansup:caninf])<2*sigmarms)
            nsmooth = floor(Int,DATADIMENSION_NOMISSING[2]*0.1)
            ninterv = floor(Int,DATADIMENSION_NOMISSING[2]/nsmooth)
            specsmooth = Array{Float64}(undef,nsmooth)
            for kx = 1:nsmooth
                specsmooth[kx] = mean(cubesource[sx,(kx-1)*ninterv+1:ninterv*kx])
            end

            wininf = Array{Float64}(undef,nsmooth)
            it = 0
            for vx=2:nsmooth
                it += 1
                rms = std(specsmooth[1:vx])
                wininf[it] = sum(specsmooth[1:vx])/sqrt(it)/sigmarms
            end
            caninf = findall(x->x==maximum(wininf[1:end-1]),wininf[1:end-1])[1]*2

            winsup = Array{Float64}(undef,nsmooth)
            it = 0
            for vx=1:nsmooth-1
                it += 1
                rms = std(specsmooth[nsmooth-vx:nsmooth])
                winsup[it] = sum(specsmooth[nsmooth-vx:nsmooth])/sqrt(it)/sigmarms
            end
            cansup = findall(x->x==maximum(winsup[1:end-1]),winsup[1:end-1])[1]*2
        end
        
        if (DATADIMENSION_NOMISSING[2]-cansup)>caninf 
            caninf = DATADIMENSION_NOMISSING[2]
            cansup = DATADIMENSION_NOMISSING[2]-1
        end

        if abs((DATADIMENSION_NOMISSING[2]-cansup)-caninf)*abs(VELOCITYVECTOR[2]-VELOCITYVECTOR[1])<0.1*sigma || (std(cubesource[sx,DATADIMENSION_NOMISSING[2]-cansup:caninf])<2*sigmarms) || (mean(cubesource[sx,DATADIMENSION_NOMISSING[2]-cansup:caninf])<2*sigmarms)
            caninf = DATADIMENSION_NOMISSING[2]
            cansup = DATADIMENSION_NOMISSING[2]-1
        end


        mask[sx,DATADIMENSION_NOMISSING[2]-cansup:caninf] .= 1

    end
    cubeout .= mask.*cubesource
    return(cubeout,mask)
end #bestsnr

#=
function petysnr(cubesource,DATADIMENSION_NOMISSING,VELOCITYVECTOR,NOISECAN)
    cubeout = similar(cubesource)
    cubeout .= 0
    mask = Array{Float64}(undef,DATADIMENSION_NOMISSING[1],DATADIMENSION_NOMISSING[2])
    mask .= 0
    step = 1
    for sx=1:DATADIMENSION_NOMISSING[1]
        sigma = moment(VELOCITYVECTOR,2,aweights(cubesource[sx,:]),0)
        sigmarms = std(cubesource[sx,NOISECAN])
        wininf = Array{Float64}(undef,DATADIMENSION_NOMISSING[2])
        it = 0
        for vx=2:DATADIMENSION_NOMISSING[2]
            it += 1
            rms = std(cubesource[sx,1:vx])
            wininf[it] = sum(cubesource[sx,1:vx])/sqrt(it)/rms
        end
        caninf = findall(x->x==maximum(wininf[1:end-1]),wininf[1:end-1])[1]

        winsup = Array{Float64}(undef,DATADIMENSION_NOMISSING[2])
        it = 0
        for vx=1:DATADIMENSION_NOMISSING[2]-1
            it += 1
            rms = std(cubesource[sx,DATADIMENSION_NOMISSING[2]-vx:DATADIMENSION_NOMISSING[2]])
            winsup[it] = sum(cubesource[sx,DATADIMENSION_NOMISSING[2]-vx:DATADIMENSION_NOMISSING[2]])/sqrt(it)/rms
        end
        cansup = findall(x->x==maximum(winsup[1:end-1]),winsup[1:end-1])[1]
        
        if ((DATADIMENSION_NOMISSING[2]-cansup)>caninf) || (abs((DATADIMENSION_NOMISSING[2]-cansup)-caninf)*abs(VELOCITYVECTOR[2]-VELOCITYVECTOR[1])<0.1*sigma) || (std(cubesource[sx,DATADIMENSION_NOMISSING[2]-cansup:caninf])<2*sigmarms) #|| (mean(cubesource[sx,DATADIMENSION_NOMISSING[2]-cansup:caninf])<2*sigmarms)
            nsmooth = floor(Int,DATADIMENSION_NOMISSING[2]*0.1)
            ninterv = floor(Int,DATADIMENSION_NOMISSING[2]/nsmooth)
            specsmooth = Array{Float64}(undef,nsmooth)
            for kx = 1:nsmooth
                specsmooth[kx] = mean(cubesource[sx,(kx-1)*ninterv+1:ninterv*kx])
            end

            wininf = Array{Float64}(undef,nsmooth)
            it = 0
            for vx=2:nsmooth
                it += 1
                rms = std(specsmooth[1:vx])
                wininf[it] = sum(specsmooth[1:vx])/sqrt(it)/rms
            end
            caninf = findall(x->x==maximum(wininf[1:end-1]),wininf[1:end-1])[1]*2

            winsup = Array{Float64}(undef,nsmooth)
            it = 0
            for vx=1:nsmooth-1
                it += 1
                rms = std(specsmooth[nsmooth-vx:nsmooth])
                winsup[it] = sum(specsmooth[nsmooth-vx:nsmooth])/sqrt(it)/rms
            end
            cansup = findall(x->x==maximum(winsup[1:end-1]),winsup[1:end-1])[1]*2
        end
        
        if (DATADIMENSION_NOMISSING[2]-cansup)>caninf 
            caninf = DATADIMENSION_NOMISSING[2]
            cansup = DATADIMENSION_NOMISSING[2]-1
        end

        if abs((DATADIMENSION_NOMISSING[2]-cansup)-caninf)*abs(VELOCITYVECTOR[2]-VELOCITYVECTOR[1])<0.1*sigma || (std(cubesource[sx,DATADIMENSION_NOMISSING[2]-cansup:caninf])<2*sigmarms) || (mean(cubesource[sx,DATADIMENSION_NOMISSING[2]-cansup:caninf])<2*sigmarms)
            caninf = DATADIMENSION_NOMISSING[2]
            cansup = DATADIMENSION_NOMISSING[2]-1
        end


        mask[sx,DATADIMENSION_NOMISSING[2]-cansup:caninf] .= 1

    end
    cubeout .= mask.*cubesource
    return(cubeout,mask)
end #bestsnr
   =#

end #module