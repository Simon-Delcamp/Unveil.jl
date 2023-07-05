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
    cubesource should be in 2D (PV), as DATADIMENSION_NOMISSING
"""
# FAIRE EN SORTE DE MOYENNER LES SPECTRES PLUTOT ; PERMETTRA DE LIMITER LES PETITES VARIATIONS 
function bestsnr(cubesource,DATADIMENSION_NOMISSING,VELOCITYVECTOR,NOISECAN)
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
            #wininf[it] = sum(cubesource[sx,1:vx])/sqrt(it+1)/rms
            wininf[it] = sum(cubesource[sx,1:vx])/sqrt(it)/sqrt(rms)
        end
        caninf = findall(x->x==maximum(wininf[1:end-1]),wininf[1:end-1])[1]

        winsup = Array{Float64}(undef,DATADIMENSION_NOMISSING[2])
        it = 0
        for vx=1:DATADIMENSION_NOMISSING[2]-1
            it += 1
            rms = std(cubesource[sx,DATADIMENSION_NOMISSING[2]-vx:DATADIMENSION_NOMISSING[2]])
            winsup[it] = sum(cubesource[sx,DATADIMENSION_NOMISSING[2]-vx:DATADIMENSION_NOMISSING[2]])/sqrt(it)/sqrt(rms)
        end
        cansup = findall(x->x==maximum(winsup[1:end-1]),winsup[1:end-1])[1]
        #println(mean(([sx,DATADIMENSION_NOMISSING[2]-cansup:caninf])))


        
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
                wininf[it] = sum(specsmooth[1:vx])/sqrt(it)/sqrt(rms)
            end
            caninf = findall(x->x==maximum(wininf[1:end-1]),wininf[1:end-1])[1]*2

            winsup = Array{Float64}(undef,nsmooth)
            it = 0
            for vx=1:nsmooth-1
                it += 1
                rms = std(specsmooth[nsmooth-vx:nsmooth])
                winsup[it] = sum(specsmooth[nsmooth-vx:nsmooth])/sqrt(it)/sqrt(rms)
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


        #println(caninf,",",cansup)
    end
    #= WORK ALMOST
    for sx=1:DATADIMENSION_NOMISSING[1]
        sigma = moment(VELOCITYVECTOR,2,aweights(cubesource[sx,:]),0)
        sigmarms = std(cubesource[sx,NOISECAN])
        wininf = Array{Float64}(undef,DATADIMENSION_NOMISSING[2])
        it = 0
        for vx=2:DATADIMENSION_NOMISSING[2]
            it += 1
            rms = std(cubesource[sx,1:vx])
            #wininf[it] = sum(cubesource[sx,1:vx])/sqrt(it+1)/rms
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
        #println(mean(([sx,DATADIMENSION_NOMISSING[2]-cansup:caninf])))


        
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


        #println(caninf,",",cansup)
    end
    =#
    
        #cubeout[sx,:] .= mask[sx,:].*cubesource[sx,:]
    
    
    cubeout .= mask.*cubesource
    return(cubeout,mask)
end #bestsnr
    # WITH SMOOTHING
    #=
    nsmooth = floor(Int,DATADIMENSION_NOMISSING[2]*0.5)
    ninterv = floor(Int,DATADIMENSION_NOMISSING[2]/nsmooth)
    specsmooth = Array{Float64}(undef,nsmooth)
    for sx=1:DATADIMENSION_NOMISSING[1]
        for kx = 1:nsmooth
            specsmooth[kx] = mean(cubesource[sx,(kx-1)*ninterv+1:ninterv*kx])
        end

        wininf = Array{Float64}(undef,nsmooth)
        it = 0
        for vx=2:nsmooth
            it += 1
            rms = std(specsmooth[1:vx])
            wininf[it] = sum(specsmooth[1:vx])/sqrt(vx)/rms
        end
        caninf = findall(x->x==maximum(wininf),wininf)[1]

        winsup = Array{Float64}(undef,nsmooth)
        it = 0
        for vx=1:nsmooth-1
            it += 1
            rms = std(specsmooth[nsmooth-vx:nsmooth])
            winsup[it] = sum(specsmooth[nsmooth-vx:nsmooth])/sqrt(nsmooth-vx)/rms
        end
        cansup = findall(x->x==maximum(winsup),winsup)[1] 
        =#

function autoconv(cubesource,DATADIMENSION_NOMISSING,VELOCITYVECTOR,NOISECAN,BLANK,RANGE; ismiss=0)
    NBRANGE = RANGE[end]-RANGE[1]+1
    minimap = Array{Int}(undef,DATADIMENSION_NOMISSING[1])
    minimap .= 0
    maskcube = Array{Float64}(undef,(DATADIMENSION_NOMISSING[1],size(VELOCITYVECTOR)[1]))
    maskcube .=0
    for px=1:DATADIMENSION_NOMISSING[1]
        specmasked = Array{Float64}(undef,NBRANGE,size(VELOCITYVECTOR)[1])
        specdif = similar(specmasked)
        nbzero = Array{Float64}(undef,NBRANGE)
        nbzerodif = similar(nbzero)
        specintdif = Array{Float64}(undef,NBRANGE)
        muc   = similar(nbzero)
        sigc  = similar(nbzero)
        gamc  = similar(nbzero)
        kapc  = similar(nbzero)
        for rx=1:NBRANGE
            specmasked[rx,:] = swo1D(cubesource[px,:],VELOCITYVECTOR,NOISECAN,BLANK,rx+RANGE[1])
            specdif[rx,:] = (cubesource[px,:].-specmasked[rx,:])
            #nbzerodif[rx] = size(findall(x->x==0,specdif[rx,:]))[1]
            #nbzero[rx] = size(findall(x->x==0,specmasked[rx,:]))[1]

            muc[rx],sigc[rx],gamc[rx],kapc[rx] = Analysis.fourmoments(specdif[rx,:],dim=1) 

        end
        
        met = sqrt.(muc.^2 .+sigc.^2 .+gamc.^2) #+(kapc.-3).^2
        met = Dataprep.replace_nantoblank(met,1e15)
        #println(met)
        #println(met)
        minimap[px] = findall(x->x==minimum(met),met)[1]
        maskcube[px,:] = specmasked[minimap[px],:] 
        #maskcube[px,:] = swo1D(cubesource[px,:],VELOCITYVECTOR,NOISECAN,BLANK,minimap[px])
    end
    return(maskcube,minimap)
end #autoconv


# Tu garde celui qui te done la dif pas plus grande que 2*sigmaT, qui te donne le moins de canaux égaux à 0 et dont le cube  reconstrit contient le plus de 0
#=
function autoconv(cubesource,DATADIMENSION_NOMISSING,VELOCITYVECTOR,NOISECAN,BLANK,RANGE; ismiss=0)
    NBRANGE = RANGE[end]-RANGE[1]+1
    minimap = Array{Int}(undef,DATADIMENSION_NOMISSING[1])
    minimap .= 0
    maskcube = Array{Float64}(undef,(DATADIMENSION_NOMISSING[1],size(VELOCITYVECTOR)[1]))
    maskcube .=0
    #metc = Array{Float64}(undef,NBRANGE)
    #mets = Array{Float64}(undef,NBRANGE)
    for px=1:DATADIMENSION_NOMISSING[1]
        specmasked = Array{Float64}(undef,NBRANGE,size(VELOCITYVECTOR)[1])
        specdif = similar(specmasked)
        nbzero = Array{Float64}(undef,NBRANGE)
        nbzerodif = similar(nbzero)
        specintdif = Array{Float64}(undef,NBRANGE)
        for rx=1:NBRANGE
            specmasked[rx,:] = swo1D(cubesource[px,:],VELOCITYVECTOR,NOISECAN,BLANK,rx+RANGE[1])
            specdif[rx,:] = (cubesource[px,:].-specmasked[rx,:])#.^2
            nbzerodif[rx] = size(findall(x->x==0,specdif[rx,:]))[1]
            nbzero[rx] = size(findall(x->x==0,specmasked[rx,:]))[1]

            #nsize = floor(Int,size(VELOCITYVECTOR)[1]/rx+RANGE[1]) 
            #integ = Array{Float64}(undef,nsize)
            #dinteg = similar(integ) 
            #for vx=1:nsize
            #    integ[vx] = sum(specdif[1:vx*INTERV].^6)
            #    if vx!=1
            #        dinteg[vx] = (integ[vx]-integ[vx-1])
            #    end#if
            #end#for
            specintdif[rx] = (sum(specdif[rx,:])/(rx+RANGE[1]))#/size(specdif[rx,:])[1]#*abs(VELOCITYVECTOR[2]-VELOCITYVECTOR[1])
            #if nbzero[rx]!=0
            #    specintdif[rx] = (specintdif[rx])#*1e-9*nbzero[rx])/(1e-9*(rx+RANGE[1])) #sqrt(nbzero[rx])/
            #else 
            #    #println("yo")
            #    specintdif[rx] = (specintdif[rx])#/(0.00000005*(rx+RANGE[1])))#*sqrt((rx+RANGE[1])))
            #end 
        end
        #println(specintdif)
        sigma= moment(cubesource[px,NOISECAN[1]:NOISECAN[2]],2)
        #println(sigma)
        #println(sigma*sqrt(size(VELOCITYVECTOR)[1]))
        #println(findall(x->x<=sigma*4,specintdif))
        #error()

        dif = findall(x->x<=sigma*10,specintdif)
        if size(dif)[1]!=0
            #println("-------",px)
            #println(size(dif)[1])
            #println(nbzerodif)
            wherezerodif = nbzerodif[dif]
            #println(size(nbzerodif)[1])
            #println(wherezerodif)
            #println(minimum(wherezerodif))
            #println(wherezerodif)
            wherezero    = nbzero[dif]
            #println(wherezero)

            min = findall(x->x==minimum(wherezerodif),wherezerodif)[1]
            max = findall(x->x==maximum(wherezero),wherezero)[1]
            
            if min==max 
                minimap[px] = min 
            else   
                minimap[px] = 100
            end 
        else 
            minimap[px]=100
        end
            #error()
        #many = wherezero[wherezerodif]
        #println(many)
        #minimap[px] = findall(x->x>=sigma,specintdif)[1]+RANGE[1]
        maskcube[px,:] = swo1D(cubesource[px,:],VELOCITYVECTOR,NOISECAN,BLANK,minimap[px])
    end
    return(maskcube,minimap)
end #autoconv
=#
"""
    convswo(cubesource,DATADIMENSION,VELOCITYVECTOR,NOISECAN,BLANK,RANGE,PATHTOSAVE)

DATADIMENSION should be 2D : PV

"""
#=
# On procède spectre par spectre. L'idée c'est d'abord de trouver les range qui donnent une différence d'intégration source-recon faible. Puis parmis ces ranges, de prendre celui qui donne que le nombre de canaux égaux à zero dans la dif du spectre source-recon soit le plus faible possible.
function convswo(cubesource,DATADIMENSION_NOMISSING,SOURCEDIMENSION,VELOCITYVECTOR,NOISECAN,BLANK,RANGE,PATHTOSAVE,missingplaces2D; ismiss=0)
    NBRANGE = RANGE[end]-RANGE[1]+1
    muc = Array{Float64}(undef,NBRANGE)
    mus = Array{Float64}(undef,NBRANGE)
    sigc = Array{Float64}(undef,NBRANGE)
    sigs = Array{Float64}(undef,NBRANGE)
    gamc = Array{Float64}(undef,NBRANGE)
    gams = Array{Float64}(undef,NBRANGE)
    kapc = Array{Float64}(undef,NBRANGE)
    kaps = Array{Float64}(undef,NBRANGE)
    SIGMAT = 0
    posi = Array{Float64}(undef,DATADIMENSION_NOMISSING[1])
    #metc = Array{Float64}(undef,NBRANGE)
    #mets = Array{Float64}(undef,NBRANGE)
    cubeintdif  = Array{Float64}(undef,DATADIMENSION_NOMISSING[1],NBRANGE) 
    sizezero    = similar(cubeintdif)

    for rx=1:NBRANGE
        cubedif     = similar(cubesource)
        sizenotzero = Array{Float64}(undef,DATADIMENSION_NOMISSING[1]) 
        cubemasked,SIGMAT,sigmamap = swo(cubesource,DATADIMENSION_NOMISSING,VELOCITYVECTOR,NOISECAN,BLANK,RANGE[rx])
        cubemasked = Dataprep.blank_equal(cubemasked,BLANK,0)
        #println(SIGMAT)
        cubeintdif[:,rx] .= (sum(cubesource,dims=2).-sum(cubemasked,dims=2)).*abs(VELOCITYVECTOR[2]-VELOCITYVECTOR[1])
        
        for px=1:DATADIMENSION_NOMISSING[1]
            #sizenotzero[px] = (DATADIMENSION_NOMISSING[2]-size(findall(x->x>0,cubemasked[px,:]))[1])/DATADIMENSION_NOMISSING[2]*100
            cubedif[px,:] .=  cubesource[px,:].-cubemasked[px,:]
            sizezero[px,rx] = size(findall(x->x==0,cubedif[px,:]))[1]
        end #px
        #=
        muc[rx],sigc[rx],gamc[rx],kapc[rx] = Analysis.fourmoments(cubeintdif./sizenotzero,dim=1) 
        mus[rx],sigs[rx],gams[rx],kaps[rx] = Analysis.fourmoments(sizezero,dim=1) 
        if ismiss == 1
            cubeintdif = Dataprep.addblank(cubeintdif,missingplaces2D[:,1],BLANK,(SOURCEDIMENSION[1],SOURCEDIMENSION[2]))
            sizenotzero = Dataprep.addblank(sizenotzero,missingplaces2D[:,1],BLANK,(SOURCEDIMENSION[1],SOURCEDIMENSION[2]))
        end
        cubeintdif    = reshape(cubeintdif,(SOURCEDIMENSION[1],SOURCEDIMENSION[2]))
        sizenotzero = reshape(sizezero,(SOURCEDIMENSION[1],SOURCEDIMENSION[2]))

        #Graphic.ploptiwind(cubedif,sizenotzero,SIGMAT,BLANK,RANGE[rx])
        #savefig("$PATHTOSAVE/convoptiwind_$(RANGE[rx]).pdf")
        =#
    end #for rx
    #=
    metc = Analysis.metricOW(muc,sigc,gamc,kapc,abs(VELOCITYVECTOR[2]-VELOCITYVECTOR[1]),SIGMAT)
    mets = Analysis.metricPCA(mus,sigs,gams,kaps,abs(VELOCITYVECTOR[2]-VELOCITYVECTOR[1]))
    return(muc,mus,sigc,sigs,gamc,gams,kapc,kaps,metc,mets,SIGMAT)=#

    for px=1:DATADIMENSION_NOMISSING[1]
        rang = findall(x->x<=SIGMAT*abs(VELOCITYVECTOR[2]-VELOCITYVECTOR[1]),abs.(cubeintdif[px,:]))
        posi[px] = findall(x->x==minimum(sizezero[px,rang]),sizezero[px,rang])[1]
    end

    Dataprep.write_fits("/home/delcamps/Data/Simulated/Pol/Construction/fbm_0-02.fits","posi","test.fits",posi,(DATADIMENSION_NOMISSING[1]),BLANK,finished=true,overwrite=true)

    error()
end #convoptiwind
=#
# WORKED BUT TRY ANOTHER VERSION BY PONDERATION OF THE NUMBER OF CANALS

function convswo(cubesource,DATADIMENSION_NOMISSING,SOURCEDIMENSION,VELOCITYVECTOR,NOISECAN,BLANK,RANGE,PATHTOSAVE,missingplaces2D; ismiss=0)
    NBRANGE = RANGE[end]-RANGE[1]+1
    muc = Array{Float64}(undef,NBRANGE)
    mus = Array{Float64}(undef,NBRANGE)
    sigc = Array{Float64}(undef,NBRANGE)
    sigs = Array{Float64}(undef,NBRANGE)
    gamc = Array{Float64}(undef,NBRANGE)
    gams = Array{Float64}(undef,NBRANGE)
    kapc = Array{Float64}(undef,NBRANGE)
    kaps = Array{Float64}(undef,NBRANGE)
    SIGMAT = 0
    #metc = Array{Float64}(undef,NBRANGE)
    #mets = Array{Float64}(undef,NBRANGE)
    for rx=1:NBRANGE
        cubeintdif = Array{Float64}(undef,DATADIMENSION_NOMISSING[1]) 
        sizenotzero = Array{Float64}(undef,DATADIMENSION_NOMISSING[1]) 
        sizezero= Array{Float64}(undef,DATADIMENSION_NOMISSING[1]) 
        cubemasked,SIGMAT = swo(cubesource,DATADIMENSION_NOMISSING,VELOCITYVECTOR,NOISECAN,BLANK,RANGE[rx])
        cubedif = similar(cubesource)
        cubemasked = Dataprep.blank_equal(cubemasked,BLANK,0)
        #println(SIGMAT)
        #cubeintdif .= (sum(cubesource,dims=2).-sum(cubemasked,dims=2)).*abs(VELOCITYVECTOR[2]-VELOCITYVECTOR[1])
        for px=1:DATADIMENSION_NOMISSING[1]
            sizenotzero[px] = (DATADIMENSION_NOMISSING[2]-size(findall(x->x>0,cubemasked[px,:]))[1])/DATADIMENSION_NOMISSING[2]*100
            cubedif[px,:] .= cubesource[px,:].-cubemasked[px,:]
            cubeintdif[px] = (sum(cubedif[px,:]))*abs(VELOCITYVECTOR[2]-VELOCITYVECTOR[1])
            sizezero[px] =size(findall(x->x==0,cubedif[px,:]))[1]
            #println(sizezero[px])
            if sizezero[px] != 0
                cubeintdif[px] = cubeintdif[px]/sizezero[px]
            end
                
        end #px
        muc[rx],sigc[rx],gamc[rx],kapc[rx] = Analysis.fourmoments(cubeintdif,dim=1) 
        mus[rx],sigs[rx],gams[rx],kaps[rx] = Analysis.fourmoments(sizezero,dim=1) 
        if ismiss == 1
            cubeintdif = Dataprep.addblank(cubeintdif,missingplaces2D[:,1],BLANK,(SOURCEDIMENSION[1],SOURCEDIMENSION[2]))
            sizenotzero = Dataprep.addblank(sizenotzero,missingplaces2D[:,1],BLANK,(SOURCEDIMENSION[1],SOURCEDIMENSION[2]))
        end
        cubeintdif    = reshape(cubeintdif,(SOURCEDIMENSION[1],SOURCEDIMENSION[2]))
        sizenotzero = reshape(sizenotzero,(SOURCEDIMENSION[1],SOURCEDIMENSION[2]))

        #Graphic.ploptiwind(cubeintdif,sizenotzero,SIGMAT,BLANK,RANGE[rx])
        #savefig("$PATHTOSAVE/convoptiwind_$(RANGE[rx]).pdf")
    end #for rx
    metc = Analysis.metricOW(muc,sigc,gamc,kapc,abs(VELOCITYVECTOR[2]-VELOCITYVECTOR[1]),SIGMAT)
    mets = Analysis.metricPCA(mus,sigs,gams,kaps,abs(VELOCITYVECTOR[2]-VELOCITYVECTOR[1]))
    return(muc,mus,sigc,sigs,gamc,gams,kapc,kaps,metc,mets,SIGMAT)
end #convoptiwind



"""
    swo(cube,DATADIMENSION,VELOCITYVECTOR,NOISECAN,BLANK,RANGE) 

INPUT DESCRIPTION
Method of window optimisation. Search for the limits of the signals in each spectrum of the cube. The cube given as input should be in 2Dimensions, with spectrum ordonnated by row. The velocity vector should also be given (obtained by Dataprep.read_fits_ppv). NOISECAN corresponds to a Vector{Int64} of 2 entries, where each entry corresponds to a position in the Velocity Vector of a noise cannal : they have to englobe only noise. Used to compute the noise rms ; ESSENTIAL IN THE COMPUTATION. The last value RANGE is the size of the increment used during the computation : lower it is and longer the code will be. 

METHOD 
A first initialisation of a mask of same size than the input cube, with 0 everywhere. For each spectra, compute the intensity area by velocity increments from each side of its CV value. When the integration is lower than sqrt{N}*NoiseRMS, the limit is kept (with N the number of velocity cannals used for the integration and NoiseRMS the noise rms). When the two limits are obtained, change values of the mask from 0 to 1 inside these limits for each spectra. Last, multiplication of the mask to the cube. 

OUTPUT DESCRIPTION  
Return the cube masked and the mean of the noise accross the map (one sigma).
"""
function swo(cube,DATADIMENSION,VELOCITYVECTOR,NOISECAN,BLANK,INTERV)
    maskcube = Array{Float64}(undef,(DATADIMENSION[1],DATADIMENSION[2]))    # PV mask. Will then multiply the source PV cube.
    maskcube .= 0
    sigmamap = Array{Float64}(undef,DATADIMENSION[1])
    nsize = floor(Int,size(VELOCITYVECTOR)[1]/INTERV) 
    integ = Array{Float64}(undef,(DATADIMENSION[1],nsize))
    #integ = Array{Float64}(undef,nsize)
    dinteg = similar(integ)
    dv = abs(VELOCITYVECTOR[2]-VELOCITYVECTOR[1])
    #lim = 1/sqrt(INTERV)#4/INTERV
    for pix=1:DATADIMENSION[1]
        sigmamap[pix] = moment(cube[pix,NOISECAN[1]:NOISECAN[2]],2)
    end #for pix
    sigmaT = moment(sigmamap,1,0)
    #count = 0
    for ix=1:size(maskcube)[1] 
        if cube[ix,2]!=BLANK
            for vx=1:nsize
                integ[ix,vx] = sum(cube[ix,1:vx*INTERV].^2)*dv
                if vx!=1
                    dinteg[ix,vx] = integ[ix,vx]-integ[ix,vx-1]
                end#if
            end#for
            if size(findall(x->x>=sigmaT*(NOISECAN[2]-NOISECAN[1])*nsize/DATADIMENSION[2],dinteg[ix,2:end]))[1]>1
                #count += 1
                borneinf = findall(x->x>=sigmaT*(NOISECAN[2]-NOISECAN[1])*nsize/DATADIMENSION[2],dinteg[ix,2:end])[1]
                bornesup = findall(x->x>=sigmaT*(NOISECAN[2]-NOISECAN[1])*nsize/DATADIMENSION[2],dinteg[ix,2:end])[end]
                #println(borneinf,",",bornesup)
                #println(borneinf*INTERV,"--",bornesup*INTERV)
                #if borneinf*(INTERV-1)<=0 && bornesup*(INTERV+1)>=DATADIMENSION[2]
                maskcube[ix,(borneinf)*(INTERV):bornesup*(INTERV)].=1
                #elseif borneinf*(INTERV-1)<=0 && bornesup*(INTERV+1)<DATADIMENSION[2]
                #    maskcube[ix,(borneinf)*(INTERV):bornesup*(INTERV+1)].=1
                #elseif bornesup*(INTERV+1)>=DATADIMENSION[2] && borneinf*(INTERV-1)>0 
                #    maskcube[ix,(borneinf)*(INTERV-1):bornesup*(INTERV)].=1
                #else 
                #    maskcube[ix,(borneinf)*(INTERV-1):bornesup*(INTERV+1)].=1
                #end
                #=
                if bornesup*(INTERV-1)<=0 && bornesup*(INTERV+1)>=DATADIMENSION[2]
                    maskcube[ix,(borneinf)*(INTERV):bornesup*(INTERV)].=1
                elseif bornesup*(INTERV-1)<=0 && bornesup*(INTERV+1)<DATADIMENSION[2]
                    maskcube[ix,(borneinf)*(INTERV):bornesup*(INTERV+1)].=1
                elseif bornesup*(INTERV+1)>=DATADIMENSION[2] && bornesup*(INTERV-1)>0 
                    maskcube[ix,(borneinf)*(INTERV-1):bornesup*(INTERV)].=1
                else 
                    maskcube[ix,(borneinf)*(INTERV-1):bornesup*(INTERV+1)].=1
                end
                =#
            end#if
        end#if
    end#for
    #println(count)
    #println(DATADIMENSION[1])
    #println(sigmaT)

    maskcube = Dataprep.blank_equal(maskcube.*cube,0,BLANK )
    return(maskcube,sigmaT,sigmamap)
end


#=
function swo(cube,DATADIMENSION,VELOCITYVECTOR,NOISECAN,BLANK,INTERV)
    maskcube = Array{Float64}(undef,(DATADIMENSION[1],DATADIMENSION[2]))    # PV mask. Will then multiply the source PV cube.
    maskcube .= 0
    sigmamap = Array{Float64}(undef,DATADIMENSION[1])
    nsize = floor(Int,size(VELOCITYVECTOR)[1]/INTERV) 
    integ = Array{Float64}(undef,(DATADIMENSION[1],nsize))
    #integ = Array{Float64}(undef,nsize)
    dinteg = similar(integ)
    dv = abs(VELOCITYVECTOR[2]-VELOCITYVECTOR[1])
    #lim = 1/sqrt(INTERV)#4/INTERV
    for pix=1:DATADIMENSION[1]
        sigmamap[pix] = moment(cube[pix,NOISECAN[1]:NOISECAN[2]],2)
    end #for pix
    sigmaT = moment(sigmamap,1,0)
    for ix=1:size(maskcube)[1] 

        if cube[ix,2]!=BLANK
            
            for vx=1:nsize
                integ[ix,vx] = sum(cube[ix,1:vx*INTERV].^2)*dv
                if vx!=1
                    dinteg[ix,vx] = integ[ix,vx]-integ[ix,vx-1]
                end#if
            end#for
            if size(findall(x->x>=sigmaT^2*INTERV*dv*(NOISECAN[2]-NOISECAN[1]),dinteg[ix,2:end]))[1]>1
                borneinf = findall(x->x>=sigmaT^2*INTERV*dv*(NOISECAN[2]-NOISECAN[1]),dinteg[ix,2:end])[1]
                bornesup = findall(x->x>=sigmaT^2*INTERV*dv*(NOISECAN[2]-NOISECAN[1]),dinteg[ix,2:end])[end]
                #println(borneinf*INTERV,"--",bornesup*INTERV)
                #if borneinf*(INTERV-1)<=0 && bornesup*(INTERV+1)>=DATADIMENSION[2]
                    maskcube[ix,(borneinf)*(INTERV):bornesup*(INTERV)].=1
                #elseif borneinf*(INTERV-1)<=0 && bornesup*(INTERV+1)<DATADIMENSION[2]
                #    maskcube[ix,(borneinf)*(INTERV):bornesup*(INTERV+1)].=1
                #elseif bornesup*(INTERV+1)>=DATADIMENSION[2] && borneinf*(INTERV-1)>0 
                #    maskcube[ix,(borneinf)*(INTERV-1):bornesup*(INTERV)].=1
                #else 
                #    maskcube[ix,(borneinf)*(INTERV-1):bornesup*(INTERV+1)].=1
                #end
                #=
                if bornesup*(INTERV-1)<=0 && bornesup*(INTERV+1)>=DATADIMENSION[2]
                    maskcube[ix,(borneinf)*(INTERV):bornesup*(INTERV)].=1
                elseif bornesup*(INTERV-1)<=0 && bornesup*(INTERV+1)<DATADIMENSION[2]
                    maskcube[ix,(borneinf)*(INTERV):bornesup*(INTERV+1)].=1
                elseif bornesup*(INTERV+1)>=DATADIMENSION[2] && bornesup*(INTERV-1)>0 
                    maskcube[ix,(borneinf)*(INTERV-1):bornesup*(INTERV)].=1
                else 
                    maskcube[ix,(borneinf)*(INTERV-1):bornesup*(INTERV+1)].=1
                end
                =#
            end#if
        end#if
    end#for
    maskcube = Dataprep.blank_equal(maskcube.*cube,0,BLANK )
    return(maskcube,sigmaT,sigmamap)
end
=#

function swo1D(spec,VELOCITYVECTOR,NOISECAN,BLANK,INTERV)
    mask = Array{Float64}(undef,size(VELOCITYVECTOR)[1])
    mask .= 0
    nsize = floor(Int,size(VELOCITYVECTOR)[1]/INTERV) 
    integ = Array{Float64}(undef,nsize)
    dinteg = similar(integ)
    dv = abs(VELOCITYVECTOR[2]-VELOCITYVECTOR[1])
    sigma= moment(spec[NOISECAN[1]:NOISECAN[2]],2)
    if spec[2]!=BLANK
        specn = spec./maximum(spec)
        for vx=1:nsize
        #    integ[vx] = sum(spec[1:vx*INTERV].^2)*dv
            #integ[vx] = sum(spec[1:vx*INTERV].^8)
            integ[vx] = sum(specn[1:vx*INTERV])
            if vx!=1
                dinteg[vx] = (integ[vx]-integ[vx-1])
            end#if
        end#for
        #println(dinteg)
        #error()
        #if size(findall(x->x>=sigma*(NOISECAN[2]-NOISECAN[1])*nsize/INTERV/size(VELOCITYVECTOR)[1],dinteg[2:end]))[1]>1
        #    #count += 1
        #    borneinf = findall(x->x>=sigma*(NOISECAN[2]-NOISECAN[1])*nsize/INTERV/size(VELOCITYVECTOR)[1],dinteg[2:end])[1]
        #    bornesup = findall(x->x>=sigma*(NOISECAN[2]-NOISECAN[1])*nsize/INTERV/size(VELOCITYVECTOR)[1],dinteg[2:end])[end]
        #if size(findall(x->x>=sigma^4*sqrt(INTERV),dinteg[1:end]))[1]>1
        #    borneinf = findall(x->x>=sigma^4*sqrt(INTERV),dinteg)[1]
        #    bornesup = findall(x->x>=sigma^4*sqrt(INTERV),dinteg)[end]
        if size(findall(x->x>=sigma,dinteg[1:end]))[1]>1
            borneinf = findall(x->x>=sigma,dinteg)[1]
            bornesup = findall(x->x>=sigma,dinteg)[end]
            mask[(borneinf-1)*(INTERV)+1:(bornesup-1)*(INTERV)].=1
            ##if bornesup!=nsize
            ##    println(bornesup,",",nsize)
            ##end
            ##println(borneinf,",",bornesup,",",nsize,",",nsize*INTERV,",",size(VELOCITYVECTOR)[1])
        else
            borneinf = 1
            bornesup = nsize
            mask[(borneinf-1)*(INTERV)+1:(bornesup-1)*(INTERV)].=1

        end#if
        
    end#if

    mask .= spec.*mask
    mask = Dataprep.blank_equal(mask,0,BLANK )

    return(mask)
end

end #module