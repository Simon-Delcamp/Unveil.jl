module Spectralwindowopti

include("Data_preparation.jl") #Read and write fits
include("Graphic.jl") #Read and write fits
include("Data_analysis.jl")
using .Data_preparation
using .Graphic
using .Data_analysis

using StatsBase, Plots

export convoptiwind
export optiwind
export otiwind1D # WORK IN PROGRESS


"""
    convoptiwind(cubesource,DATADIMENSION,VELOCITYVECTOR,NOISECAN,BLANK,RANGE,PATHTOSAVE)

DATADIMENSION should be 2D : PV

"""
function convoptiwind(cubesource,DATADIMENSION_NOMISSING,SOURCEDIMENSION,VELOCITYVECTOR,NOISECAN,BLANK,RANGE,PATHTOSAVE,missingplaces2D; ismiss=0)
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
        cubedif = Array{Float64}(undef,DATADIMENSION_NOMISSING[1]) 
        sizenotzero = Array{Float64}(undef,DATADIMENSION_NOMISSING[1]) 
        cubemasked,SIGMAT = optiwind(cubesource,DATADIMENSION_NOMISSING,VELOCITYVECTOR,NOISECAN,BLANK,RANGE[rx])
        #println(SIGMAT)
        cubedif .= (sum(cubesource,dims=2).-sum(cubemasked,dims=2)).*abs(VELOCITYVECTOR[2]-VELOCITYVECTOR[1])
        for px=1:DATADIMENSION_NOMISSING[1]
            sizenotzero[px] = (DATADIMENSION_NOMISSING[2]-size(findall(x->x>0,cubemasked[px,:]))[1])/DATADIMENSION_NOMISSING[2]*100
        end #px
        muc[rx],sigc[rx],gamc[rx],kapc[rx] = Data_analysis.fourmoments(cubedif,dim=1) 
        mus[rx],sigs[rx],gams[rx],kaps[rx] = Data_analysis.fourmoments(sizenotzero,dim=1) 
        if ismiss == 1
            cubedif = Data_preparation.addblank(cubedif,missingplaces2D[:,1],BLANK,(SOURCEDIMENSION[1],SOURCEDIMENSION[2]))
            sizenotzero = Data_preparation.addblank(sizenotzero,missingplaces2D[:,1],BLANK,(SOURCEDIMENSION[1],SOURCEDIMENSION[2]))
        end
        cubedif    = reshape(cubedif,(SOURCEDIMENSION[1],SOURCEDIMENSION[2]))
        sizenotzero = reshape(sizenotzero,(SOURCEDIMENSION[1],SOURCEDIMENSION[2]))

        #Graphic.ploptiwind(cubedif,sizenotzero,SIGMAT,BLANK,RANGE[rx])
        #savefig("$PATHTOSAVE/convoptiwind_$(RANGE[rx]).pdf")
    end #for rx
    metc = Data_analysis.calcmetricOW(muc,sigc,gamc,kapc,abs(VELOCITYVECTOR[2]-VELOCITYVECTOR[1]),SIGMAT)
    mets = Data_analysis.calcmetric(mus,sigs,gams,kaps,abs(VELOCITYVECTOR[2]-VELOCITYVECTOR[1]))
    return(muc,mus,sigc,sigs,gamc,gams,kapc,kaps,metc,mets,SIGMAT)
end #convoptiwind








 
#= function optiwind(cube,DATADIMENSION,cvmap,VELOCITYVECTOR,NOISECAN,BLANK,RANGE)
    maskcube = Array{Float64}(undef,(DATADIMENSION[1],DATADIMENSION[2]))
    maskcube .= 0    # Initiate the cube with 0 everywhere
    integ = 0
    dv = abs(VELOCITYVECTOR[2]-VELOCITYVECTOR[1])
    interv = RANGE   # Increment for integration process given in velocity cannal numbers
    #convert cvmap en indice d'array
    for ix=1:size(maskcube)[1] #previous integ : integral sur un interval assez grand.
        sigmaT = moment(cube[ix,NOISECAN[1]:NOISECAN[2]],2)   # Second moment order of noise canals. Give the noise RMS.
        if cvmap[ix]!=BLANK
            countwhilel = 1    # Counter for integration to the left of the cv value
            countwhiler = 1    # Counter for integration to the right of the cv value
            previousinteg = 0
            arrcvmap = findall(x->x<=cvmap[ix],VELOCITYVECTOR)[1]   # Look for the positions of the cv values in the velocity vector.

            # Optimisation to the left
            integ = sum(cube[ix,arrcvmap-countwhilel*interv:arrcvmap])*dv     # Intensity area from the cv value position to -1 times the interval
            while ((integ-previousinteg)>(sqrt(countwhilel*interv)*sigmaT) && arrcvmap-(countwhilel+1)*interv>1)  # While the newly integrated part is GT sqrt(N*interval)*NoiseRMS AND that the next integrand position is not outside the spectra , THEN do integration on a largest range
                countwhilel += 1
                previousinteg = integ
                integ = sum(cube[ix,arrcvmap-countwhilel*interv:arrcvmap])*dv 
            end

            # Optimisation to the right
            previousinteg = 0
            integ = sum(cube[ix,arrcvmap:arrcvmap+countwhiler*interv])*dv 
            while ((integ-previousinteg)>(sqrt(countwhiler*interv)*sigmaT) && arrcvmap+(countwhiler+1)*interv<size(VELOCITYVECTOR)[1]) # While the newly integrated part is GT sqrt(N*interval)*NoiseRMS AND that the next integrand position is not outside the spectra , THEN do integration on a largest range
                countwhiler += 1
                previousinteg = integ
                integ = sum(cube[ix,arrcvmap:arrcvmap+countwhiler*interv])*dv 
            end
            maskcube[ix,arrcvmap-countwhilel*interv:arrcvmap+countwhiler*interv] .= 1 # When previous conditions are satisfied, put a 1 to every positions inside the limits obtained above. 
        end
    end
    maskcube .= cube.*maskcube   # Apply the mask to the cube
    return(maskcube)
end =#
"""
    optiwind(cube,DATADIMENSION,VELOCITYVECTOR,NOISECAN,BLANK,RANGE) 

INPUT DESCRIPTION
Method of window optimisation. Search for the limits of the signals in each spectrum of the cube. The cube given as input should be in 2Dimensions, with spectrum ordonnated by row. A first CV computation on each spectra should also be given as input : cvmap. Cvmap is in 1D, with each position corresponding to the spectrum of the same position in the 2D cube. The velocity vector should also be given (obtained by Data_preparation.read_fits_ppv). NOISECAN corresponds to a Vector{Int64} of 2 entries, where each entry corresponds to a position in the Velocity Vector of a noise cannal : they have to englobe only noise. Used to compute the noise rms ; ESSENTIAL IN THE COMPUTATION. The last value RANGE is the size of the increment used during the computation : lower it is and longer the code will be. 

METHOD 
A first initialisation of a mask of same size than the input cube, with 0 everywhere. For each spectra, compute the intensity area by velocity increments from each side of its CV value. When the integration is lower than sqrt{N}*NoiseRMS, the limit is kept (with N the number of velocity cannals used for the integration and NoiseRMS the noise rms). When the two limits are obtained, change values of the mask from 0 to 1 inside these limits for each spectra. Last, multiplication of the mask to the cube. 

OUTPUT DESCRIPTION  
Return the cube masked and the mean of the noise accross the map (one sigma).
"""
function optiwind(cube,DATADIMENSION,VELOCITYVECTOR,NOISECAN,BLANK,INTERV)
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
                if bornesup*(INTERV-1)<=0 && bornesup*(INTERV+1)>=DATADIMENSION[2]
                    maskcube[ix,(borneinf)*(INTERV):bornesup*(INTERV)].=1
                elseif bornesup*(INTERV-1)<=0 && bornesup*(INTERV+1)<DATADIMENSION[2]
                    maskcube[ix,(borneinf)*(INTERV):bornesup*(INTERV+1)].=1
                elseif bornesup*(INTERV+1)>=DATADIMENSION[2] && bornesup*(INTERV-1)>0 
                    maskcube[ix,(borneinf)*(INTERV-1):bornesup*(INTERV)].=1
                else 
                    maskcube[ix,(borneinf)*(INTERV-1):bornesup*(INTERV+1)].=1

                end
            end#if
        end#if
    end#for
    return(maskcube.*cube,sigmaT)
end
#= function optiwindTEST(cube,DATADIMENSION,cvmap,VELOCITYVECTOR,NOISECAN,BLANK,INTERV)
    maskcube = Array{Float64}(undef,(DATADIMENSION[1],DATADIMENSION[2]))    # PV mask. Will then multiply the source PV cube.
    maskcube .= 0                                                           # Initiate the cube with 0 everywhere
    dv = abs(VELOCITYVECTOR[2]-VELOCITYVECTOR[1])                                                          # Increment for integration process given in velocity cannal numbers
    wp = Array{Float64}(undef,size(maskcube)[1],2)
    wp[:,1] .= 1
    wp[:,2] .= size(VELOCITYVECTOR)[1]
    #nsize = floor(Int,size(VELOCITYVECTOR)[1]/2) 
    integ = Array{Float64}(undef,size(maskcube)[1],size(VELOCITYVECTOR)[1],2)   # 
    dif = similar(integ)
    integ .= 0
    for ix=1:size(maskcube)[1]  #previous integ : integral sur un interval assez grand.
        #if cvmap[ix]!=BLANK
            arrcvmap = findall(x->x<=cvmap[ix],VELOCITYVECTOR)[1]   # Look for the positions of the cv values in the velocity vector.

            # Detect how many channels on right and left of the cvmap position
            nintright = floor(Int,(size(VELOCITYVECTOR)[1]-arrcvmap)/INTERV)
            nintleft = floor(Int,(arrcvmap-1)/INTERV) 
            (nintleft%2==0) && (nintleft+=1)
            (nintright%2==0) && (nintright+=1)

            nleft = floor(Int,nintleft)    
            nright = floor(Int,nintright)  
            integleft = Array{Float64}(undef,nleft)
            integright = Array{Float64}(undef,nright)
            limleft = Array{Float64}(undef,nleft)
            limright = Array{Float64}(undef,nright)

            sigmaT = moment(cube[ix,NOISECAN[1]:NOISECAN[2]],2,0)      # Second moment order of noise canals. Give the noise RMS.
            #println(sigmaT)
            println(sigmaT)
            # To the left
            for kx=1:nleft
                npos = arrcvmap-kx*INTERV
                (arrcvmap-kx*INTERV <=0) && (npos = 1)
                integleft[kx] = sum(cube[ix,npos:arrcvmap])*dv     # Intensity area from the cv value position to -1 times the interval
                #limleft[kx]=(sqrt(kx*INTERV)*sigmaT)
                integ[ix,kx,1] = integleft[kx]
                #lim[ix,kx,1]=(sigmaT)
                if kx!=nleft && kx!=1
                    dif[ix,kx,1] = integleft[kx]-integleft[kx-1]
                end
            end
            
            #while ((integ[ix,kx,1]-integ[ix,kx-1,1])<(sqrt(kx*INTERV)*sigmaT) && kx>2)
            kx = nleft
            kx = INTERV+1
            
            #while ((abs(dif[ix,kx,1])<abs(sigmaT/5)) && kx>2)
            #    count +=1
            #    kx -= 1
            #    wp[ix,1]=kx
            #end
            while ((kx-INTERV>0) && abs(std(dif[ix,kx-INTERV:kx,1]))<abs(sigmaT/20) && kx>2)
                kx += 1
                wp[ix,1]=kx
            end

            # To the right
            for kx=1:nright
                npos = arrcvmap+kx*INTERV
                (arrcvmap+kx*INTERV >= size(VELOCITYVECTOR)[1]) && (npos=nright)
                integright[kx] = sum(cube[ix,arrcvmap:npos])*dv     # Intensity area from the cv value position to -1 times the interval
                #limright[kx]=(sqrt(kx*INTERV)*sigmaT)
                integ[ix,kx+nleft,2] = integright[kx]
                if kx!=nright && kx!=1
                    dif[ix,kx,2] = integright[kx]-integright[kx-1]
                end
                #lim[ix,kx,2]=(sigmaT)

            end

            kx = nright-2*INTERV-1
            #while ((integ[ix,kx,2]-integ[ix,kx-1,2])<(sqrt(kx*INTERV)*sigmaT) && kx>2)
            while (abs(std(dif[ix,kx:kx+2*INTERV,2]))<abs(sigmaT/20) && kx>2)
                kx -= 1
                wp[ix,2]=kx
            end

            if isnan(arrcvmap-wp[ix,1]*INTERV)==false
                posu = floor(Int,arrcvmap-wp[ix,1]*INTERV)
            else
                posu = arrcvmap
            end
            if isnan(arrcvmap+wp[ix,2]*INTERV)==false
                 posd = floor(Int,arrcvmap+wp[ix,2]*INTERV)
            else
                posd = arrcvmap
            end
            #println(posu,"--",posd)
            (posd>DATADIMENSION[2]) && (posd=DATADIMENSION[2])
            (posu<1               ) && (posu=1)

            maskcube[ix,posu:posd] .= 1
            #(wp[ix,2]>DATADIMENSION[2]) && (wp[ix,2]=DATADIMENSION[2])
            #(wp[ix,1]<1               ) && (wp[ix,1]=1)
            #(wp[ix,1]==wp[ix,2]) && (wp[ix,1]=wp[ix,1]-1)
            #maskcube[ix,floor(Int,wp[ix,1]):floor(Int,wp[ix,2])].= 1
    end
    maskcube .= cube.*maskcube   # Apply the mask to the cube
    return(integ,wp,maskcube,dif)
end
 =#

function optiwind1D(spec,cvmap,VELOCITYVECTOR,NOISECAN,BLANK,RANGE)
    mask = Array{Float64}(undef,size(VELOCITYVECTOR)[1])
    integ = 0
    dv = abs(VELOCITYVECTOR[2]-VELOCITYVECTOR[1])
    #interv = 10
    interv = RANGE #18
    #convert cvmap en indice d'array
    sigmaT = moment(spec[NOISECAN[1]:NOISECAN[2]],2)
    countwhilel = 1
    countwhiler = 1
    previousinteg = 0
    arrcvmap = findall(x->x<=cvmap,VELOCITYVECTOR)[1]
    if arrcvmap==1 || arrcvmap==size(VELOCITYVECTOR)[1]
        arrcvmap = trunc(Int,size(VELOCITYVECTOR)[1])
    end
    # d'abord optimiser à gauche
    integ = sum(spec[arrcvmap-countwhilel*interv:arrcvmap])*dv 
    #while integ-previousinteg>=THRESHOLD && arrcvmap-(countwhilel+1)*interv>1 #integ+previous integ < seuil, continue à additioner
    while (integ-previousinteg)>(sqrt(countwhilel*interv)*sigmaT) && arrcvmap-(countwhilel+1)*interv>1
        countwhilel += 1
        previousinteg = integ
        integ = sum(spec[arrcvmap-countwhilel*interv:arrcvmap])*dv 
        # Condition : on s'arrête quand integ-previousinteg<sqrt(N)*sigma(bruit) , avec N la quantité de canaux que l'on a intégré
    end

    # Puis optimiser à droite
    integ = sum(spec[arrcvmap-countwhilel*interv:arrcvmap+countwhiler*interv])*dv 
    #while integ-previousinteg>=THRESHOLD && arrcvmap+(countwhiler+1)*interv<DATADIMENSION[2] #integ+previous integ < seuil, continue à additioner
    while (integ-previousinteg)>(sqrt(countwhiler*interv)*sigmaT) && (arrcvmap+(countwhiler+1)*interv)<size(VELOCITYVECTOR)[1]
        countwhiler += 1
        previousinteg = integ
        integ = sum(spec[arrcvmap-countwhilel*interv:arrcvmap+countwhiler*interv])*dv 
    end
    mask[arrcvmap-countwhilel*interv:arrcvmap+countwhiler*interv] .= 1 # quand l'addition n'apporte plus rien, on sort les numéros d'intervals que l'on souhaite conserver.
    mask .= spec.*mask
    return(mask)
end










end #module