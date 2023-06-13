###################################################################
# STILL WORK IN PROGRESS
# Function used for sorting spectra on a cube. Multiple methods are tested here.  
# Called these functions by calling (or writting on a script ):
#       >include("../src/Functionforsorting.jl")
#       >using .Functionforsorting
#       >output = Functionforsorting.NameOfTheFunction(input)
###################################################################

module Functionforsorting


using FITSIO, MultivariateStats, StatsBase

export divide
export integ
export meanspec
export sort


function divide(RANGE,NUMBEROFMAPS)
    return((RANGE[2]-RANGE[1])/(NUMBEROFMAPS))
end


function integ(spec,RANGE,SIZEINTERVAL,NUMBEROFMAPS,VELOCITYVECTOR)
    DV = abs(VELOCITYVECTOR[2]-VELOCITYVECTOR[1])
    valinteg = Array{Float64}(undef,NUMBEROFMAPS)
    for ix=1:NUMBEROFMAPS-2
        valinteg[ix] = sum(spec[RANGE[1]*ix:RANGE[1]*(ix)+SIZEINTERVAL])*DV
    end
    return(valinteg)
end

function meanspec(cube,VELOCITYVECTOR,DATADIMENSION)
    meanspec = Array{Float64}(undef,DATADIMENSION[3])
    for ix=1:DATADIMENSION[3]
        meanspec[ix] = moment(cube[:,ix],1,0)
    end

    return(meanspec)
end


function smooth(cube::Matrix{Float64},NBINTERVAL,DATADIMENSIONNOMISSING,VELOCITYVECTOR)
    cubesmoothed = Array{Float64}(undef,DATADIMENSIONNOMISSING[1],NBINTERVAL)
    NEWVEL = Array{Float64}(undef,NBINTERVAL)
    INTERVALSIZE = floor(size(VELOCITYVECTOR)[1]/NBINTERVAL) |> Int64
    for ix=1:NBINTERVAL
        for jx=1:DATADIMENSIONNOMISSING[1]
            cubesmoothed[jx,ix] = moment(cube[jx,(ix-1)*INTERVALSIZE+1:(ix)*INTERVALSIZE],1,0)
            NEWVEL[ix] = moment(VELOCITYVECTOR[(ix-1)*INTERVALSIZE+1:ix*INTERVALSIZE],1,0)
        end
    end
    return(cubesmoothed,NEWVEL,INTERVALSIZE)
end

function smooth(cube::Vector{Float64},NBINTERVAL,DATADIMENSIONNOMISSING,VELOCITYVECTOR)
    cubesmoothed = Array{Float64}(undef,NBINTERVAL)
    NEWVEL = Array{Float64}(undef,NBINTERVAL)
    INTERVALSIZE = floor(size(VELOCITYVECTOR)[1]/NBINTERVAL) |> Int64
    for ix=1:NBINTERVAL
        cubesmoothed[ix] = moment(cube[(ix-1)*INTERVALSIZE+1:(ix)*INTERVALSIZE],1,0)
        NEWVEL[ix] = moment(VELOCITYVECTOR[(ix-1)*INTERVALSIZE+1:ix*INTERVALSIZE],1,0)
    end
    return(cubesmoothed,NEWVEL,INTERVALSIZE)
end


function sorts(cubesmoothed,meanspecsmooth,NEWVEL,DATADIMENSION_NOMISSING,INTERVALSIZE,cubesource,BLANK)
    cubeout = Array{Float64}(undef,size(cubesmoothed)[1],DATADIMENSION_NOMISSING[2],size(NEWVEL)[1]+1)
    cubeout .=BLANK
    dif = Array{Float64}(undef,size(cubesmoothed)[1],size(NEWVEL)[1])
    dif .= BLANK
    count = 0
    for ix=1:size(cubesmoothed)[1]
        dif[ix,:] .= meanspecsmooth.-cubesmoothed[ix,:]
        for jx = 1:size(NEWVEL)[1]
            if dif[ix,jx]>0.5
                count += 1
                cubeout[ix,:,jx] = cubesource[ix,:] 
            else
                count += 1
                cubeout[ix,:,end] = cubesource[ix,:]
            end
        end
    end

    return(cubeout)

end


function sort(cubesource,cubereconstructed,DATADIMENSION_NOMISSING,INTEG,RANGE,SIZEINTERVAL,NUMBEROFMAPS,VELOCITYVECTOR,BLANK)
    cubeout = Array{Float64}(undef,size(cubesource)[1],DATADIMENSION_NOMISSING[2],NUMBEROFMAPS)
    cubeout .= BLANK
    dif = Array{Float64}(undef,NUMBEROFMAPS)
    for ix=1:size(cubesource)[1]
        valinteg = integ(cubereconstructed[ix,:],RANGE,SIZEINTERVAL,NUMBEROFMAPS,VELOCITYVECTOR)
        for jx=1:size(INTEG)[1]
            dif[jx] = valinteg[jx]-INTEG[jx]
            #= if valinteg[jx]>=INTEG[jx]
               cubeout[ix,:,jx] = cubesource[ix,:] 
            else   
               cubeout[ix,:,size(INTEG)[1]] = cubesource[ix,:] 
            end =#
        end
        maxipos = findall(x->x==maximum(dif),dif)
        if size(maxipos)[1]>0
            cubeout[ix,:,maxipos[1]] = cubesource[ix,:]
        else
            cubeout[ix,:,NUMBEROFMAPS] = cubesource[ix,:] 
        end
    end

    #= cubeout1 = Array{Float64}(undef,size(cubeout1t)[1],size(VELOCITYVECTOR)[1])
    for ix=1:size(cubeout1)[1]
        cubeout1[ix,:] = cubeout1t[ix] 
    end

    cubeout2 = Array{Float64}(undef,size(cubeout2t)[1],size(VELOCITYVECTOR)[1])
    for ix=1:size(cubeout2)[1]
        cubeout2[ix,:] = cubeout2t[ix] 
    end

    cubeout3 = Array{Float64}(undef,size(cubeout3t)[1],size(VELOCITYVECTOR)[1])
    for ix=1:size(cubeout3)[1]
        cubeout3[ix,:] = cubeout3t[ix] 
    end

    cubeout4 = Array{Float64}(undef,size(cubeout4t)[1],size(VELOCITYVECTOR)[1])
    for ix=1:size(cubeout4)[1]
        cubeout4[ix,:] = cubeout4t[ix] 
    end

    cubeout5 = Array{Float64}(undef,size(cubeout5t)[1],size(VELOCITYVECTOR)[1])
    for ix=1:size(cubeout5)[1]
        cubeout5[ix,:] = cubeout5t[ix] 
    end

    cubeout6 = Array{Float64}(undef,size(cubeout6t)[1],size(VELOCITYVECTOR)[1])
    for ix=1:size(cubeout6)[1]
        cubeout6[ix,:] = cubeout6t[ix] 
    end
    return(cubeout1,cubeout2,cubeout3,cubeout4,cubeout5,cubeout6)
    
end =#
    return(cubeout)
end







function test()
    for ix=1:10
        nmap$ix = Array{Float64}(undef,10)
    end
    return(nmap$ix)
end

end #module