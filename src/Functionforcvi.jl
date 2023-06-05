###################################################################
# Function used for the CV and CVI computations processes. 
# Called these functions by calling (or writting on a script ):
#       >include("../src/Functionforcvi.jl")
#       >using .Functionforcvi
#       >output = Functionforcvi.NameOfTheFunction(input)
###################################################################


module Functionforcvi

include("Data_preparation.jl") #Read and write fits


using FITSIO, MultivariateStats, StatsPlots, Plots, ShiftedArrays, Distributions,StatsBase
using Mmap
using .Data_preparation

export construct_cvmap
export construct_cvimap
export construct_cvimap_abs
export cv_increment!
export cv_increment
export moment_four
export moment_one
export moment_one_field
export moment_three
export moment_two
export moment_two_field
export multiple_moment


"""
    construct_cvmap(xyarr,velvector)

Will construct a cv map by calculating the first velocity moment order of the array and replace the inf values by NaN. xyarr have to be in 2D (pixel*pixel), while velvector have to be in 1D (velocities).
"""
function construct_cvmap(xyarr,velvector)
    cvpc = reshape(moment_one_field(xyarr,velvector),size(xyarr)[1],size(xyarr)[2])
    cvpc = Data_preparation.replace_inf_in_nan(cvpc)
    #Conversion type, in case of missing values in the original dataset
    cvpc = convert(Array{Float64},cvpc)
    return(cvpc)
end


"""
    construct_cvimap!(xyarr,Lag::Vector{Int64},mapdim,cvi_averaged_alllag::Array{Union{Missing, Float64},2},cvi_allangle_alllag::Array{Union{Missing, Float64}, 3},cvi_allangle::Array{Union{Missing,Float64},3}; diff="relative",keepmissing=true)

Will construct a cvi map using preallocated array. Need an array preallocated for the storage of the average on all cvi angle (cvimean), and another for the cvi calculated with all angle and lag (cvi_allangle_alllag). The cvi map is constructed by taking the mean of all rotations of the cv increment calculation at each pixel. The Lag is the increment. xyarr have to be in 2D (pixel*pixel). Mapdim is the dimension of your 2Dmap. The differences can be absolute or relative. 
"""
function construct_cvimap!(xyarr,Lag::Vector{Int64},nangle,mapdim,cvi_averaged_alllag::Array{Union{Missing, Float64},2},cvi_allangle_alllag::Array{Union{Missing, Float64}, 3},cvi_allangle::Array{Union{Missing,Float64},3}; diff="relative",keepmissing=true)
    nangle = floor(Int,2pi ./(atan.(1 ./Lag)))
    
    cvi_allangle_alllag = cv_increment!(xyarr,Lag,nangle,cvi_allangle_alllag,cvi_allangle,diff=diff)
    @inbounds @views for lagstep=1:size(Lag)[1]
        if keepmissing==true
            missing1D = findall(ismissing,cvi_allangle_alllag[:,:,lagstep])
            cvi_allangle_alllag[missing1D,lagstep] .= missing
        end
        @inbounds @views for pix=1:size(cvi_allangle_alllag)[1]
            keepmissing==false && (cvi_averaged_alllag[pix,lagstep] = mean(skipmissing(cvi_allangle_alllag[pix,:,lagstep])))
            keepmissing==true  && (cvi_averaged_alllag[pix,lagstep] = mean(cvi_allangle_alllag[pix,:,lagstep]))
        end
    end
    cvi_averaged_alllag = reshape(cvi_averaged_alllag,mapdim[1],mapdim[2],size(Lag)[1])
    return(cvi_averaged_alllag,cvi_allangle_alllag,nangle)
end



"""
    construct_cvimap!(xyarr,Lag::Int64,mapdim,cvi_averaged::Array{Union{Missing, Float64}},cvi_allangle::Array{Union{Missing,Float64},3}; diff="relative",keepmissing=true) 

Will construct a cvi map using preallocated array. Need an array preallocated for the storage of the average on all cvi angle (cvi_averaged), and another for the cvi calculated with all angle and lag (cvi_allangle). The cvi map is constructed by taking the mean of all rotations of the cv increment calculation at each pixel. The Lag is the increment. xyarr have to be in 2D (pixel*pixel). Mapdim is the dimension of your 2Dmap. The differences can be absolute or relative. 
"""
function construct_cvimap!(xyarr,Lag::Int64,mapdim,cvi_averaged::Array{Union{Missing, Float64}},cvi_allangle::Array{Union{Missing,Float64},3}; diff="relative",keepmissing=true) 
    nangle = floor.(Int,2pi ./(atan.(1 ./Lag)))
    cvi_allangle = reshape(cv_increment!(xyarr,Lag,nangle,cvi_allangle,diff=diff),mapdim[1]*mapdim[2],nangle) 
    if keepmissing==true
        missing1D               = findall(ismissing,cvi_allangle)
        cvi_allangle[missing1D] .= missing
    end
    @inbounds @views for pix in eachindex(cvi_averaged)
        keepmissing==false && (cvi_averaged[pix] = mean(skipmissing(cvi_allangle[pix,:])))
        keepmissing==true  && (cvi_averaged[pix] = mean(cvi_allangle[pix,:]))
    end
    cvi_averaged = reshape(cvi_averaged,mapdim[1],mapdim[2])
    return(cvi_averaged,cvi_allangle,nangle)
end



"""
    construct_cvimap(cvmap,Lag::Vector{Int64},mapdim; diff="relative",keepmissing=true,BLANK=-1000)

Construct a Centroid Velocity Increment map based on a 'cvmap'. The cvi map is constructed by taking the mean of all rotations of the Centroid Velocity Increment calculation at each pixel. The Lag is the increment. xyarr have to be in 2D (pixel*pixel). Mapdim is the dimension of your 2Dmap. The differences can be absolute or relative. 
"""
function construct_cvimap(cvmap,Lag::Vector{Int64},mapdim; diff="relative",keepmissing=true,BLANK=-1000)
    nangle = floor.(Int,2pi ./(atan.(1 ./Lag)) )
    cvi_allangle_alllag = cv_increment(cvmap,Lag,nangle,diff=diff)
    #println("CV inc done") 
    cvmap = 0.0
    GC.gc()  # CLEANING MEMORY

    cvi_allangle_alllag = Data_preparation.replace_nantomissing(cvi_allangle_alllag)   # Replace NaN into missing for the average
    cvi_allangle_alllag = Data_preparation.replace_blanktomissing(cvi_allangle_alllag,BLANK)   # Replace blank into missing for the average

    cvi_averaged_alllag = Array{Union{Missing,Float64},2}(undef,size(cvi_allangle_alllag)[1],size(Lag)[1])
    cvi_averaged_alllag .= BLANK
    @inbounds @views for lagstep=1:size(Lag)[1]
        if keepmissing==true
            missing1D = findall(ismissing,cvi_allangle_alllag[:,:,lagstep])
            cvi_allangle_alllag[missing1D,lagstep] .= missing
        end
        @inbounds @views for pix=1:size(cvi_allangle_alllag)[1]
            #(ismissing(cvi_allangle_alllag[pix,1,lagstep])) || (cvi_averaged_alllag[pix,lagstep] = mean(skipmissing((cvi_allangle_alllag[pix,:,lagstep]))))
            (cvi_averaged_alllag[pix,lagstep] = mean(skipmissing((cvi_allangle_alllag[pix,:,lagstep]))))
        end
        #cvi_averaged_alllag = Data_preparation.addblank(cvi_averaged_alllag,missing1D,BLANK,(mapdim[1]*mapdim[2],size(Lag)[1]))

    end
    
    cvi_averaged_alllag = reshape(cvi_averaged_alllag,mapdim[1],mapdim[2],size(Lag)[1])
    return(cvi_allangle_alllag,cvi_averaged_alllag,nangle)
end

"""
    construct_cvimap(cvmap,Lag::Int64,mapdim; diff="relative",keepmissing=true) 

Same as construct_cvimap if Lag is a Int64 of one lag. Mapdim is the dimension of your 2Dmap. The differences can be absolute or relative.
"""
function construct_cvimap(cvmap,Lag::Int64,mapdim; diff="relative",keepmissing=true)
    nangle = floor(Int,2pi/(atan(1/Lag)))
    cvi_allangle = reshape(cv_increment(cvmap,Lag,nangle,diff=diff,mapdim),mapdim[1]*mapdim[2],nangle) 
    #ss = open("/tmp/mmap_cvi.bin","w+")
    #cvi_averaged = Mmap.mmap(ss,BitArray,size(cvi_allangle)[1])
    cvi_averaged = Array{Union{Missing,Float64}}(undef,size(cvi_allangle)[1])
    if keepmissing==true
        missing1D               = findall(ismissing,cvi_allangle)
        cvi_allangle[missing1D] .= missing
    end
    @inbounds @views for pix in eachindex(cvi_averaged)
        keepmissing==false && (cvi_averaged[pix] = mean(skipmissing(cvi_allangle[pix,:])))
        keepmissing==true  && (cvi_averaged[pix] = mean(cvi_allangle[pix,:]))
    end
    cvi_averaged = reshape(cvi_averaged,mapdim[1],mapdim[2])
    #return(cvi_averaged,cvi_allangle)
    return(cvi_averaged,nangle)

end








"""
    cv_increment!(xyarr,Lag::Vector{Int64},nangle,cvi_allangle_alllag::Array{Union{Missing, Float64}, 3},cvi_allangle ; diff="relative",periodic=false)

Same as cv_increment(xyarr,Lag::Vector{Int64},nangle; diff="relative",periodic=false) but on preallocated array (array produced outside this function). To produce the array, use : cvi_allangle = Array{Union{Missing,Float64},3}(undef,DataDimension[1],DataDimension[2],nangle) and cvi_allangle_alllag = Array{Union{Missing,Float64},3}(undef,DataDimension[1]*DataDimension[2],nangle,size(Lag)[1])
"""
function cv_increment!(xyarr,Lag::Vector{Int64},nangle,cvi_allangle_alllag::Array{Union{Missing, Float64}, 3},cvi_allangle ; diff="relative",periodic=false)
        @inbounds @views for lagstep=1:size(Lag)[1]
            # Iteration for angles
                 @inbounds @views for angl=1:nangle
                    alpha = angl*2.0*pi/nangle
                    periodic==false && (xyarr_shifted = ShiftedArray(xyarr,(trunc(Int,Lag[lagstep]*cos(alpha)),trunc(Int,Lag[lagstep]*sin(alpha)))))
                    periodic==true  && (xyarr_shifted = circshift(xyarr,(trunc(Int,Lag[lagstep]*cos(alpha)),trunc(Int,Lag[lagstep]*sin(alpha)))))                    
                    # Iteration in columns
                    @inbounds @views for col=1:size(xyarr)[2]
                        # Iteration in rows
                        @inbounds @views for row=1:size(xyarr)[1]
                            cvi_allangle[row,col,angl] = xyarr_shifted[row,col]-xyarr[row,col]
                        end
                    end
                end
                cvi_allangle_alllag[:,:,lagstep] = reshape(cvi_allangle,size(xyarr)[1]*size(xyarr)[2],nangle)
            end
        diff == "absolute" && return(abs.(cvi_allangle_alllag))
        return(cvi_allangle_alllag)
end



"""
    cv_increment!(xyarr,Lag::Int64,blank,cvimap_lag::Array{Union{Missing,Float64},2})

Same as (xyarr,Lag::Int64,nangle; diff="relative",periodic=false) but on preallocated array (array produced outside this function). To produce the array, use : cvi_allangle = Array{Union{Missing,Float64},3}(undef,DataDimension[1],DataDimension[2],nangle).

"""
function cv_increment!(xyarr,Lag::Int64,nangle,cvi_allangle::Array{Union{Missing, Float64}, 3}; diff="relative",periodic=false)
        # Iteration for angles
        @inbounds @views for angl=1:nangle
            alpha = angl*2.0*pi/nangle
            periodic==false && (xyarr_shifted = ShiftedArray(xyarr,(trunc(Int,Lag*cos(alpha)),trunc(Int,Lag*sin(alpha)))))
            periodic==true  && (xyarr_shifted = circshift(xyarr,(trunc(Int,Lag*cos(alpha)),trunc(Int,Lag*sin(alpha)))))
            # Iteration in columns
            @inbounds @views for col=1:size(xyarr)[2]
                # Iteration in rows
                @inbounds @views for row=1:size(xyarr)[1]
                    cvi_allangle[row,col,angl] = xyarr_shifted[row,col]-xyarr[row,col]
                end
            end
        end
        cvi_allangle = reshape(cvi_allangle,size(xyarr)[1]*size(xyarr)[2],nangle)
        diff == "absolute" && return(abs.(cvi_allangle))
        return(cvi_allangle)
end


"""
    cv_increment(xyarr,Lag::Vector{Int64},nangle; diff="relative",periodic=false)

Compute the centroid velocity increment of xyarr at multiple Lag values. Nangle is the number of angle using to compute the differences (it's a value in the parameter file, equal to 192). Diff (default relative) is for differences between two pixels : absolute or relative. Periodic=true (default=false) is for working on periodic data (from simulations like fbm). The returned array will have the first dimension equal to the size of the map (pixel square), the second dimension is the cvi computed at each angle, and the third dimension is for each value of Lag.
"""
function cv_increment(xyarr,Lag::Vector{Int64},nangle; diff="relative",periodic=false, BLANK=-1000) 


    cvi_allangle             = convert(Array{Union{Missing,Float64}},zeros(Float64,size(xyarr)[1],size(xyarr)[2],maximum(nangle)))
    cvi_allangle_alllag      = convert(Array{Union{Missing,Float64}},zeros(Float64,size(xyarr)[1]*size(xyarr)[2],maximum(nangle),size(Lag)[1]))
    cvi_allangle            .= BLANK
    cvi_allangle_alllag     .= BLANK 
    for lagstep=1:size(Lag)[1]
    # Iteration for angles
        for angl=1:nangle[lagstep]
            alpha = angl*2.0*pi/nangle[lagstep]
            periodic==false && (xyarr_shifted = ShiftedArray(xyarr,(trunc(Int,Lag[lagstep]*cos(alpha)),trunc(Int,Lag[lagstep]*sin(alpha)))))
            periodic==true  && (xyarr_shifted = circshift(xyarr,(trunc(Int,Lag[lagstep]*cos(alpha)),trunc(Int,Lag[lagstep]*sin(alpha)))))
            #xyarr_shifted = Data_preparation.replace_missingtoblank(xyarr_shifted,BLANK)
            #xyarr = Data_preparation.replace_missingtoblank(xyarr,BLANK)

            if diff=="absolute"
                # Iteration in columns
                @inbounds @views for col=1:size(xyarr)[2]
                    # Iteration in rows
                    @inbounds @views for row=1:size(xyarr)[1]
                        # First, check if the value is missing. If so, do not continue and go to row+1
                        # Second, if previous step is false, then check if it's a blank value. If not continue and do the difference.
                        
                        ((ismissing(xyarr_shifted[row,col]) || ismissing(xyarr[row,col])) || (xyarr_shifted[row,col] != BLANK || xyarr[row,col] != BLANK)) && (cvi_allangle[row,col,angl] = abs(xyarr_shifted[row,col]-xyarr[row,col]))
                        ((ismissing(xyarr_shifted[row,col]) || ismissing(xyarr[row,col])) || (xyarr_shifted[row,col] == BLANK || xyarr[row,col] == BLANK)) && (cvi_allangle[row,col,angl] = BLANK)
                        #((diff == "relative") && (ismissing(xyarr_shifted[row,col]) || ismissing(xyarr[row,col])) || (xyarr_shifted[row,col] != BLANK || xyarr[row,col] != BLANK)) && (cvi_allangle[row,col,angl] = xyarr_shifted[row,col]-xyarr[row,col])
                        #cvi_allangle[row,col,angl] = xyarr_shifted[row,col]-xyarr[row,col]
                    end
                end
            else
                 for col=1:size(xyarr)[2]
                    # Iteration in rows
                    for row=1:size(xyarr)[1]
                        (ismissing(xyarr_shifted[row,col]) || ismissing(xyarr[row,col]) || xyarr_shifted[row,col] != BLANK || xyarr[row,col] != BLANK) && (cvi_allangle[row,col,angl] = xyarr_shifted[row,col]-xyarr[row,col])
                    end
                end
            end
        end
        cvi_allangle = Data_preparation.replace_blanktomissing(cvi_allangle,BLANK)
        cvi_allangle_alllag[:,:,lagstep] = reshape(cvi_allangle,size(xyarr)[1]*size(xyarr)[2],maximum(nangle))   # (P*P,ROT,Lag)  
    end
    #diff == "absolute" && return(abs.(cvi_allangle_alllag))
    return(cvi_allangle_alllag)
end



"""
    cv_increment(xyarr,Lag::Int64,nangle; diff="relative",periodic=false)

Compute the centroid velocity increment of xyarr at one Lag. Nangle is the number of angle using to compute the differences (it's a value in the parameter file, equal to 192). Diff (default relative) is for differences between two pixels : absolute or relative. Periodic=true (default=false) is for working on periodic data (from simulations like fbm). The returned array will have the first two dimensions equal to the size of the map, and the third dimension is the cvi computed at each angle.
"""
function cv_increment(xyarr,Lag::Int64,nangle,DataDimension; diff="relative",periodic=false)
    cvi_allangle  = convert(Array{Union{Missing,Float64}},zeros(Float64,size(xyarr)[1],size(xyarr)[2],size(nangle)[1]))
    # Iteration for angles
    @inbounds @views for angl=1:nangle[lagstep]
        alpha = angl*2.0*pi/nangle[lagstep]
        periodic==false && (xyarr_shifted = ShiftedArray(xyarr,(trunc(Int,Lag*cos(alpha)),trunc(Int,Lag*sin(alpha)))))
        periodic==true  && (xyarr_shifted = circshift(xyarr,(trunc(Int,Lag*cos(alpha)),trunc(Int,Lag*sin(alpha)))))
        # Iteration in columns
        @inbounds @views for col=1:size(xyarr)[2]
            # Iteration in rows
            @inbounds @views for row=1:size(xyarr)[1]
                cvi_allangle[row,col,angl] = xyarr_shifted[row,col]-xyarr[row,col]
            end
        end
    end
    diff == "absolute" && return(abs.(cvi_allangle))
    return(cvi_allangle)
end



"""
    moment_four(yarr,xarr)

Return the fourth moment order of xarr weighted by yarr.
"""
function moment_four(yarr,xarr)
    return(moment(xarr,4,aweights(yarr)))
end



"""
    moment_one(yarr,xarr)

Compute the first moment order of the array xarr weighted by yarr.
Centroid = int(T*v)dv/int(T)dv
"""
function moment_one(yarr,xarr)
    return(moment(xarr,1,aweights(yarr),0))
    #return(mean(xarr,aweights(yarr)))
end



"""
    moment_three(yarr,xarr)

Return the third moment order of xarr weighted by yarr.
"""
function moment_three(yarr,xarr)
    return(moment(xarr,3,aweights(yarr)))
end



"""
    moment_one_field(arr,velvector)

Return the first moment order of all pixels in an entire field, weighted by velvector. Arr can be a 3D array (ppv) or a 2D array (pv). The function replace all missing values by 
"""
function moment_one_field(arr,SIGMAT,THRESHOLD,velvector,BLANK)
    arr = Data_preparation.blank_inf(arr,SIGMAT*THRESHOLD,0)
    typeof(size(arr)) == Tuple{Int64,Int64,Int64} && (arr = reshape(arr,size(arr)[1]*size(arr)[2],size(arr)[3]))
    eltype(arr)       == Union{Missing,Float64}   && (arr = convert(Array{Float64,2},Data_preparation.replace_missingtoblank(arr,0)))
    #arr = Data_preparation.replace_blanktomissing(arr,BLANK)
    arr         = Data_preparation.permcolrow(arr) # Allow to do the for loop per column (optimized in Julia)
    momentfield = zeros(Float64,size(arr)[2])
    momentfield .= BLANK
    maxvel = maximum(velvector)
    minvel = minimum(velvector)
    for ix=1:size(arr)[2] # For loop on the pixel number (column and row permuted two lines before)
        #println("DAH")
        #println(maximum(arr[:,ix]))
        #(maximum(arr[:,ix])!=BLANK) && (momentfield[ix] = moment_one(arr[:,ix],velvector))
        (momentfield[ix] = moment_one(arr[:,ix],velvector))
        ((momentfield[ix]>maxvel) || (momentfield[ix]<minvel)) && (momentfield[ix]=BLANK)
    end
    return(momentfield)
end



"""
    moment_two(yarr,xarr)

Return the second moment order of xarr weighted by yarr.
"""
function moment_two(yarr,xarr)
    return(moment(xarr,2,aweights(yarr)))
end



"""
    moment_two_field(arr,velvector)

Return the order 2 moment of an entire field.
"""
function moment_two_field(arr,velvector)
    typeof(size(arr)) == Tuple{Int64,Int64,Int64} && (arr = reshape(arr,size(arr)[1]*size(arr)[2],size(arr)[3]))
    eltype(arr)       == Union{Missing,Float64}   && (arr = convert(Array{Float64,2},Data_preparation.replace_missingtoblank(arr,NaN)))
    arr         = Data_preparation.permcolrow(arr)
    momentfield = zeros(Float64,size(arr)[1])
    for ix=1:size(momentfield)[1]
        momentfield[ix] = moment_two(arr[:,ix],velvector)
    end
    return(momentfield)
end




"""
    multiple_moment(cube::Array{Float64,3},velvector)

Return the first or second velocity moment order of a 3D array. The 3D array in input have its first dimension giving the pixel position, the second dimension giving spectra, and the third dimension giving other data with the same dimensions. This function was built in order to obtain in the third dimensions multiple centroid velocity maps with data reconstructed from different number of principal component (see Principal Component Analysis)
"""
function multiple_moment(cube::Array{Float64,3},velvector)
    momentone = zeros(Float64,size(cube)[1],size(cube)[3])
    momenttwo = similar(momentone)
    for ix = 1:size(cube)[3]
        for pix = 1:size(cube)[1]
            momentone[pix,ix] = moment_one(cube[pix,:,ix],velvector)
            momenttwo[pix,ix] = moment_two(cube[pix,:,ix],velvector)
        end
    end
    return(momentone,momenttwo)
end




"""
    delete_aberations(arr1,arr2)

Replace values in arr1 which are higher or lower than values in arr2 by missing. Don't use an arr1 with missing values !
"""
function delete_aberations(arr1,arr2)
    min = minimum(arr2)
    max = maximum(arr2)
    for jx in eachindex(arr1)
        arr1[jx]<=min && (arr1[jx]=-100000)
        arr1[jx]>=max && (arr1[jx]=-100000)
    end
    arr1 = replace_blanktomissing(arr1,-100000)
    return(arr1)

end






end #module
