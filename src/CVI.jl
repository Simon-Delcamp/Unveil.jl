module CVI
include("Dataprep.jl") #Read and write fits
using .Dataprep

using ShiftedArrays, StatsBase
using ProgressBars

export construct_cvimap!
export construct_cvimap
export cv_increment!
export cv_increment
export moment_one_field 





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
   # nangle = fill(floor.(Int,2pi ./atan(1)),size(Lag)[1]) #Array{Float64}(undef,size(Lag)[1]).*0 .+floor.(Int,2pi ./atan(1))    #floor.(Int,2pi ./(atan.(1 ./Lag)) )
    nangle = floor.(Int,2pi ./(atan.(1 ./Lag)) )
    cvi_allangle_alllag = cv_increment(cvmap,Lag,nangle,diff=diff)
    #println("CV inc done") 
    cvmap = 0.0
    GC.gc()  # CLEANING MEMORY
    cvi_allangle_alllag = Dataprep.replace_nantomissing(cvi_allangle_alllag)   # Replace NaN into missing for the average
    cvi_allangle_alllag = Dataprep.replace_blanktomissing(cvi_allangle_alllag,BLANK)   # Replace blank into missing for the average

    #cvi_averaged_alllag = Array{Union{Missing,Float64},2}(undef,size(cvi_allangle_alllag)[1],size(Lag)[1])
    #cvi_averaged_alllag .= BLANK
    cvi_averaged_alllag = Array{Union{Missing,Float64},3}(undef,mapdim[1],mapdim[2],size(Lag)[1])
    cvi_averaged_alllag .= BLANK

    @inbounds @views for lagstep=1:size(Lag)[1]
        if keepmissing==true
            missing1D = findall(ismissing,cvi_allangle_alllag[:,:,lagstep])
            cvi_allangle_alllag[missing1D,lagstep] .= missing
        end
        @inbounds @views for col=1:size(cvi_allangle_alllag)[2]
            # Iteration in rows
            @inbounds @views for row=1:size(cvi_allangle_alllag)[1]
            #(ismissing(cvi_allangle_alllag[pix,1,lagstep])) || (cvi_averaged_alllag[pix,lagstep] = mean(skipmissing((cvi_allangle_alllag[pix,:,lagstep]))))
                (cvi_averaged_alllag[row,col,lagstep] = mean(skipmissing((cvi_allangle_alllag[row,col,:,lagstep]))))
            end
        end
        #@inbounds @views for pix=1:size(cvi_allangle_alllag)[1]
        #    #(ismissing(cvi_allangle_alllag[pix,1,lagstep])) || (cvi_averaged_alllag[pix,lagstep] = mean(skipmissing((cvi_allangle_alllag[pix,:,lagstep]))))
        #    (cvi_averaged_alllag[pix,lagstep] = mean(skipmissing((cvi_allangle_alllag[pix,:,lagstep]))))
        #end
        #cvi_averaged_alllag = Dataprep.addblank(cvi_averaged_alllag,missing1D,BLANK,(mapdim[1]*mapdim[2],size(Lag)[1]))

    end
    
    #cvi_averaged_alllag = reshape(cvi_averaged_alllag,mapdim[1],mapdim[2],size(Lag)[1])
    return(cvi_allangle_alllag,cvi_averaged_alllag,nangle)
end



"""
    construct_cvimap(cvmap,Lag::Int64,mapdim; diff="relative",keepmissing=true) 

Same as construct_cvimap if Lag is a Int64 of one lag. Mapdim is the dimension of your 2Dmap. The differences can be absolute or relative.
"""
function construct_cvimap(cvmap,Lag::Int64,mapdim; diff="relative",keepmissing=true,BLANK=-1000)
    nangle = floor(Int,2pi/(atan(1/Lag)))
    cvi_allangle_alllag = cv_increment(cvmap,Lag,nangle,mapdim,diff=diff)

    cvmap = 0.0
    GC.gc()  # CLEANING MEMORY
    cvi_allangle_alllag = Dataprep.replace_nantomissing(cvi_allangle_alllag)   # Replace NaN into missing for the average
    cvi_allangle_alllag = Dataprep.replace_blanktomissing(cvi_allangle_alllag,BLANK)   # Replace blank into missing for the average

    #cvi_averaged_alllag = Array{Union{Missing,Float64},2}(undef,size(cvi_allangle_alllag)[1],size(Lag)[1])
    #cvi_averaged_alllag .= BLANK
    cvi_averaged_alllag = Array{Union{Missing,Float64},2}(undef,mapdim[1],mapdim[2])
    cvi_averaged_alllag .= BLANK

    #println(size(cvi_allangle_alllag))
    #println(size(cvi_averaged_alllag))

    if keepmissing==true
        missing1D = findall(ismissing,cvi_allangle_alllag)
        cvi_allangle_alllag[missing1D] .= missing
    end
    #for lagstep in ProgressBar(1:size(Lag)[1])
    #@inbounds @views for col=1:size(cvi_allangle_alllag)[2]
    for col in ProgressBar(1:size(cvi_allangle_alllag)[2])
        # Iteration in rows
        @inbounds @views for row=1:size(cvi_allangle_alllag)[1]
                (cvi_averaged_alllag[row,col] = mean(skipmissing((cvi_allangle_alllag[row,col,:]))))
        end
    end


    return(cvi_allangle_alllag,cvi_averaged_alllag,nangle)
    
    ####cvi_allangle = reshape(cv_increment(cvmap,Lag,nangle,diff=diff,mapdim),mapdim[1]*mapdim[2],nangle) 
    #####ss = open("/tmp/mmap_cvi.bin","w+")
    #####cvi_averaged = Mmap.mmap(ss,BitArray,size(cvi_allangle)[1])
    ####cvi_averaged = Array{Union{Missing,Float64}}(undef,size(cvi_allangle)[1])
    ####if keepmissing==true
    ####    missing1D               = findall(ismissing,cvi_allangle)
    ####    cvi_allangle[missing1D] .= missing
    ####end
    ####@inbounds @views for pix in eachindex(cvi_averaged)
    ####    keepmissing==false && (cvi_averaged[pix] = mean(skipmissing(cvi_allangle[pix,:])))
    ####    keepmissing==true  && (cvi_averaged[pix] = mean(cvi_allangle[pix,:]))
    ####end
    ####cvi_averaged = reshape(cvi_averaged,mapdim[1],mapdim[2])
    #####return(cvi_averaged,cvi_allangle)
    ####return(cvi_allangle,cvi_averaged,nangle)

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
    #cvi_allangle             = convert(Array{Union{Missing,Float64}},zeros(Float64,size(xyarr)[1],size(xyarr)[2],maximum(nangle)))
    cvi_allangle_alllag      = convert(Array{Union{Missing,Float64}},zeros(Float64,size(xyarr)[1],size(xyarr)[2],maximum(nangle),size(Lag)[1]))
    #cvi_allangle            .= BLANK
    cvi_allangle_alllag     .= BLANK 
    for lagstep in ProgressBar(1:size(Lag)[1]) #WORKED # for lagstep=1:size(Lag)[1]
    # Iteration for angles
        for angl=1:nangle[lagstep]
            alpha = angl*2.0*pi/nangle[lagstep]
            periodic==false && (xyarr_shifted = ShiftedArray(xyarr,(trunc(Int,Lag[lagstep]*cos(alpha)),trunc(Int,Lag[lagstep]*sin(alpha)))))
            periodic==true  && (xyarr_shifted = circshift(xyarr,(trunc(Int,Lag[lagstep]*cos(alpha)),trunc(Int,Lag[lagstep]*sin(alpha)))))
            #xyarr_shifted = Dataprep.replace_missingtoblank(xyarr_shifted,BLANK)
            #xyarr = Dataprep.replace_missingtoblank(xyarr,BLANK)

            if diff=="absolute"
                # Iteration in columns
                @inbounds @views for col=1:size(xyarr)[2]
                    # Iteration in rows
                    @inbounds @views for row=1:size(xyarr)[1]
                        # First, check if the value is missing. If so, do not continue and go to row+1
                        # Second, if previous step is false, then check if it's a blank value. If not continue and do the difference.
                        
                        ((ismissing(xyarr_shifted[row,col]) || ismissing(xyarr[row,col])) || (xyarr_shifted[row,col] != BLANK || xyarr[row,col] != BLANK)) && (cvi_allangle_alllag[row,col,angl] = abs(xyarr_shifted[row,col]-xyarr[row,col]))
                        ((ismissing(xyarr_shifted[row,col]) || ismissing(xyarr[row,col])) || (xyarr_shifted[row,col] == BLANK || xyarr[row,col] == BLANK)) && (cvi_allangle_alllag[row,col,angl] = BLANK)
                        #((diff == "relative") && (ismissing(xyarr_shifted[row,col]) || ismissing(xyarr[row,col])) || (xyarr_shifted[row,col] != BLANK || xyarr[row,col] != BLANK)) && (cvi_allangle[row,col,angl] = xyarr_shifted[row,col]-xyarr[row,col])
                        #cvi_allangle[row,col,angl] = xyarr_shifted[row,col]-xyarr[row,col]
                    end
                end
            else
                @inbounds @views for col=1:size(xyarr)[2]
                    # Iteration in rows
                    @inbounds @views for row=1:size(xyarr)[1]
                        ((ismissing(xyarr_shifted[row,col])) || ismissing(xyarr[row,col]) || (xyarr_shifted[row,col] != BLANK || xyarr[row,col] != BLANK)) && (cvi_allangle_alllag[row,col,angl,lagstep] = xyarr_shifted[row,col]-xyarr[row,col])
                    end
                end
            end
        end
        #cvi_allangle = Dataprep.replace_blanktomissing(cvi_allangle,BLANK)
        #cvi_allangle_alllag[:,:,lagstep] = reshape(cvi_allangle,size(xyarr)[1]*size(xyarr)[2],maximum(nangle))   # (P*P,ROT,Lag)  
    end
    #diff == "absolute" && return(abs.(cvi_allangle_alllag))
    return(cvi_allangle_alllag)
end



"""
    cv_increment(xyarr,Lag::Int64,nangle; diff="relative",periodic=false)

Compute the centroid velocity increment of xyarr at one Lag. Nangle is the number of angle using to compute the differences (it's a value in the parameter file, equal to 192). Diff (default relative) is for differences between two pixels : absolute or relative. Periodic=true (default=false) is for working on periodic data (from simulations like fbm). The returned array will have the first two dimensions equal to the size of the map, and the third dimension is the cvi computed at each angle.
"""
function cv_increment(xyarr,Lag::Int64,nangle,DataDimension; diff="relative",periodic=false,BLANK=-1000)
    cvi_allangle_alllag      = convert(Array{Union{Missing,Float64}},zeros(Float64,size(xyarr)[1],size(xyarr)[2],nangle))
    cvi_allangle_alllag     .= BLANK 
    for angl=1:nangle
    
        alpha = angl*2.0*pi/nangle
        periodic==false && (xyarr_shifted = ShiftedArray(xyarr,(trunc(Int,Lag*cos(alpha)),trunc(Int,Lag*sin(alpha)))))
        periodic==true  && (xyarr_shifted = circshift(xyarr,(trunc(Int,Lag*cos(alpha)),trunc(Int,Lag*sin(alpha)))))
        #xyarr_shifted = Dataprep.replace_missingtoblank(xyarr_shifted,BLANK)
        #xyarr = Dataprep.replace_missingtoblank(xyarr,BLANK)

        if diff=="absolute"
            # Iteration in columns
        @inbounds @views for col=1:size(xyarr)[2]
            # Iteration in rows
                @inbounds @views for row=1:size(xyarr)[1]
                    # First, check if the value is missing. If so, do not continue and go to row+1
                    # Second, if previous step is false, then check if it's a blank value. If not continue and do the difference.
                    
                    ((ismissing(xyarr_shifted[row,col]) || ismissing(xyarr[row,col])) || (xyarr_shifted[row,col] != BLANK || xyarr[row,col] != BLANK)) && (cvi_allangle_alllag[row,col,angl] = abs(xyarr_shifted[row,col]-xyarr[row,col]))
                    ((ismissing(xyarr_shifted[row,col]) || ismissing(xyarr[row,col])) || (xyarr_shifted[row,col] == BLANK || xyarr[row,col] == BLANK)) && (cvi_allangle_alllag[row,col,angl] = BLANK)
                    #((diff == "relative") && (ismissing(xyarr_shifted[row,col]) || ismissing(xyarr[row,col])) || (xyarr_shifted[row,col] != BLANK || xyarr[row,col] != BLANK)) && (cvi_allangle[row,col,angl] = xyarr_shifted[row,col]-xyarr[row,col])
                    #cvi_allangle[row,col,angl] = xyarr_shifted[row,col]-xyarr[row,col]
                end
            end
        else
            @inbounds @views for col=1:size(xyarr)[2]
                # Iteration in rows
                @inbounds @views for row=1:size(xyarr)[1]
                    ((ismissing(xyarr_shifted[row,col])) || ismissing(xyarr[row,col]) || (xyarr_shifted[row,col] != BLANK || xyarr[row,col] != BLANK)) && (cvi_allangle_alllag[row,col,angl] = xyarr_shifted[row,col]-xyarr[row,col])
                end
            end
        end
    end
        #cvi_allangle = Dataprep.replace_blanktomissing(cvi_allangle,BLANK)
        #cvi_allangle_alllag[:,:,lagstep] = reshape(cvi_allangle,size(xyarr)[1]*size(xyarr)[2],maximum(nangle))   # (P*P,ROT,Lag)  
    
    #diff == "absolute" && return(abs.(cvi_allangle_alllag))
    return(cvi_allangle_alllag)

    ############"cvi_allangle  = convert(Array{Union{Missing,Float64}},zeros(Float64,size(xyarr)[1],size(xyarr)[2],nangle))
    ############"# Iteration for angles
    ############"@inbounds @views for angl=1:nangle#[lagstep]
    ############"    alpha = angl*2.0*pi/nangle#[lagstep]
    ############"    periodic==false && (xyarr_shifted = ShiftedArray(xyarr,(trunc(Int,Lag*cos(alpha)),trunc(Int,Lag*sin(alpha)))))
    ############"    periodic==true  && (xyarr_shifted = circshift(xyarr,(trunc(Int,Lag*cos(alpha)),trunc(Int,Lag*sin(alpha)))))
    ############"    # Iteration in columns
    ############"    @inbounds @views for col=1:size(xyarr)[2]
    ############"        # Iteration in rows
    ############"        @inbounds @views for row=1:size(xyarr)[1]
    ############"            cvi_allangle[row,col,angl] = xyarr_shifted[row,col]-xyarr[row,col]
    ############"        end
    ############"    end
    ############"end
    ############"diff == "absolute" && return(abs.(cvi_allangle))
    ############"return(cvi_allangle)
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
    moment_one_field(arr,velvector)

Return the first moment order of all pixels in an entire field, weighted by velvector. Arr can be a 3D array (ppv) or a 2D array (pv). The function replace all missing values by 
"""
function moment_one_field(arr,SIGMAT,THRESHOLD,velvector,BLANK)
    arr = Dataprep.blank_inf(arr,SIGMAT*THRESHOLD,0)
    typeof(size(arr)) == Tuple{Int64,Int64,Int64} && (arr = reshape(arr,size(arr)[1]*size(arr)[2],size(arr)[3]))
    #eltype(arr)       == Union{Missing,Float64}   && (arr = convert(Array{Float64,2},Dataprep.replace_missingtoblank(arr,0)))
    eltype(arr)       == Union{Missing,Float64}   && (arr = convert(Array{Float64,2},Dataprep.replace_missingtoblank(arr,BLANK)))
    #arr = Dataprep.replace_blanktomissing(arr,BLANK)
    arr         = Dataprep.permcolrow(arr) # Allow to do the for loop per column (optimized in Julia)
    momentfield = zeros(Float64,size(arr)[2])
    momentfield .= BLANK
    maxvel = maximum(velvector)
    minvel = minimum(velvector)
    for ix=1:size(arr)[2] # For loop on the pixel number (column and row permuted two lines before)
        #println("DAH")
        #println(maximum(arr[:,ix]))
        #(maximum(arr[:,ix])!=BLANK) && (momentfield[ix] = moment_one(arr[:,ix],velvector))
        (momentfield[ix] = moment_one(arr[:,ix],velvector))
        #(momentfield[ix] = moment_one(skipmissing(arr[:,ix]),velvector))  # MOMENTS can't works with skipmissing
        ((momentfield[ix]>maxvel) || (momentfield[ix]<minvel)) && (momentfield[ix]=BLANK)
    end
    return(momentfield)
end




end #module