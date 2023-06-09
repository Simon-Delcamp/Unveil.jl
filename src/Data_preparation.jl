###################################################################
# Functions used to work and prepare data. For example, writing fits, remove blanking values, change blanking values to missing one, reading ppv...
# Called these functions by calling (or writting on a script ):
#       >include("../src/Data_preparation.jl")
#       >using .Data_preparation
#       >output = Data_preparation.NameOfTheFunction(input)
###################################################################

module Data_preparation

export addblank
export addmask
export allmissing
export blank_equal
export blank_inf
export blank_inf_abs
export blank_snr
export blank_sup
export boolmatrix_missing
export clear_var
export coordtopix
export data_path
export deletemissing
export delete_allnotvalue
export directory_prep
export header_entries
export lookingfor_interval
export noise_canals
export pca_prep
export permcolrow
export pixtocoord
export read_dim
export read_fits_pp
export read_fits_pv
export read_fits_ppv
export read_fits_ppf
export read_fits_spect
export read_var_files
export regular_blanking
export replace_blanktomissing
export replace_inf_in_nan
export replace_missingtoblank
export replace_nantoblank
export replace_nantomissing
export replace_nosignal
export valid_header
export wait_for_key
export write_cv_fits
export write_cvi_fits
export write_dat
export write_fits
export write_pca_fits


using FITSIO, MultivariateStats, StatsBase, Distributions, DelimitedFiles
using Makie, GLMakie, StaticArrays, PolygonOps # For add mask funtcion

import GLMakie.heatmap!# For add mask funtcion
import GLMakie.heatmap # For add mask funtcion
import GLMakie.Scene
import GLMakie.cgrad
import GLMakie.scatter!
import GLMakie.events
import GLMakie.Mouse
import GLMakie.mouseposition



"""
    addblank(data,missing_arr,blank,data_dimension)

Reconstruct a data with missing values after they where deleted. Usually deleted in order to compute a PCA. Can reconstruct 3D cubes (e.g. PPV) as well as 2D maps (PP). In the first case, the input missing_arr should be in 2D, while in the second case it should be in 1D. The input missing_arr is produces with the function Data_preparation.pca_prep. The input data should be a 1D vector. Data_dimension has to be equal to the dimension you want to obtain. 
    To reproduce a 3D cube, it should have blank values at the exact same positions across the 3rd dimension. 
"""
function addblank(data,missing_arr,blank,data_dimension)
    typeof(data_dimension)==Tuple{Int64,Int64,Int64} && (datawithmissing = zeros(Float64,data_dimension[1],data_dimension[2],data_dimension[3]))
    typeof(data_dimension)==Tuple{Int64,Int64}       && (datawithmissing = zeros(Float64,data_dimension[1],data_dimension[2]))
    count = 1 # for iteration on the reconstructed data without missing values
    for ix in eachindex(missing_arr)
        missing_arr[ix]==0 && (datawithmissing[ix]=data[count])
        missing_arr[ix]==0 && (count +=1)
        missing_arr[ix]==1 && (datawithmissing[ix]=blank)
    end

    return(datawithmissing)
end


"""

    addmask(cube,colorscale,nbmask,data_dimension,delta_xvec,delta_yvec)

Blank all the pixels in the third dimension selected manually by the user (using a contour). Can do multiple mask.
"""
function addmask(cube,colorscale,nbmask,data_dimension,delta_xvec,delta_yvec,BLANK,PATHTOSAVE ; alreadydone=false)
    cube = convert(Array{Float64,3},replace_missingtoblank(cube,BLANK))
    cube = permcolrow(cube)
    pixelpos = Array{Float64}(undef,data_dimension[1]*data_dimension[2],data_dimension[3])
    ins = ones(data_dimension[1]*data_dimension[2])
    if alreadydone==false
        println(" ")
        println(" This process can be long ; wait for the next instruction.")
        scene = Scene()
        for kx=1:nbmask
            println(" Mask number $(kx) ")
            heatmap!(scene,delta_xvec,delta_yvec,cube[:,:,100],scale_plot=false,colormap=cgrad(:curl),colorrange=colorscale)
            display(scene)
            clicks = Node(Point2f0[(0,0)])
            on(events(scene).mousebutton,priority=0) do event
                if event.button == Mouse.left
                    if event.action == Mouse.press
                        pos = mouseposition(scene)
                        clicks[] = push!(clicks[],pos)
                    end
                end
            end
            scatter!(scene, clicks, color = :red, marker = '+', markersize = 4)
            println("You can circle the pixels where mask have to be used. When done, please write the name of the file you want this mask to be saved, then press any key.")
            nametosave = readline()
            xpoly = []
            ypoly = []
            for ix=2:size(clicks.val)[1]
                push!(xpoly,-clicks.val[ix][1])
                push!(ypoly,clicks.val[ix][2])
            end
            push!(xpoly,xpoly[1])
            push!(ypoly,ypoly[1])
            polygon = SVector.(xpoly,ypoly)
            xa = delta_xvec
            ya = delta_yvec
            points = vec(SVector.(xa',ya))
            inside = [inpolygon(p, polygon; in=0., on=0., out=1.) for p in points]
            ins .*= inside
            open("$(PATHTOSAVE)/Data/$(DIRECTORYNAME)/$nametosave.txt","w") do io
                writedlm(io, ins)
            end
        end
        for jx =1:data_dimension[3]
            img = permcolrow(cube[:,:,jx])
            for ix in eachindex(ins)
                pixelpos[ix,jx] = ins[ix]*img[ix]
            end
        end
        pixelpos = replace_blanktomissing(pixelpos,0.)
        pixelpos = replace_blanktomissing(pixelpos,BLANK)
        pixelpos = replace_blanktomissing(pixelpos,-0.)
        return(reshape(pixelpos,data_dimension[1],data_dimension[2],data_dimension[3]))
    elseif alreadydone==true
        for masknb = 1:nbmask
            println(" What is the name of the txt file (without the extension) ?")
            filename = readline()
            readed = readdlm("sims/","$filename.txt")
            ins = ins.*readed
        end
        for jx =1:data_dimension[3]
            img = permcolrow(cube[:,:,jx])
            for ix in eachindex(ins)
                pixelpos[ix,jx] = ins[ix]*img[ix]
            end
        end
        pixelpos = replace_blanktomissing(pixelpos,0.)
        pixelpos = replace_blanktomissing(pixelpos,-1000.0)
        pixelpos = replace_blanktomissing(pixelpos,-0.)
        return(reshape(pixelpos,data_dimension[1],data_dimension[2],data_dimension[3]))
    end
end



"""
    allmissing(data)

Return true if the data contain any missing value.
"""
function allmissing(data)
        return(all(ismissing,data))
end



"""
    blank_equal(array,oldvalue,newvalue)

Change the values equals to a value by another one in an array. This is useful if change of blanking is desired. You have to use the function "replace_nantoblank" if the old value is NaN.
"""
function blank_equal(array,oldvalue,newvalue)
    new_array=array.*1
    for i in eachindex(array)
        (array[i]==oldvalue) && (new_array[i]=newvalue)
    end
    return new_array
end



"""
    blank_inf(array,oldvalue,newvalue)

Change the values tinier than an old value to a new value in an array. This is useful if a change of blanking is desired
"""
function blank_inf(array,oldvalue,newvalue)
    new_array=array.*1
    for i in eachindex(array)
        array[i]<=oldvalue && (new_array[i]=newvalue)
    end
    return new_array
end



"""
    blank_inf_abs(array,oldvalue,newvalue)

Change the values tinier than an old value (its absolute value) to a new value in an array. This is useful if a change of blanking is desired
"""
function blank_inf_abs(array,oldvalue,newvalue)
    new_array=array.*1
    for i in eachindex(array)
        abs(array[i])<=oldvalue && (new_array[i]=newvalue)
    end
    return new_array
end



"""
    blank_snr(data2D,NoiseCan,thresh)

Blank all the spectra which have a signal to noise ratio inferior than a threshold. Work on 2D data, the first dimension giving the pixel position and the second dimension the spectra referring to this position.
"""
function blank_snr(data2D,NoiseCan,thresh)
    data2D = replace_missingtoblank(data2D,-10000)
    snr_field = similar(data2D[:,1])    # Will be multiplied by the data2D at the end.
    for ix=1:size(snr_field)[1]
        snr_field[ix] = maximum(data2D[ix,:])/(std(data2D[ix,NoiseCan[1]])+std(data2D[ix,NoiseCan[2]]))/2
    end
    snr_field = blank_inf(snr_field,thresh,-10000)  # Replace the bad pixel by -100000, which will be replace by missing after.
    snr_field = blank_sup(snr_field,thresh,1.0)     # Replace the good pixels by 1 (will be multiplied by the data2D)
    snr_field = blank_equal(snr_field,thresh,1.)
    snr_field = replace_blanktomissing(snr_field,-10000)
    snr_field = replace_nantomissing(snr_field)
    #(snr_field.<thresh) && replace(snr_field,(snr_field.<thresh)=>missing)
    #replace(snr_field,(snr_field.>thresh)=>1.)

    return(data2D.*snr_field)
end



"""
    blank_sup(array,oldvalue,newvalue)

Change the values bigger than an old value to a new value of an array . This is useful if a change of blanking is desired
"""
function blank_sup(array,oldvalue,newvalue)
    for i in eachindex(array)
        array[i]>oldvalue && (array[i]=newvalue)
    end
    return array
end



"""
    boolmatrix_missing(array)

Return a 1D and a 2D boolean matrices indicating the indexes of missing values in array. Array should be a 2 dimension data (pixels*spectra for example)
"""
function boolmatrix_missing(array)
        missing2D = broadcast(isequal,array,missing)
        missing1D = findall(ismissing,array[:])
        return(missing1D, missing2D)
end



"""
    clear_var(var1,var2,var3,var4)

Replace all variable in input by 0 to "free" space.
"""
function clear_var(var1,var2,var3,var4)
    var1 = 0
    var2 = 0
    var3 = 0
    var4 = 0
    return(var1, var2, var3, var4)
end



"""
    coordtopix(header,coord,unit,ax ; off=true)

Do the conversion between physical coordinate and pixel number thanks to the header of the data. Works only for small-size fields. Need the unit of the coordinate given (arcsec,arcmin,rad,deg only) and the axis. Consider by default that the coordinate given is an offset (place at the center of the image). The values of CRVAL1... are considered in deg.
"""
function coordtopix(header,coord,unit,ax,xvec ; off=true)
    unit == "arcsec" && (convfactor = 1/3600)
    unit == "arcmin" && (convfactor = 1/60)
    unit == "deg"    && (convfactor = 1)
    unit == "rad"    && (convfactor = 180/pi)
    off  == true     && (coord += xvec[trunc(Int,size(xvec)[1]/2)])
    ax   == "x"      && (pixpos=trunc(Int,(coord*convfactor-header["CRVAL1"])/header["CDELT1"]+header["CRPIX1"]))
    ax   == "y"      && (pixpos=trunc(Int,(coord*convfactor-header["CRVAL2"])/header["CDELT2"]+header["CRPIX2"]))

    return(pixpos)
end



"""
    deletemissing(data,missing1D)

Delete all rows where there is at least one missing value using a boolean vector of the missing values.
"""
function deletemissing(data,missing1D)
    DATASIZETYPE = eltype(size(data))
    if typeof(size(data))==Tuple{DATASIZETYPE,DATASIZETYPE}
        return(data[setdiff(1:end,missing1D),:])
    else
        error("Please send a 2D map (pixel*pixel,velocity)")
    end
end


"""
    delete_allnotvalue(data,blank)

Delete from data all value which are not a valid value, e.g. missing, NaN and blank.
"""
function delete_allnotvalue(data,blank)
    data = replace_nantomissing(data)
    data = replace_blanktomissing(data,blank)
    missing1D, missing2D = boolmatrix_missing(data)
    data = convert(Array{Float64},deleteat!(data[:],missing1D))
    return(data)
end



"""
    directory_prep(PATHTOSAVE::String)

Construct a new directory to save plots of the convergence criteria (Convergence_Criteria) process and another to save the fits (Data).
"""
function directory_prep(PATHTOSAVE)
    #(isdir("$(PATHTOSAVE)/Plots/"))==0                                            && mkdir("$(PATHTOSAVE)/Plots/")
    #(isdir("$(PATHTOSAVE)/Plots/$(DIRECTORYNAME)"))==0                            && mkdir("$(PATHTOSAVE)/Plots/$(DIRECTORYNAME)")
    #(isdir("$(PATHTOSAVE)/Plots/$(DIRECTORYNAME)/CVI"))==0                        && mkdir("$(PATHTOSAVE)/Plots/$(DIRECTORYNAME)/CVI")
    #(isdir("$(PATHTOSAVE)/Plots/$(DIRECTORYNAME)/CV"))==0                         && mkdir(("$(PATHTOSAVE)/Plots/$(DIRECTORYNAME)/CV"))
    #(isdir("$(PATHTOSAVE)/Plots/$(DIRECTORYNAME)/convergence_criteria"))==0       && mkdir("$(PATHTOSAVE)/Plots/$(DIRECTORYNAME)/convergence_criteria")

    (isdir("$(PATHTOSAVE)/Data/"))==0                                            && mkdir("$(PATHTOSAVE)/Data/")
    (isdir("$(PATHTOSAVE)/Figures/"))==0                                            && mkdir("$(PATHTOSAVE)/Figures/")


    #(isdir("$(PATHTOSAVE)/Data/$(DIRECTORYNAME)"))==0                            && mkdir("$(PATHTOSAVE)/Data/$(DIRECTORYNAME)")
    #(isdir("$(PATHTOSAVE)/Data/$(DIRECTORYNAME)/CVI"))==0                        && mkdir("$(PATHTOSAVE)/Data/$(DIRECTORYNAME)/CVI")
    #(isdir("$(PATHTOSAVE)/Data/$(DIRECTORYNAME)/CV"))==0                         && mkdir(("$(PATHTOSAVE)/Data/$(DIRECTORYNAME)/CV"))
    #(isdir("$(PATHTOSAVE)/Data/$(DIRECTORYNAME)/convergence_criteria"))==0       && mkdir("$(PATHTOSAVE)/Data/$(DIRECTORYNAME)/convergence_criteria")
    #(isdir("$(PATHTOSAVE)/Convergence_Criteria"))==0       && mkdir("$(PATHTOSAVE)/Convergence_Criteria")
    return()
end

"""
    flatiterator(iterator)

Transform an iterator to a matrix
"""

function flatiterator(iterator)
    n = length(first(iterator))
    return(reshape(collect(Iterators.flatten(iterator)), :, n))
end#flatiterator

"""
    lookingfor_interval(arr,value,sigma)

Return the array (1D) cut at [value-sigma:value+sigma]. Return also the
indices of these values.
"""
function lookingfor_interval(arr,value,sigma)
    deltav = arr[2]-arr[1]
    if deltav<0
        # xmin and xmax are the limits of the interval
        xmin = findall(x->x>=(value+sigma),arr[:])
        xmax = findall(x->x<=(value-sigma),arr[:])
        return(arr[xmax[1]:xmin[end]],xmax[1]:xmin[end])
    elseif deltav>0
        xmin = findall(x->x>=(value-sigma),arr[:])
        xmax = findall(x->x<=(value+sigma),arr[:])
        return arr[xmin[1]:xmax[end]],xmin[1]:xmax[end]
    end
end


"""
    noise_channels()

Construct a vector with two subranges. 
"""
function noise_channels()
    println("")
    println("What are the canals with only noise ?")
    println("How many segments ? 1 or 2 ?")
    segment         = parse(Int64,readline())
    println("FIRST SEGMENT : First canal with only noise :")
    firstcanalnoise = parse(Int64,readline())
    println("FIRST SEGMENT : Last canal with only noise (without the emission):")
    lastcanalnoise  = parse(Int64,readline())
    noisecanals1=firstcanalnoise:lastcanalnoise
    if segment==2 
        println("SECOND SEGMENT : First canal with only noise :")
        firstcanalnoise = parse(Int64,readline())
        println("SECOND SEGMENT : Last canal with only noise (without the emission):")
        lastcanalnoise  = parse(Int64,readline())
        noisecanals2=firstcanalnoise:lastcanalnoise
        return(noisecanals1,noisecanals2)
    elseif segment==1
        return(noisecanals1)
    else
        error("Not a valid number of segment (1 or 2)")
    end
end



"""
    pca_prep(arr,arraydimension)

Prepare data (given in 2D (pv) or 3D (ppv)) in order to conduct PCA on them by deleting missing values and reshape in 2D.

Will also check if any missing value still exist in the data, showing that the data are not regularly blanked. Return data in 2 dimensions without missing value, a 1D boolean matrix to check if any missing value still exist in the data (if 0 element then no missing value) and a 2D matrix with booleans corresponding to missing values in the dataset (used to reconstruct the data with the missing values).
"""
function pca_prep(arr,arraydimension)
        DATASIZETYPE = eltype(size(arr))
        (typeof(size(arr))==Tuple{DATASIZETYPE,DATASIZETYPE,DATASIZETYPE}) && (arr = reshape(arr,arraydimension[1]*arraydimension[2],arraydimension[3])) #3D -> 2D

        #println("produce boolean matrix")
        missing1D, missing2D = Data_preparation.boolmatrix_missing(arr) #Boolean matrix to catch missing values
        #println("convert and delete")
        arr2D_nomissing     = Data_preparation.deletemissing(arr,missing1D) #Data without missing value
        arr = 0.0
        #println("regular blanking")
        regular_blanking(arr2D_nomissing)
        return(arr2D_nomissing,missing1D,missing2D)
end



"""
    permcolrow(arr)

Return the array with the first and second dimensions permuted. 2D and 3D arrays accepted.
"""
function permcolrow(arr)
    typeof(size(arr))==Tuple{Int64,Int64}             && return(permutedims(arr,(2,1)))
    typeof(size(arr))==Tuple{Int64,Int64,Int64}       && return(permutedims(arr,(2,1,3)))
    typeof(size(arr))==Tuple{Int64,Int64,Int64,Int64} && error("Not a valid dimension ; Need 3D or 2D.")

end



"""
    pixtocoord(header,unit)

Do the conversion between pixel number and physical coordinate thanks to the header of the data. Return the X and Y coordinates in the physical coordinates you want, then the X and Y coordinate centered at the middle of the plot.
"""
function pixtocoord(header,unit ; factor=1)
    unit == "arcsec" && (convfactor=3600)
    unit == "arcmin" && (convfactor=60)
    unit == "other"    && (convfactor=factor)#80/pi*60)

    xvec       = Array(header["CRVAL1"]+(1-header["CRPIX1"])*header["CDELT1"] : header["CDELT1"] : header["CRVAL1"]+(header["NAXIS1"]-header["CRPIX1"])*header["CDELT1"]).*convfactor
    yvec       = Array(header["CRVAL2"]+(1-header["CRPIX2"])*header["CDELT2"] : header["CDELT2"] : header["CRVAL2"]+(header["NAXIS2"]-header["CRPIX2"])*header["CDELT2"]).*convfactor
    
    delta_xvec = xvec.-xvec[trunc(Int,size(xvec)[1]/2)]
    delta_yvec = yvec.-yvec[trunc(Int,size(yvec)[1]/2)]
    return(xvec,yvec,delta_xvec,delta_yvec)
end




"""
    read_dat(DATFILEPATH ; com='#')

Return values of the .dat file given as input. By default, will consider as comment rows starting with a #. Can be changed.
"""
function read_dat(DATFILEPATH ; com='#')
    return(readdlm("$DATFILEPATH",'\t',Float64,comments=true,comment_char=com))

end


"""
    read_dim(arr)

Return a vector with dimensions of the data.
"""
function read_dim(arr)
    ndims(arr)<=3 && return(size(arr))
    return("Not valid dimensions (>3D)")
end



""" 
    read_fits_cvi(path ; check = true)

Return the cube, its dimension and its header. Use to read CVI cube, with this order : Pixel position, Rotation, LAG.
"""
function read_fits_cvi(path ; check = true)
    fitname    =  FITS("$(path)")
    header     =  read_header(fitname[1])
    cube       =  permcolrow(read(fitname[1]))
    cube       =  convert(Array{Float64,3},cube)
    check == true && valid_header(header)
    haskey(header,"BLANK")==1 && haskey(header,"BSCALE")==1 && haskey(header,"BZERO")==1 && (cube = replace_blanktomissing(cube,header["BLANK"]*header["BSCALE"]+header["BZERO"]))
    dimens     =  read_dim(cube)
    close(fitname)

    return(cube,dimens,header)
end



"""
    read_fits_pp(path)

Read a 2D fits (pixel*pixel). Return values (matrix), header and dimensions.
"""
function read_fits_pp(path)
    fitname     = FITS("$(path)")
    header      = read_header(fitname[1])
    map         = read(fitname[1])
    map         = permcolrow(map)
    dimens      = read_dim(map)
    haskey(header,"BLANK")==1 && haskey(header,"BSCALE")==1 && haskey(header,"BZERO")==1 && (map = replace_blanktomissing(map,header["BLANK"]*header["BSCALE"]+header["BZERO"]))
    close(fitname)

    return(map,header,dimens)
end



function read_fits_pv(path,vel_units ; check = true)
    fitname    =  FITS("$(path)")
    header     =  read_header(fitname[1])
    cube       =  read(fitname[1])
    cube       =  convert(Array{Float64,2},cube)
    check == true && valid_header(header)
    haskey(header,"BLANK")==1 && haskey(header,"BSCALE")==1 && haskey(header,"BZERO")==1 && (cube = replace_blanktomissing(cube,header["BLANK"]*header["BSCALE"]+header["BZERO"]))
    dimens     =  read_dim(cube)
    vel_units == "m/s"  && (convfactor = 1e-3)
    vel_units == "km/s" && (convfactor = 1)
    if haskey(header,"CDELT2")==true
        deltav = header["CDELT2"].*convfactor
        veloname = Array(header["CRVAL2"]+(1-header["CRPIX2"])*header["CDELT2"] : header["CDELT2"] : header["CRVAL2"]+(header["NAXIS2"]-trunc(Int,header["CRPIX2"]))*header["CDELT2"]).*convfactor
        return(cube,veloname,dimens,deltav,header)
    end
    close(fitname)

    return(cube,dimens,header)
end



"""
    read_fits_ppv(path,vel_units ; check = true)

Read data of a PPV fits from its path, test if the fits is conform, then return :
[1] Data in an array
[2] Range of velocity in an array (in km/s)
[3] Dimensions of the data (one vector with each element a dimension)
[4] The velocity resolution (in km/s)
[5] The header of the fits
If the header does not contain any third dimension indication, the function will still return the data, the dimension and the header.
Velocities need to be in the third dimension.
This function conduct a sanity check by default. Give check = false if you don't want the test being conducted.
"""
function read_fits_ppv(path,vel_units ; check = true)
    fitname    =  FITS("$(path)")
    header     =  read_header(fitname[1])
    cube       =  permcolrow(read(fitname[1]))
    cube       =  convert(Array{Float64,3},cube)
    check == true && valid_header(header)
    haskey(header,"BLANK")==1 && haskey(header,"BSCALE")==1 && haskey(header,"BZERO")==1 && (cube = replace_blanktomissing(cube,header["BLANK"]*header["BSCALE"]+header["BZERO"]))
    dimens     =  read_dim(cube)
    vel_units == "m/s"  && (convfactor = 1e-3)
    vel_units == "km/s" && (convfactor = 1)
    if haskey(header,"CDELT3")==true
        deltav = header["CDELT3"].*convfactor
        veloname = Array(header["CRVAL3"]+(1-header["CRPIX3"])*header["CDELT3"] : header["CDELT3"] : header["CRVAL3"]+(header["NAXIS3"]-trunc(Int,header["CRPIX3"]))*header["CDELT3"]).*convfactor
        return(cube,veloname,dimens,deltav,header)
    end
    close(fitname)

    return(cube,dimens,header)
end



"""
    read_fits_ppf(field ; path="none", check = true) WIP

Read data of a PPF fits from its path, test if the fits is conform, then return :
[1] Data in an array
[2] Range of velocity in an array (in km/s)
[3] X coordinate in arcsecond
[4] Y coordinate in arcsecond
[5] X coordinates centered in the middle of the image
[6] Y coordinates centered in the middle of the image
[7] Dimensions of the data (one vector with each element a dimension)
[8] The velocity resolution (in km/s)
[9] The header of the fits
[10] The path to the fits
Make the conversion of frequency into velocity.
If the header does not contain any third dimension indication, the function will still return the data, the dimension and the header.
Velocities need to be in the third dimension.
This function conduct a sanity check by default. Give check = false if you don't want the test being conducted.
"""
function read_fits_ppf(path ; check = true)
    fitname    =  FITS(fitsfile)
    header     =  read_header(fitname[1])
    cube       =  permcolrow(read(fitname[1]))
    cube       =  convert(Array{Float64,3},cube)
    check == true && valid_header(header)
    haskey(header,"BLANK")==1 && (cube = replace_blanktomissing(cube,header["BLANK"]*header["BSCALE"]+header["BZERO"]))
    dimens     =  read_dim(cube)
    #xvec,yvec,delta_xvec,delta_yvec = pixtocoord(header)

    freq       =  Array(header["CRVAL3"]+(1-header["CRPIX3"])*1 : 1 : header["CRVAL3"]+(header["NAXIS3"]-header["CRPIX3"])*1)
    dfreq      =  abs(freq[1]-freq[2])*1e9
    velo_vec   =  (freq.*1e9.-header["RESTFRQ"]).*-866.96338 ./dfreq .+ (header["VELOSYS"])
    # println("In which units are the velocity dimension in the fits ? km/s or m/s ? ")
    # velocityunits = readline()
    # velocityunits == "m/s"  && (convfactor = 1e-3)
    # velocityunits == "km/s" && (convfactor = 1)
    #haskey(header,"CDELT3") && (deltav = header["CDELT3"])
    #haskey(header,"CDELT3") && (veloname = Array(header["CRVAL3"]+(1-header["CRPIX3"])*header["CDELT3"] : header["CDELT3"] : header["CRVAL3"]+(header["NAXIS3"]-header["CRPIX3"])*header["CDELT3"]))
    #haskey(header,"NAXIS3") && (return(cube,veloname,xvec,yvec,delta_xvec,delta_yvec,dimens,deltav,header,fitsfile))
    close(fitname)

    return(cube,velo_vec,dimens,header)
end



"""
    read_fits_spect(path,vel_units )

Read data of a spectra PV fits from the path of the field and return :
[1] an array of the data
[2] an array of the range of velocity (in km/s)
[3] the velocity resolution (in km/s)
[4] the header of the fits
[5] the file of the fits

Velocities need to be in the second dimension.
By default, the path of the fits is "none", because we usually use local path listed in the function data_path.

"""
function read_fits_spect(path,vel_units )
    fitname  =  FITS("$(path)")
    header   =  read_header(fitname[1])
    data     =  permcolrow(read(fitname[1]))
    data     =  replace_blanktomissing(data,(header["BLANK"]*header["BSCALE"]+header["BZERO"]))

    vel_units == "m/s"  && (convfactor = 1e-3)
    vel_units == "km/s" && (convfactor = 1)
    deltav = header["CDELT2"].*convfactor
    veloname = Array(header["CRVAL2"]+(1-header["CRPIX2"])*header["CDELT2"] : header["CDELT2"] : header["CRVAL2"]+(header["NAXIS2"]-header["CRPIX2"])*header["CDELT2"]).*convfactor #m/s
    close(fitname)
    return(data,veloname,deltav,header)
end



""" 
    read_fits_velocityvector(path,vel_units)

Reconstruct the velocity vector of a fits file. Have to specify the units of the velocity vector of the fits.
"""
function read_fits_velocityvector(path,vel_units)
    fitname  = FITS("$(path)")
    header   = read_header(fitname[1])
    vel_units == "m/s"  && (convfactor = 1e-3)
    vel_units == "km/s" && (convfactor = 1)
    deltav = header["CDELT2"].*convfactor
    VelocityVectory = Array(header["CRVAL2"]+(1-header["CRPIX2"])*header["CDELT2"] : header["CDELT2"] : header["CRVAL2"]+(header["NAXIS2"]-header["CRPIX2"])*header["CDELT2"]).*convfactor #m/s
    close(fitname)
    return(VelocityVectory)
end



"""
    read_var_files(varfile_path)

A function to read and import values of a variable files. The variable files can be anywhere in the machine, but preferentially localised in the "var_file" folders.
"""
function read_var_files(varfile_path)
    #if isfile("$(varfile_path)")==false
    #    println("Your file does not exist. Would you like to create one ? (Y/N)")
    #    answer = readline()
    #    answer == "N" && (error("No .txt file at that location"))
    #    create_var_file(varfile_path)
    #end
    (isfile("$(varfile_path)")==false) && (error("No .txt file at that location"))
    var_data = readdlm("$(varfile_path)",comments=true)
    toreturn = Array{Any,2}(undef, size(var_data)[1], 1)
    for ix=1:size(var_data)[1] # Through the line
        toreturn[ix] = var_data[ix,2]
    end

    return(toreturn)
end



"""
    regular_blanking(data)

Check if the data of a fits file contains any irregular blanking e.g. a spectra with one
blanking value (or more, but less than the spectra size).
"""
function regular_blanking(data)
        any(broadcast(isequal,data,missing))==1 && (error("The sample data contains
                 irregular blanking, e.g. a pixel not blanked at all velocity channels."))
end



"""
    replace_blanktomissing(arr,blanktomissing)

Replace specific values (blanktomissing, generally blank value) by missing value in an array.
"""
function replace_blanktomissing(arr,blanktomissing)
    return(replace(arr,blanktomissing=>missing))
end



"""
    replace_inf_in_nan(arr)

Replace the Inf values by NaN values. Useful when searching for max and min values in the array.
"""
function replace_inf_in_nan(arr)
    arr = replace(arr,(-Inf)=>NaN)
    arr = replace(arr,Inf=>NaN)
    return(arr)
end



"""
    replace_missingtoblank(arr,missingtoblank)

Replace missing values by a 'missingtoblank' value in an array (generally a blank value).
"""
function replace_missingtoblank(arr,missingtoblank)
    return(replace(arr,missing=>missingtoblank))
end



"""
    replace_nantoblank(data,blank)

Replace the NaN values in an array by a new blank value.
"""
function replace_nantoblank(arr,blank)
    return(replace(arr,NaN=>blank))
end



"""
    replace_nantomissing(data)

Replace the NaN values in an array by missing value.
"""
function replace_nantomissing(arr)
    arr = replace_nantoblank(arr,-10000)
    arr = replace_blanktomissing(arr,-10000)
    return(arr)
end


"""
    function replace_nosignal(cube,DATADIMENSION,VELOCITYVECTOR,BLANK,SIGMAMAP)

Replace spectra without signal by blank value. We recover the spectra without signal by compute the intensity integration of each of them. If it is less than 2 times the dispersion of the noise, then it is blanked. SIGMAMAP is the 2D map of each rms noise of your cube ; use Data_analysis.rms_cube(cube,NOISECAN)[1] to compute it.
"""
function replace_nosignal(cube,DATADIMENSION,VELOCITYVECTOR,BLANK,SIGMAMAP)
    DELTAV = abs(VELOCITYVECTOR[2]-VELOCITYVECTOR[1])
    #integ = Array{Float64}(undef,DATADIMENSION[1]*DATADIMENSION[2])
    for pixc=1:DATADIMENSION[2]
        for pixr=1:DATADIMENSION[1]
            integ = sum(skipmissing(cube[pixr,pixc,:]))*DELTAV
            if integ<(2*SIGMAMAP[pixr,pixc])   # if integ<(2*SIGMAMAP[pixr,pixc])
                cube[pixr,pixc,:] .= missing
            end
        end
           
    end
    return(cube)

end


"""
        valid_header(header)

    Check if the header of a fits file contains a blank, a bzero, and a bscale values.
If not, return an error.
"""
function valid_header(header)
        haskey(header,"BLANK")==0  && (error("Not a valid header ; missing BLANK value"))
        haskey(header,"BZERO")==0  && (error("Not a valid header ; missing BZERO value"))
        haskey(header,"BSCALE")==0 && (error("Not a valid header ; missing BSCALE value"))
end


"""
    wait_for_key(; prompt = "Press any key", io = stdin))

Pause the script until the user press a key. The text prompt can be changed.
"""
function wait_for_key(; prompt = "Press any key", io = stdin)
    setraw!(raw) = ccall(:jl_tty_set_mode, Int32, (Ptr{Cvoid},Int32), io.handle, raw)
    println(io, prompt)
    setraw!(true)
    read(io, 1)
    setraw!(false)
    nothing
end



"""
    write_cv_fits(fitstocopy::String,newname::String,data,datadim,blank ; pc=0,varpercent=100)

Create a fits file used for cv reconstructed data. Contains on its header the number of pc used (default=0), the percentage of variance explained (default=100), the lag and the angle. Need a fits to copy its header entries.
"""
function write_cv_fits(fitstocopy::String,newname::String,pathtosave::String,data,datadim,blank ; pc=0,varpercent=100)
    # create a new fits, with a header based on fitstocopy
    eltype(data)==Union{Missing, Float64} && (data=convert(Array{Float64,2},replace_missingtoblank(data,blank)))
    secondfits = write_fits(fitstocopy,newname,pathtosave,data,datadim,blank;finished=false)
    # write in the header of the new fits
    write_key(secondfits[1],"#PC","$(pc)")
    write_key(secondfits[1],"PERCENTVAR","$(varpercent)")

    close(secondfits)
end



"""
    write_cvi_fits(lag,fitstocopy::String,newname::String,data,datadim,blank; pc=0,varpercent=100)

Create a fits used for cvi map. Contains on its header the lag used to do the cvi. Possibility to add in the header the number of PC used if data reconstructed from a PCA and the percentage of variance reproduced thank's to the PCA.
"""
function write_cvi_fits(lag,fitstocopy::String,newname::String,pathtosave::String,data,datadim,blank; pc=0,varpercent=100)
    # create a new fits, with a header based on fitstocopy
    eltype(data)==Union{Missing, Float64} && (data=convert(Array{Float64},replace_missingtoblank(data,blank)))
    secondfits = write_fits(fitstocopy,newname,pathtosave,data,datadim,blank;finished=false)
    # write in the header of the new fits
    write_key(secondfits[1],"#PC","$(pc)")
    write_key(secondfits[1],"PERCENTVAR","$(varpercent)")
    write_key(secondfits[1],"LAG",lag[1])
    close(secondfits)
end



""" 
   write_dat(matrix,PATHTOSAVE,NAME ; more=[""], overwrite=false)

Write matrix in a .dat file. Each column of the matrix is writen as new column in the file. Delimiters are spaces by default.
The option overwrite is false by default. If want to overwrite file with the same name as given in input, change it to true.
If want to add more entries on the header of the dat, use the option 'more'. Will be added as commentaries with a # at the begining of the line.
"""
function write_dat(matrix,PATHTOSAVE,NAME ; more=[""], overwrite=false,del='\t')
    # Check if there is already a file with the same name
    newname = NAME
    if (overwrite==false && isfile("$(PATHTOSAVE)/$(NAME).dat")==true)
        println("THE GIVEN FILE NAME ALREADY EXIST. AN INDICE WILL BE ADDED AT THE END OF THE GIVEN NAME, EQUAL TO THE NUMBER OF FILES WITH THE SAME NAME +1 ")
        count = 0
        for ix=1:size((findall.("$(NAME)",readdir("$(PATHTOSAVE)"))))[1]
            if size(findall("$(NAME)",readdir("$(PATHTOSAVE)")[ix]))[1]!=0
                count += 1
            end
        end
        newname = "$(NAME)_$(count)"
    end #if
    open("$(PATHTOSAVE)/$(newname).dat","w") do io
        if size(more)[1]>=1
            for kx = 1:size(more)[1]
                writedlm(io,["# $(more[kx])"],del)
            end
        end
        writedlm(io,matrix,del)
    end #open io
    println("File saved in $PATHTOSAVE/$newname.dat")
end #write_dat



"""
    write_fits(fitstocopy::String,newname::String,pathtosave::String,datatosave,datadim,blank ; finished=true)

Create a fits with its own header based on 'fitstocopy', and with the name 'newfits'. DO NOT ADD THE EXTENSION (e.g ".fits") FOR THE NAME OF THE NEWFITS). Need to give the blanking value. The fits can be constructed by data with two dimensions or three dimensions.
Variable 'finished' is by default true, so the function will by default close the new fits file at the end of the execution of this function. If 'finished' is false, the function will not close the fits file, allowing it to be modified again in another function (create_cvi_fits for example).
If want to add more entries on the header, use the option 'more'. Alternate between the name printed in the header and the value. Example, if want to add an entry for 'imsize' equals to 50, and an entry for 'powerlaw' equals to -3, write : more=["imsize","50","powerlaw","-3"].
The option overwrite is false by default. If want to overwrite file with the same name as given in input, change it to true.
"""
function write_fits(fitstocopy::String,newname::String,pathtosave::String,datatosave,datadim,blank ; finished=true , cvi=false , lags=0 , more=[""], overwrite=false)
    #Permute columns and rows ; needed to have the good orientation of data in result
    datatosave = permcolrow(datatosave)
    #Assign the header your new file will be based on to a variable
    fitname    = FITS(fitstocopy)
    head       = read_header(fitname[1])
    # Check if there is already a file with the same name
    if (overwrite==false && isfile("$(pathtosave)/$(newname).fits")==true)
        println("THE GIVEN FILE NAME ALREADY EXIST. AN INDICE WILL BE ADDED AT THE END OF THE GIVEN NAME, EQUAL TO THE NUMBER OF FILES WITH THE SAME NAME +1 ")
        count = 1
        for ix=1:size((findall.("$(newname)",readdir("$(pathtosave)"))))[1]
            if size(findall("$(newname)",readdir("$(pathtosave)")[ix]))[1]!=0
                count += 1
            end
        end
        newname = "$(newname)_$(count)"
    end
    #Create an empty fits file and assign the header to a variable
    secondfits = FITS("$pathtosave/$(newname).fits","w")
    write(secondfits,datatosave)
    #Copy the header in your new fits file
    #Don't take terms we don't want (e.g. dimension 3 description if 2D data)
    @inbounds @views for i = 1:length(head)
        typeof(datadim)==Tuple{Int64,Int64} && (read_key(fitname[1],i)[1]=="NAXIS3") && continue #will be written outside of this loop, depending on the dimension of the data
        typeof(datadim)==Tuple{Int64,Int64} && (read_key(fitname[1],i)[1]=="CROTA3") && continue
        typeof(datadim)==Tuple{Int64,Int64} && (read_key(fitname[1],i)[1]=="CTYPE3") && continue
        typeof(datadim)==Tuple{Int64,Int64} && (read_key(fitname[1],i)[1]=="CRVAL3") && continue
        typeof(datadim)==Tuple{Int64,Int64} && (read_key(fitname[1],i)[1]=="CDELT3") && continue
        typeof(datadim)==Tuple{Int64,Int64} && (read_key(fitname[1],i)[1]=="CRPIX3") && continue
        (read_key(fitname[1],i)[1]=="BSCALE") && continue
        (read_key(fitname[1],i)[1]=="BZERO")  && continue
        (read_key(fitname[1],i)[1]=="BITPIX") && continue
        (read_key(fitname[1],i)[1]=="END")    && break
        write_key(secondfits[1],read_key(fitname[1],i)[1],head[read_key(fitname[1],i)[1]])
    end
    replace_inf_in_nan(datatosave)
    #Write in the header specific attributes of the data
    #If 2 dims for the data, and reading NAXIS in the header, then assign 2 in that key.
    #If 3 dims for the data, and reading NAXIS in the header, then assign 3 in that key
    typeof(datadim)==Tuple{Int64,Int64} && haskey(head,"NAXIS") && write_key(secondfits[1],"NAXIS",2)
    typeof(datadim)==Tuple{Int64,Int64,Int64} && haskey(head,"NAXIS") && write_key(secondfits[1],"NAXIS",3)
    typeof(datadim)==Tuple{Int64,Int64,Int64} && haskey(head,"NAXIS3")  &&  write_key(secondfits[1],"NAXIS3",size(datatosave)[3])
    #haskey(head,"DATAMIN") &&  write_key(secondfits[1],"DATAMIN",minimum(datatosave))
    #haskey(head,"DATAMAX") &&  write_key(secondfits[1],"DATAMAX",maximum(datatosave))
    haskey(head,"NAXIS2")  &&  write_key(secondfits[1],"NAXIS2",size(datatosave)[2])
    haskey(head,"NAXIS1")  &&  write_key(secondfits[1],"NAXIS1",size(datatosave)[1])
    write_key(secondfits[1],"BLANK",blank)
    if size(more)[1]>=2
        for k = 1:2:size(more)[1]
            write_key(secondfits[1],"$(more[k])","$(more[k+1])")
        end
    end
    if cvi==true
        write_key(secondfits[1],"LAGS",join(lags,","))
    end
    close(fitname)
    if finished==true
        close(secondfits)
    else
        return(secondfits)
    end
end



"""
    write_pca_fits(pc,varpercent,fitstocopy::String,newname::String,data,datadim,blank)

DEPRECATED, USE WRITE_FITS INSTEAD, WITH OPTION MORE, SEE THE DOC. 
Create a fits used for pca reconstructed data. Contains on its header the number of pc used and the percentage of variance explained. Need a fits to copy its header entries.
"""
function write_pca_fits(pc,varpercent,fitstocopy::String,newname::String,pathtosave::String,data,datadim,blank)
    secondfits = write_fits(fitstocopy,newname,pathtosave,data,datadim,blank;finished=false)
    write_key(secondfits[1],"PC","$(pc)")
    write_key(secondfits[1],"PERCENTVAR","$(varpercent)")
    close(secondfits)
end




end #module
