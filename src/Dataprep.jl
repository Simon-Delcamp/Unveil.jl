module Dataprep 


export addblank
export blank_equal 
export blank_inf
export boolmatrix_missing
export directory_prep
export flatiterator
export pca_prep
export read_dat
export read_dim
export read_fits_cvi
export read_fits_pp
export read_fits_ppv
export regular_blanking
export replace_blanktomissing
export replace_inf_in_nan
export replace_missingtoblank
export replace_nantoblank
export replace_nantomissing
export replace_nosignal
export read_var_files
export write_dat
export write_fits

using FITSIO #, MultivariateStats, StatsBase, DelimitedFiles, StaticArrays
using DelimitedFiles

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
    boolmatrix_missing(array)

Return a 1D and a 2D boolean matrices indicating the indexes of missing values in array. Array should be a 2 dimension data (pixels*spectra for example)
"""
function boolmatrix_missing(array)
        missing2D = broadcast(isequal,array,missing)
        missing1D = findall(ismissing,array[:])
        return(missing1D, missing2D)
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
    (isdir("$(PATHTOSAVE)/Data/"))==0     &&  mkdir("$(PATHTOSAVE)/Data/")
    (isdir("$(PATHTOSAVE)/Figures/"))==0  &&  mkdir("$(PATHTOSAVE)/Figures/")
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
    pca_prep(arr,arraydimension)

Prepare data (given in 2D (pv) or 3D (ppv)) in order to conduct PCA on them by deleting missing values and reshape in 2D.

Will also check if any missing value still exist in the data, showing that the data are not regularly blanked. Return data in 2 dimensions without missing value, a 1D boolean matrix to check if any missing value still exist in the data (if 0 element then no missing value) and a 2D matrix with booleans corresponding to missing values in the dataset (used to reconstruct the data with the missing values).
"""
function pca_prep(arr,arraydimension)
        DATASIZETYPE = eltype(size(arr))
        (typeof(size(arr))==Tuple{DATASIZETYPE,DATASIZETYPE,DATASIZETYPE}) && (arr = reshape(arr,arraydimension[1]*arraydimension[2],arraydimension[3])) #3D -> 2D

        #println("produce boolean matrix")
        missing1D, missing2D = Dataprep.boolmatrix_missing(arr) #Boolean matrix to catch missing values
        #println("convert and delete")
        arr2D_nomissing     = Dataprep.deletemissing(arr,missing1D) #Data without missing value
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
    prodvarfile(;PATH=".",com=true)

Create the .txt varfiles needed as input of Unveil code. By default at the directory where Unveil is run. Can change the directory with option PATH. Comments can be removed with option 'com' set to false.
""" 
function prodvarfile(;PATH=".",com=true)
    if com==true
        open("$PATH/pca.txt","w") do io
            towrite = ["# FILE USED TO RUN FUNCTION Unveil.pca(). Construct a cube by using N PCs with the Principal Component Analysis method (PCA) on a noisy PPV cube .",
            "# COMMENTS ARE #",
            "#---------------------------------------------------------------------------------------------------------------#",
            "FITSPATH     \"PATH/TO/FITS\"   # Path of the fits. Its header will be used for the writting of the output fits file.",
            "FITSNAME         \"FITSNAME\"                                                # Name of the fits at the FITSPATH",
            "PATHTOSAVE       \"PATH/FOR/SAVING\"                                             # Path where Data and Plots will be saved",
            "SAVENAME	 \"SAVENAME\"					    # Generic name for the outputs. Extensions and attributes will be added by Unveil.",
            "UNITVELOCITY      m/s									         # Velocity units of the fits file",
            "NBPC              64                                              # Number of PC used for the reconstruction. ",
            "BLANK             -1000                                          # Blanking data",
            "NOISECAN	      1,25						 # Noise velocity canals",
            "OVERWRITE         false                                             # Would you like to overwrite output files with the same name ? true or false ",
            "#---------------------------------------------------------------------------------------------------------------#",
            ]
            writedlm(io,towrite, quotes=false)
        end #open io

        open("$PATH/convpca.txt","w") do io
            towrite = ["# FILE USED TO RUN FUNCTION Unveil.convpca(). Calculate the optimum number of PCs to used for the data reconstruction of a PPV cube (keeping most of the signals but removing  most of the noise)",
            "# COMMENTS ARE #",
            "#---------------------------------------------------------------------------------------------------------------#",
            "FITSPATH     \"PATH/TO/FITS\"   # Path of the fits. Its header will be used for the writting of the output fits file.",
            "FITSNAME         \"FITSNAME\"                                                # Name of the fits at the FITSPATH",
            "PATHTOSAVE       \"PATH/FOR/SAVING\"                                             # Path where Data and Plots will be saved",
            "SAVENAME	 \"SAVENAME\"					    # Generic name for the outputs. Extensions and attributes will be added by Unveil.",
            "NOISECAN          1,30                                           # Positions of the noise channels",
            "UNITVELOCITY      m/s									            # Velocity units of the fits file",
            "HIGHESTPC          20                                                # Maximum number of PCs to test the convergence",
            "BLANK             -1000                                             # Blanking data",
            "OVERWRITE         false                                             # Would you like to overwrite output files with the same name ? true or false ",
            "#---------------------------------------------------------------------------------------------------------------#",
            ]
            writedlm(io,towrite, quotes=false)
        end #open io

        open("$PATH/cv.txt","w") do io
            towrite = ["# FILE USED TO RUN FUNCTION Unveil.cv(). Calculate the Centroid Velocities (CV) on a PPV cube.",
            "# COMMENTS ARE #",
            "#---------------------------------------------------------------------------------------------------------------#",
            "FITSPATH     \"PATH/TO/FITS\"   # Path of the fits. Its header will be used for the writting of the output fits file.",
            "FITSNAME         \"FITSNAME\"                                                # Name of the fits at the FITSPATH",
            "PATHTOSAVE       \"PATH/FOR/SAVING\"                                             # Path where Data and Plots will be saved",
            "FITSOURCE       \"PATH/OF/FITS\"              # Path to the original fits. This is the cube without PCA nor SWO treatment, so with full noise.",
            "SAVENAME	 \"SAVENAME\"					    # Generic name for the outputs. Extensions and attributes will be added by Unveil.",
            "UNITVELOCITY      m/s									                    # Velocity units of the fits file",
            "THREHSOLD         0.2                                                       # Values lower than this factor multiplied by the noise RMS will be blanked. ",
            "NOISECANTXT       1,25                                                      # Noise channels ",
            "BLANK             -1000                                                     # Blanking data",
            "OVERWRITE         false                                             # Would you like to overwrite output files with the same name ? true or false ",
            "#---------------------------------------------------------------------------------------------------------------#",
            ]
            writedlm(io,towrite, quotes=false)
        end #open io

        open("$PATH/cvi.txt","w") do io
            towrite = ["# FILE USED TO RUN FUNCTION Unveil.cvi(). Compute the Centroid Velocity Increments (CVI) on a Cendroid Velocity (CV) map, at multiple lags.",
            "# COMMENTS ARE #",
            "#---------------------------------------------------------------------------------------------------------------#",
            "FITSPATH     \"PATH/TO/FITS\"   # Path of the fits. Its header will be used for the writting of the output fits file.",
            "FITSNAME         \"FITSNAME\"                                                # Name of the CVMAP at the FITSPATH",
            "PATHTOSAVE       \"PATH/FOR/SAVING\"                                             # Path where Data and Plots will be saved",
            "SAVENAME	 \"SAVENAME\"					    # Generic name for the outputs. Extensions and attributes will be added by Unveil.",
            "BLANK             -1000                                                     # Blanking data",
            "LAG               2,3,5,7,10,20,50                                          # Values of the Lags used for CVI calculations",
            "DIFFTYPE          relative                                                  # Type of differences computed during CVI calculation (relative or abs)",
            "OVERWRITE         false                                             # Would you like to overwrite output files with the same name ? true or false ",
            "#---------------------------------------------------------------------------------------------------------------#",
            ]
            writedlm(io,towrite, quotes=false)
        end #open io

        open("$PATH/cvcvi.txt","w") do io
            towrite = ["# FILE USED TO RUN FUNCTION Unveil.cvcvi(). Compute the Centroid Velocity (CV) of a PPV cube, and then the Centroid Velocity Increments (CVI) at multiple lags.",
            "# COMMENTS ARE #",
            "#---------------------------------------------------------------------------------------------------------------#",
            "FITSPATH     \"PATH/TO/FITS\"   # Path of the fits. Its header will be used for the writting of the output fits file.",
            "FITSNAME         \"FITSNAME\"                                                # Name of the fits at the FITSPATH",
            "PATHTOSAVE       \"PATH/FOR/SAVING\"                                             # Path where Data and Plots will be saved",
            "FITSOURCE       \"PATH/OF/FITS\"              # Path to the original fits. This is the cube without PCA nor SWO treatment, so with full noise.",
            "SAVENAME	 \"SAVENAME\"					    # Generic name for the outputs. Extensions and attributes will be added by Unveil.",
            "THREHSOLD         0.2                                                       # Values lower than this factor multiplied by the noise RMS will be blanked. ",
            "NOISECANTXT       1,25                                                      # Noise channels ",
            "UNITVELOCITY      m/s									            # Velocity units of the fits file",
            "REMOVE            false                                             # Remove spectra with very low intensity. ",
            "BLANK             -1000                                             # Blanking data",
            "LAG               2,3,5,7,10,20,50                                  # Values of the Lags used for CVI calculations",
            "DIFFTYPE          relative                                          # Type of differences computed during CVI calculation (relative or abs)",
            "OVERWRITE         false                                             # Would you like to overwrite output files with the same name ? true or false ",
            "#---------------------------------------------------------------------------------------------------------------#",
            ]
            writedlm(io,towrite, quotes=false)

            
        end #open io

        open("$PATH/swo.txt","w") do io
            towrite = ["# FILE USED TO RUN FUNCTION Unveil.swo(). Construct a cube by using the Spectral Window Optimisation method (SWO) on a noisy cube.",
            "# COMMENTS ARE #",
            "#---------------------------------------------------------------------------------------------------------------#",
            "FITSPATH     \"PATH/TO/FITS\"   # Path of the fits. Its header will be used for the writting of the output fits file.",
            "FITSNAME         \"FITSNAME\"                                                # Name of the fits at the FITSPATH",
            "PATHTOSAVE       \"PATH/FOR/SAVING\"                                             # Path where Data and Plots will be saved",
            "SAVENAME	 \"SAVENAME\"					    # Generic name for the outputs. Extensions and attributes will be added by Unveil.",
            "UNITVELOCITY      m/s									         # Velocity units of the fits file",
            "BLANK             -1000                                          # Blanking data",
            "NOISECAN          1,30                                           # Positions of the noise channels",
            "EXAMPLES          YES                                            # YES or NO. If yes, will plot 3 figures of 16 randomly chosen spectra each of the new cube in front of the source cube.",
            "OVERWRITE         false                                             # Would you like to overwrite output files with the same name ? true or false ",
            
            "#---------------------------------------------------------------------------------------------------------------#",
            ]
            writedlm(io,towrite, quotes=false)
        end #open io


        open("$PATH/compmethod_stcfct.txt","w") do io
            towrite = ["# FILE USED TO RUN FUNCTION Unveil.compmethod_stcfct(). Plot exponents of the structures functions computed on cube treated by multiple methods : PCA, SWO, Noise Free (if coming from a simulation) and Raw (if without data pre-treatment).",
            "# COMMENTS ARE #",
            "#---------------------------------------------------------------------------------------------------------------#",
            "DATPATH		\"/PATH/TO/FITS\" # Path to the FITS",
            "PCNAME		\"PCDATA.dat\" #Name of the structure functions dat file computed from the cube reconstructed by PCA. Let blank if none.",
            "SWONAME 	\"SWODATA.dat\" #Name of the structure functions dat file computed from the cube reconstructed by SWO. Let blank if none.",
            "NFNAME      \"NOISEFREEDATA.dat\" #Name of the structure functions dat file computed from the cube without noise (from simulation). Let blank if none.",
            "RAWNAME		\"RAWDATA\" #Name of the structure functions dat file computed from a cube reconstructed by PCA. Let blank if none.",
            "PATHTOSAVE      \"PATH/FOR/SAVING\" # Path where Data and Plots will be saved",  
            "SAVENAME	\"NAME\"  # Generic name for the outputs. Extensions and attributes will be added by Unveil.",
            "OVERWRITE	true # Would you like to overwrite output files with the same name ? true or false ",
            "#---------------------------------------------------------------------------------------------------------------#",
            ]
            writedlm(io,towrite, quotes=false)
        end #open io

        open("$PATH/structure_functions.txt","w") do io
            towrite = ["# FILE USED TO RUN FUNCTION Unveil.structure_functions(). Compute structure functions and exponents from a CVI cube (better to consider all rotations, not an azimutal average)",
            "# COMMENTS ARE #",
            "#---------------------------------------------------------------------------------------------------------------#",
            "FITSPATH     \"PATH/TO/FITS\"   # Path of the fits. Its header will be used for the writting of the output fits file.",
            "FITSNAME         \"FITSNAME\"                                                # Name of the fits at the FITSPATH",
            "PATHTOSAVE       \"PATH/FOR/SAVING\"                                             # Path where Data and Plots will be saved",
            "SAVENAME	 \"SAVENAME\"					    # Generic name for the outputs. Extensions and attributes will be added by Unveil.",
            "ORDERS		  1,2,3,4,5,6                   # Orders of the structures functions",
            "BLANK             -1000                                                     # Blanking data",
            "OVERWRITE         false                                             # Would you like to overwrite output files with the same name ? true or false ",
            "#---------------------------------------------------------------------------------------------------------------#",
            ]
            writedlm(io,towrite, quotes=false)
        end #open io


    else 
        open("$PATH/pca.txt","w") do io
            towrite = [
            "FITSPATH     \"PATH/TO/FITS\"   ",
            "FITSNAME         \"FITSNAME\"     ",
            "PATHTOSAVE       \"PATH/FOR/SAVING\"        ",
            "SAVENAME	 \"SAVENAME\"				",
            "UNITVELOCITY      m/s					",
            "NBPC              64                 ",
            "BLANK             -1000              ",
            "NOISECAN	      1,25						",
            "OVERWRITE         false          ",
           ]
            writedlm(io,towrite, quotes=false)
        end #open io

        open("$PATH/convpca.txt","w") do io
            towrite = [
            "FITSPATH     \"PATH/TO/FITS\"  ",
            "FITSNAME         \"FITSNAME\"   ",
            "PATHTOSAVE       \"PATH/FOR/SAVING\"     ",
            "SAVENAME	 \"SAVENAME\"				",
            "NOISECANTXT       1,25         ",
            "UNITVELOCITY      m/s				",
            "HIGHESTPC          20              ",
            "BLANK             -1000            ",
            "OVERWRITE         false            ",
            ]
            writedlm(io,towrite, quotes=false)
        end #open io

        open("$PATH/cv.txt","w") do io
            towrite = [
            "FITSPATH     \"PATH/TO/FITS\"   ",
            "FITSNAME         \"FITSNAME\"    ",
            "PATHTOSAVE       \"PATH/FOR/SAVING\"     ",
            "FITSOURCE       \"PATH/OF/FITS\"       ",
            "SAVENAME	 \"SAVENAME\"	",
            "UNITVELOCITY      m/s		",
            "THREHSOLD         0.2    ",
            "NOISECANTXT       1,25            ",
            "BLANK             -1000        ",
            "OVERWRITE         false      ",
            ]
            writedlm(io,towrite, quotes=false)
        end #open io

        open("$PATH/cvi.txt","w") do io
            towrite = [
            "FITSPATH     \"PATH/TO/FITS\"  ",
            "FITSNAME         \"FITSNAME\"                  ",
            "PATHTOSAVE       \"PATH/FOR/SAVING\"         ",
            "SAVENAME	 \"SAVENAME\"					   ",
            "BLANK             -1000                        ",
            "LAG               2,3,5,7,10,20,50          ",
            "DIFFTYPE          relative              ",
            "OVERWRITE         false            ",
            ]
            writedlm(io,towrite, quotes=false)
        end #open io

        open("$PATH/cvcvi.txt","w") do io
            towrite = [
            "FITSPATH     \"PATH/TO/FITS\" ",
            "FITSNAME         \"FITSNAME\"  ",
            "PATHTOSAVE       \"PATH/FOR/SAVING\"    ",
            "FITSOURCE       \"PATH/OF/FITS\"       ",
            "SAVENAME	 \"SAVENAME\"			",
            "THREHSOLD         0.2         ",
            "NOISECANTXT       1,25         ",
            "UNITVELOCITY      m/s			",
            "REMOVE            false       ",
            "BLANK             -1000           ",
            "LAG               2,3,5,7,10,20,50        ",
            "DIFFTYPE          relative                 ",
            "OVERWRITE         false          ",
            ]
            writedlm(io,towrite, quotes=false)

            
        end #open io

        open("$PATH/swo.txt","w") do io
            towrite = [
            "FITSPATH     \"PATH/TO/FITS\"  ",
            "FITSNAME         \"FITSNAME\"      ",
            "PATHTOSAVE       \"PATH/FOR/SAVING\"     ",
            "SAVENAME	 \"SAVENAME\"					",
            "UNITVELOCITY      m/s					",
            "BLANK             -1000                ",
            "NOISECAN          1,30                 ",
            "EXAMPLES          YES                 ",
            "OVERWRITE         false               ",
            ]
            writedlm(io,towrite, quotes=false)
        end #open io


        open("$PATH/compmethod_stcfct.txt","w") do io
            towrite = [
            "DATPATH		\"/PATH/TO/FITS\" ",
            "PCNAME		\"PCDATA.dat\" ",
            "SWONAME 	\"SWODATA.dat\" ",
            "NFNAME      \"NOISEFREEDATA.dat\" ",
            "RAWNAME		\"RAWDATA\" ",
            "PATHTOSAVE      \"PATH/FOR/SAVING\" ",  
            "SAVENAME	\"NAME\" ",
            "OVERWRITE	true ",
            ]
            writedlm(io,towrite, quotes=false)
        end #open io

        open("$PATH/structure_functions.txt","w") do io
            towrite = [
            "FITSPATH     \"PATH/TO/FITS\"  ",
            "FITSNAME         \"FITSNAME\"    ",
            "PATHTOSAVE       \"PATH/FOR/SAVING\"      ",
            "SAVENAME	 \"SAVENAME\"		",
            "ORDERS		  1,2,3,4,5,6    ",
            "BLANK             -1000      ",
            "OVERWRITE         false       ",
            ]
            writedlm(io,towrite, quotes=false)
        end #open io
    end



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
        #write_key(secondfits[1],"LAGS",join(lags,","))
        haskey(head,"CTYPE1")  &&  write_key(secondfits[1],"CTYPE1","$(read_key(fitname[1],"CTYPE1")[1]),$(read_key(fitname[1],"CTYPE2")[1])")
        haskey(head,"CTYPE2")  &&  write_key(secondfits[1],"CTYPE2","ROTATION -- DEPEND OF THE LAG (cf comments)")
        haskey(head,"CTYPE3")  &&  write_key(secondfits[1],"CTYPE3","LAG -- PIXEL")
        haskey(head,"CRVAL1")  &&  write_key(secondfits[1],"CRVAL1",0)
        haskey(head,"CRPIX1")  &&  write_key(secondfits[1],"CRPIX1",0)
        haskey(head,"CRDELT1")  &&  write_key(secondfits[1],"CRDEL1",0)
        haskey(head,"CRVAL2")  &&  write_key(secondfits[1],"CRVAL2",0)
        haskey(head,"CRPIX2")  &&  write_key(secondfits[1],"CRPIX2",0)
        haskey(head,"CRDELT2")  &&  write_key(secondfits[1],"CRDEL2",0)
        haskey(head,"CRVAL3")  &&  write_key(secondfits[1],"CRVAL3",0)
        haskey(head,"CRPIX3")  &&  write_key(secondfits[1],"CRPIX3",0)
        haskey(head,"CRDELT3")  &&  write_key(secondfits[1],"CRDEL3",0)
        haskey(head,"BUNIT")  &&  haskey(head,"CTYPE3") && write_key(secondfits[1],"BUNIT","$(read_key(fitname[1],"CTYPE3")[1])")

        write_key(secondfits[1],"COMMENTS1","Each lag gives a specific number of rotations ")
        write_key(secondfits[1],"COMMENTS2","(e.g l=3 needs less rotations than l=50).")
        write_key(secondfits[1],"COMMENTS3","Spectra are not ordered simply")
        write_key(secondfits[1],"COMMENTS4","thus can't reproduce the image in 2D with that file.")

    end
    close(fitname)
    if finished==true
        close(secondfits)
    else
        return(secondfits)
    end
end



end #module