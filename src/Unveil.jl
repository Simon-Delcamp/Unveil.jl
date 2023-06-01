module Unveil

import Pkg
import Plots
import DelimitedFiles
import FITSIO
import StatsBase
import Distributions
import Measures
import MultivariateStats
import Makie
import GLMakie
import StaticArrays
import LaTeXStrings
import StatsPlots
import ShiftedArrays
import KernelDensity
import LsqFit
import Formatting
import PyPlot
import Format
import CSV
import DataFrames
import Interact
import CurveFit
import Documenter
import Statistics
import Colors
import Profile

include("Data_preparation.jl") # Read and write fits
include("Data_analysis.jl") 
include("Functionforpca.jl")   # Functions for method PCA
include("Spectralwindowopti.jl")  # Functions for method SWO
include("Graphic.jl")   # Plotting 
include("Functionforcvi.jl")   # Functions for CVI computations

using .Data_preparation
using .Functionforpca
using .Spectralwindowopti
using .Graphic
using .Data_analysis

export pca
export swo
export convpca
export convswo
export cvcvi
export cvi 
export multipca
export combinecv



"""
    convpca(VARFILEPATH)

Produce calculations to find the PCA convergence criteria based on the matrix projection from PCA. No PCA reconstructed cube will be saved. A '.txt' file should be used accordingly as an input (see models inside folders '/varfiles/convpca.txt').

OUTPUTS : A plots with every moments of the projection matrix + the metric (see the doc). Also add a .dat file with moments, metric and number of PC used for each.

Use this script in a julia terminal with :
    julia>Unveil.convpca(VARFILEPATH)
"""
function convpca(VARFILEPATH)
    FITSPATH,FILENAME,PATHTOSAVE,UNITVELOCITY,HIGHESTPC,BLANK = read_var_files(VARFILEPATH)

    # Read the fits from the path. Return the data, the VelocityVector, the dimension, the velocity_increment, and the header.
    cube,VELOCITYVECTOR,DATADIMENSION,VELOCITYINCREMENT,HEAD = Data_preparation.read_fits_ppv("$(FITSPATH)/$(FILENAME)",UNITVELOCITY ; check=false)

    # Prepare directories where plots and data will be saved.
    Data_preparation.directory_prep(PATHTOSAVE)


    # Replace any NaN value into a missing value and deleted them (can't do PCA on missing values with that package).
    cube = Data_preparation.replace_nantomissing(cube)
    ismis = 0
    if any(ismissing,cube) 
        ismis = 1
        cube,missingplaces1D,missingplaces2D  = Data_preparation.pca_prep(cube,DATADIMENSION)
        cube                                  = convert(Array{Float64},cube)
        DATADIMENSION_NOMISSING               = Data_preparation.read_dim(cube)
    else
        cube                                 = reshape(cube,DATADIMENSION[1]*DATADIMENSION[2],DATADIMENSION[3])
        cube                                 = convert(Array{Float64},cube)
        DATADIMENSION_NOMISSING              = (DATADIMENSION[1]*DATADIMENSION[2],DATADIMENSION[3])
    end

    #First PCA
    println("Perform PCA")
    M, Yt, VARPERCENT,cubereconstructed = Functionforpca.pca(cube,HIGHESTPC)

    proj = MultivariateStats.projection(M)
    mom1,mom2,mom3,mom4 = Data_analysis.fourmoments(proj,dim=2)
    metric = Data_analysis.calcmetric(mom1,mom2,mom3,mom4,abs(VELOCITYINCREMENT))

    if ismis == 1
        proj = Data_preparation.addblank(proj,missingplaces2D[:,1:HIGHESTPC],BLANK,(DATADIMENSION[1],DATADIMENSION[2],HIGHESTPC))
    end
    proj = reshape(proj,(DATADIMENSION[1],DATADIMENSION[2],HIGHESTPC))
    Data_preparation.write_fits("$(FITSPATH)/$FILENAME","projectionmatrix","$(PATHTOSAVE)/Data/",proj,(DATADIMENSION_NOMISSING[1],DATADIMENSION_NOMISSING[2],HIGHESTPC),BLANK,finished=true,overwrite=true)
    xvector = range(1,HIGHESTPC)#[1:HIGHESTPC]
    Graphic.distribcv_multipc(mom1[2:end-1],mom2[2:end-1],mom3[2:end-1],mom4[2:end-1],metric[2:end-1],xvector[2:end-1])
    Plots.savefig("$(PATHTOSAVE)/pca_metric.pdf")

    Data_preparation.write_dat([metric mom1 mom2 mom3 mom4 xvector],"$PATHTOSAVE","metricPCA",overwrite=true,more=["$FILENAME","Metric  Mom1   Mom2   Mom3   Mom4   PCs"])
    minimetr = minimum(metric)
    minipc   = xvector[findall(x->x==minimetr,metric)]
    println("Metric minimum=$(minimetr)")
    println("#PC of metric minimum=$(minipc)")
end #convpca




"""
    convswo(VARFILEPATH)

Look for the best interval to integrate the spectra for the window optimisation process. A '.txt' file should be used accordingly as an input (see models inside folders '/varfiles/convswo.txt').

Use this script in a julia terminal with :
    julia>Unveil.convswo(VARFILEPATH)
"""
function convswo(VARFILEPATH)
    FITSPATH,FILENAME,PATHTOSAVE,UNITVELOCITY,RANGETXT,BLANK,NOISECANTXT = read_var_files(VARFILEPATH)
    NOISECAN = [parse(Int, ss) for ss in split(NOISECANTXT,",")]
    RANGE = [parse(Int, ss) for ss in split(RANGETXT,",")]
    RANGE = [RANGE[1]:RANGE[2]]


    # Read the fits from the path. Return the data, the VelocityVector, the dimension, the velocity_increment, and the header.
    cube,VELOCITYVECTOR,DATADIMENSION,VELOCITYINCREMENT,HEAD = Data_preparation.read_fits_ppv("$(FITSPATH)/$(FILENAME)",UNITVELOCITY ; check=false)

    # Prepare directories where plots and data will be saved.
    Data_preparation.directory_prep(PATHTOSAVE)


    # Replace any NaN value into a missing value and deleted them (can't do PCA on missing values with that package).
    cube = Data_preparation.replace_nantomissing(cube)
    ismis = 0
    if any(ismissing,cube) 
        ismis = 1
        cube,missingplaces1D,missingplaces2D  = Data_preparation.pca_prep(cube,DATADIMENSION)
        cube                                  = convert(Array{Float64},cube)
        DATADIMENSION_NOMISSING               = Data_preparation.read_dim(cube)
    else
        missingplaces2D = 0
        cube                                 = reshape(cube,DATADIMENSION[1]*DATADIMENSION[2],DATADIMENSION[3])
        cube                                 = convert(Array{Float64},cube)
        DATADIMENSION_NOMISSING              = (DATADIMENSION[1]*DATADIMENSION[2],DATADIMENSION[3])
    end


    muc,mus,sigc,sigs,gamc,gams,kapc,kaps,metc,mets,SIGMAT=Spectralwindowopti.convoptiwind(cube,DATADIMENSION_NOMISSING,DATADIMENSION,VELOCITYVECTOR,NOISECAN,BLANK,RANGE[1],PATHTOSAVE,missingplaces2D,ismiss=ismis)

    Graphic.distribcv_multiow(muc[1:end],sigc[1:end],gamc[1:end],kapc[1:end],metc[1:end],RANGE[1][1:end],"Parameters of the distribution of the differences \n of Tdv between source cube and window optimized cube",SIGMAT=SIGMAT)
    Plots.savefig("$PATHTOSAVE/swo_difTDV.pdf")

    Graphic.distribcv_multiow(mus[1:end],sigs[1:end],gams[1:end],kaps[1:end],mets[1:end],RANGE[1][1:end],"Parameters of distribution of the percentage of flagged velocity canals",SIGMAT=0)
    Plots.savefig("$PATHTOSAVE/swo_percFLAGGED.pdf")

    Data_preparation.write_dat([mets mus sigs gams kaps RANGE[1][:]],"$PATHTOSAVE","metricSWO",overwrite=true,more=["$FILENAME","Metric  Mom1   Mom2   Mom3   Mom4   RANGE"])

end #convswo




"""
    cv(VARFILEPATH)

Calculate the CV from a cube given in input. A '.txt' file should be used accordingly as an input (see models inside folders '/varfiles/cv.txt').

Use this script in a julia terminal with :
    julia>Unveil.cv(VARFILEPATH)
"""
function cv(VARFILEPATH)
    FITSPATH,FITSNAME,PATHTOSAVE,UNITVELOCITY,NBPC,BLANK= read_var_files(VARFILEPATH)
    LAG = [parse(Int, ss) for ss in split(LAG,",")]
    if NBPC == 0
        NBPC="raw"
    elseif NBPC == -1
        NBPC="SWO"
    elseif NBPC > 0
        NBPC="$(NBPC)PC"
    end


    # Read the fits from the path. Return the data, the VelocityVector, the dimension, the velocity_increment, and the header.
    cube,VELOCITYVECTOR,DATADIMENSION,VELOCITYINCREMENT,HEAD = Data_preparation.read_fits_ppv("$FITSPATH/$FITSNAME",UNITVELOCITY ; check=false)

    # Prepare directories where plots and data will be saved.
    Data_preparation.directory_prep(PATHTOSAVE)


    # Replace any NaN value into a missing value and deleted them (can't do PCA on missing values with that package).
    cube = Data_preparation.replace_nantomissing(cube)
    ismis = 0
    if any(ismissing,cube) 
        ismis = 1
        cube,missingplaces1D,missingplaces2D  = Data_preparation.pca_prep(cube,DATADIMENSION)
        cube                                  = convert(Array{Float64},cube)
        DATADIMENSION_NOMISSING               = Data_preparation.read_dim(cube)
    else
        cube                                 = reshape(cube,DATADIMENSION[1]*DATADIMENSION[2],DATADIMENSION[3])
        cube                                 = convert(Array{Float64},cube)
        DATADIMENSION_NOMISSING              = (DATADIMENSION[1]*DATADIMENSION[2],DATADIMENSION[3])
    end


    println("------ CV CALCULATION ------")
    cvmap = Functionforcvi.moment_one_field(cube,VELOCITYVECTOR,BLANK) # Calculate the first velocity moment order on data reconstructed
    cube  = 0.0 
    GC.gc()

    if ismis == 1
        cvmap = Data_preparation.addblank(cvmap,missingplaces2D[:,1],BLANK,(DATADIMENSION[1],DATADIMENSION[2]))
    end
    cvmap = reshape(cvmap,(DATADIMENSION[1],DATADIMENSION[2]))


    Data_preparation.write_fits("$(FITSPATH)/$FITSNAME","cv_$(NBPC)","$(PATHTOSAVE)/Data/",cvmap,(DATADIMENSION_NOMISSING[1],DATADIMENSION_NOMISSING[2]),BLANK,finished=true)
    #cvmap = 0.0

    println("CV map saved in $(PATHTOSAVE)/Data/cv_$(NBPC)_NumberOfFilesWithTheSameNameAsPrefixe.fits as a fits.")


end #cv


"""
    cvcvi(VARFILEPATH)

Calculate the CV and CVI from a cube given in input. A '.txt' file should be used accordingly as an input (see models inside folders '/varfiles/cvcvi.txt').

Use this script in a julia terminal with :
    julia>Unveil.cvcvi(VARFILEPATH)
"""
function cvcvi(VARFILEPATH)
    FITSPATH,FITSNAME,PATHTOSAVE,UNITVELOCITY,NBPC,BLANK,LAG,DIFFTYPE = read_var_files(VARFILEPATH)
    LAG = [parse(Int, ss) for ss in split(LAG,",")]
    if NBPC == 0
        NBPC="raw"
    elseif NBPC == -1
        NBPC="SWO"
    elseif NBPC > 0
        NBPC="$(NBPC)PC"
    end


    # Read the fits from the path. Return the data, the VelocityVector, the dimension, the velocity_increment, and the header.
    cube,VELOCITYVECTOR,DATADIMENSION,VELOCITYINCREMENT,HEAD = Data_preparation.read_fits_ppv("$FITSPATH/$FITSNAME",UNITVELOCITY ; check=false)

    # Prepare directories where plots and data will be saved.
    Data_preparation.directory_prep(PATHTOSAVE)


    # Replace any NaN value into a missing value and deleted them (can't do PCA on missing values with that package).
    cube = Data_preparation.replace_nantomissing(cube)
    ismis = 0
    if any(ismissing,cube) 
        ismis = 1
        cube,missingplaces1D,missingplaces2D  = Data_preparation.pca_prep(cube,DATADIMENSION)
        cube                                  = convert(Array{Float64},cube)
        DATADIMENSION_NOMISSING               = Data_preparation.read_dim(cube)
    else
        cube                                 = reshape(cube,DATADIMENSION[1]*DATADIMENSION[2],DATADIMENSION[3])
        cube                                 = convert(Array{Float64},cube)
        DATADIMENSION_NOMISSING              = (DATADIMENSION[1]*DATADIMENSION[2],DATADIMENSION[3])
    end


    println("------ CV CALCULATION ------")
    cvmap = Functionforcvi.moment_one_field(cube,VELOCITYVECTOR,BLANK) # Calculate the first velocity moment order on data reconstructed
    cube  = 0.0 
    GC.gc()

    if ismis == 1
        cvmap = Data_preparation.addblank(cvmap,missingplaces2D[:,1],BLANK,(DATADIMENSION[1],DATADIMENSION[2]))
    end
    cvmap = reshape(cvmap,(DATADIMENSION[1],DATADIMENSION[2]))


    Data_preparation.write_fits("$(FITSPATH)/$FITSNAME","cv_$(NBPC)","$(PATHTOSAVE)/Data/",cvmap,(DATADIMENSION_NOMISSING[1],DATADIMENSION_NOMISSING[2]),BLANK,finished=true)
    #cvmap = 0.0

    GC.gc()
    println("CV map saved in $(PATHTOSAVE)/Data/cv_$(NBPC)_NumberOfFilesWithTheSameNameAsPrefixe.fits as a fits.")

    println("------ CVI CALCULATION ------")
    if DIFFTYPE=="relative" 
        cviallangle,cvimap_averaged = Functionforcvi.construct_cvimap(cvmap,LAG,(DATADIMENSION[1],DATADIMENSION[2]),diff="relative")
        DIFFTYPE = "rel"

    elseif DIFFTYPE=="abs" 
        cviallangle,cvimap_averaged = Functionforcvi.construct_cvimap(cvmap,LAG,(DATADIMENSION[1],DATADIMENSION[2]),diff="absolute")
        DIFFTYPE = "abs"

    else
        error("Not good argument in DIFFTYPE (should be abs or relative)")
    end


    cvimap_averaged = reshape(cvimap_averaged,DATADIMENSION[1],DATADIMENSION[2],size(LAG)[1])
    cvimap_averaged = Data_preparation.replace_missingtoblank(cvimap_averaged,BLANK)
    cvimap_averaged = Data_preparation.blank_equal(cvimap_averaged,0.0,BLANK)
    cvimap_averaged = convert(Array{Float64},cvimap_averaged)

    ## DO NOT NEED TO RECONSTRUCT THE MAP WITH THE BLANK VALUES. ONLY USED TO COMPUTE PDF, SO ONLY Statistics
    cviallangle = Data_preparation.replace_nantoblank(cviallangle,BLANK)
    cviallangle = Data_preparation.replace_missingtoblank(cviallangle,BLANK)
    cviallangle = Data_preparation.blank_equal(cviallangle,0.0,BLANK)
    cviallangle = convert(Array{Float64},cviallangle)

    Data_preparation.write_fits("$(FITSPATH)/$FITSNAME","cvi$(DIFFTYPE)_multlag_METH$(NBPC)","$(PATHTOSAVE)/Data/",cvimap_averaged,(DATADIMENSION[1],DATADIMENSION[2],size(LAG)[1]),BLANK,finished=true,cvi=true,lags=LAG)
    println("CVI map reconstructed from PCA saved in $(PATHTOSAVE)/Data/cvi$(DIFFTYPE)_multlag_$(NBPC)PC_NumberOfFilesWithTheSameNameAsPrefixe.fits as a fits.")

    Data_preparation.write_fits("$(FITSPATH)/$FITSNAME","cvi$(DIFFTYPE)_multlag_allangle_METH$(NBPC)","$(PATHTOSAVE)/Data",cviallangle,(DATADIMENSION[1]*DATADIMENSION[2],NANGLE,size(LAG)[1]),BLANK,finished=true)
    println("CVI map with all angles values reconstructed from PCA saved in the $(PATHTOSAVE)/Data/cvi$(DIFFTYPE)_multlag_allangle_METH$(NBPC)_NumberOfFilesWithTheSameNameAsPrefixe.fits as a fits.")

end #function cvcvi




"""
    cvi(VARFILEPATH)

Calculate the CVI from a CV map. A '.txt' file should be used accordingly as an input (see models inside folders '/varfiles/cvi.txt'). 

INPUT : path to the variable '.txt' file
OUTPUTS : Save the CV fits, 

Use this script in a julia terminal with :
    julia>Unveil.cvi(VARFILEPATH)
"""
function cvi(VARFILEPATH)
    FITSPATH,FITSNAME,PATHTOSAVE,UNITVELOCITY,NBPC,BLANK,LAG,DIFFTYPE = read_var_files(VARFILEPATH)
    LAG = [parse(Int, ss) for ss in split(LAG,",")]
    if NBPC == 0
        NBPC="raw"
    elseif NBPC == -1
        NBPC="SWO"
    elseif NBPC > 0
        NBPC="$(NBPC)PC"
    end

    cvmap,HEAD,DATADIMENSION = Data_preparation.read_fits_pp("$FITSPATH/$FITSNAME")

    # Prepare directories where plots and data will be saved.
    Data_preparation.directory_prep(PATHTOSAVE)


    println("------ CVI CALCULATION ------")
    if DIFFTYPE=="relative" 
        cviallangle,cvimap_averaged = Functionforcvi.construct_cvimap(cvmap,LAG,(DATADIMENSION[1],DATADIMENSION[2]),diff="relative")
        DIFFTYPE = "rel"

    elseif DIFFTYPE=="abs" 
        cviallangle,cvimap_averaged = Functionforcvi.construct_cvimap(cvmap,LAG,(DATADIMENSION[1],DATADIMENSION[2]),diff="absolute")
        DIFFTYPE = "abs"

    else
        error("Not good argument in DIFFTYPE (should be abs or relative)")
    end


    cvimap_averaged = reshape(cvimap_averaged,DATADIMENSION[1],DATADIMENSION[2],size(LAG)[1])
    cvimap_averaged = Data_preparation.replace_missingtoblank(cvimap_averaged,BLANK)
    cvimap_averaged = Data_preparation.blank_equal(cvimap_averaged,0.0,BLANK)
    cvimap_averaged = convert(Array{Float64},cvimap_averaged)

    ## DO NOT NEED TO RECONSTRUCT THE MAP WITH THE BLANK VALUES. ONLY USED TO COMPUTE PDF, SO ONLY Statistics
    cviallangle = Data_preparation.replace_nantoblank(cviallangle,BLANK)
    cviallangle = Data_preparation.replace_missingtoblank(cviallangle,BLANK)
    cviallangle = Data_preparation.blank_equal(cviallangle,0.0,BLANK)
    cviallangle = convert(Array{Float64},cviallangle)

    Data_preparation.write_fits("$(FITSPATH)/$FITSNAME","cvi$(DIFFTYPE)_multlag_METH$(NBPC)","$(PATHTOSAVE)/Data/",cvimap_averaged,(DATADIMENSION[1],DATADIMENSION[2],size(LAG)[1]),BLANK,finished=true,cvi=true,lags=LAG)
    println("CVI map reconstructed from PCA saved in $(PATHTOSAVE)/Data/cvi$(DIFFTYPE)_multlag_$(NBPC)PC_NumberOfFilesWithTheSameNameAsPrefixe.fits as a fits.")

    Data_preparation.write_fits("$(FITSPATH)/$FITSNAME","cvi$(DIFFTYPE)_multlag_allangle_METH$(NBPC)","$(PATHTOSAVE)/Data",cviallangle,(DATADIMENSION[1]*DATADIMENSION[2],NANGLE,size(LAG)[1]),BLANK,finished=true)
    println("CVI map with all angles values reconstructed from PCA saved in the $(PATHTOSAVE)/Data/cvi$(DIFFTYPE)_multlag_allangle_METH$(NBPC)_NumberOfFilesWithTheSameNameAsPrefixe.fits as a fits.")
end #function cvi


"""
    multipca(VARFILEPATH)

Use multiple PCA processes on a cube and with multiple numbers of PCs given as input. A '.txt' file should be used accordingly as an input (see models inside folders '/varfiles/multipca.txt').

OUTPUTS : One file per cube reconstructed with one of the number of PC given as an input.

Use this script in a julia terminal with :
    julia>Unveil.multipca(VARFILEPATH)
"""
function multipca(VARFILEPATH)
    FITSPATH,FILENAME,PATHTOSAVE,UNITVELOCITY,NBPC,BLANK = read_var_files(VARFILEPATH)
    NBPC = [parse(Int, ss) for ss in split(NBPC,",")]

    # Read the fits from the path. Return the data, the VelocityVector, the dimension, the velocity_increment, and the header.
    cube,VELOCITYVECTOR,DATADIMENSION,VELOCITYINCREMENT,HEAD = Data_preparation.read_fits_ppv("$FITSPATH/$FILENAME",UNITVELOCITY ; check=false)


    # Prepare directories where plots and data will be saved.
    Data_preparation.directory_prep(PATHTOSAVE)

    # Replace any NaN value into a missing value and deleted them (can't do PCA on missing values with that package)
    cube = Data_preparation.replace_nantomissing(cube)
    cube = Data_preparation.replace_blanktomissing(cube,BLANK)

    ismis = 0
    if any(ismissing,cube) 
        ismis = 1
        cube,missingplaces1D,missingplaces2D  = Data_preparation.pca_prep(cube,DATADIMENSION)
        cube                                  = convert(Array{Float64},cube)
        DATADIMENSION_NOMISSING               = Data_preparation.read_dim(cube)
    else
        cube                                 = reshape(cube,DATADIMENSION[1]*DATADIMENSION[2],DATADIMENSION[3])
        cube                                 = convert(Array{Float64},cube)
        DATADIMENSION_NOMISSING              = (DATADIMENSION[1]*DATADIMENSION[2],DATADIMENSION[3])
    end


    # Perform the first PCA analysis (same notation as in the MultivariateStats doc)
    println("Perform PCA")
    M, Yt, VARPERCENT,cubereconstructed = Functionforpca.pca(cube,DATADIMENSION[2]-1)

    
    s = open("$(PATHTOSAVE)/Data/Yt_$(NBPC)PC.bin", "w+")
    write(s,Yt)
    close(s)
    Ytpath = "$(PATHTOSAVE)/Data/Yt_$(NBPC)PC.bin"

    for ix=1:size(NBPC)[1]
        cubereconstructed = Functionforpca.pca_nomorecalc(M,Yt,1,NBPC[ix],DATADIMENSION)
        if ismis == 1
            cubereconstructed = Data_preparation.addblank(cubereconstructed,missingplaces2D,BLANK,DATADIMENSION)
        end 


        cubereconstructed = reshape(cubereconstructed,DATADIMENSION)


        println("Saving Fits")
        Data_preparation.write_fits("$FITSPATH/$FILENAME","RECONSTRUCTED_$(NBPC[ix])PC","$PATHTOSAVE/Data",cubereconstructed,DATADIMENSION,BLANK,overwrite=true,more=["NBPC",NBPC[ix],"VARPERCENT",VARPERCENT[NBPC[ix]]])
        println("Data reconstructed from PCA saved in $(PATHTOSAVE)/Data/RECONSTRUCTED_$(NBPC[ix])PC_NumberOfFilesWithTheSameNameAsPrefixe.fits as a fits.")

    end




end #function multipca




"""
    pca(VARFILEPATH)

Use a PCA (Principal Component Analysis) process on a cube and with N PCs given as input. A '.txt' file should be used accordingly as an input (see models inside folders '/varfiles/pca.txt'). 

INPUTS : path to the variable '.txt' file
OUTPUTS : save a cube reconstructed by the number of PC asked in the '.txt' file

Use this function in a julia terminal with :
    julia> Unveil.pca(VARFILEPATH)
"""
function pca(VARFILEPATH)

    #println("")
    #println(" Path to the variable file ? (txt file containing all the informations relevant to read and work on the data)") 
    ##VARFILEPATH = "varfiles/pca.txt"
    #GC.gc()
    FITSPATH,FILENAME,PATHTOSAVE,UNITVELOCITY,NBPC,BLANK = read_var_files(VARFILEPATH)
    (NBPC == 0) && (NBPC="raw")

    # Read the fits from the path. Return the data, the VelocityVector, the dimension, the velocity_increment, and the header.
    cube,VELOCITYVECTOR,DATADIMENSION,VELOCITYINCREMENT,HEAD = Data_preparation.read_fits_ppv("$(FITSPATH)/$(FILENAME)",UNITVELOCITY ; check=true)

    # Prepare directories where plots and data will be saved.
    Data_preparation.directory_prep(PATHTOSAVE)

    # Replace any NaN value into a missing value and deleted them (can't do PCA on missing values with that package)
    cube = Data_preparation.replace_nantomissing(cube)
    cube = Data_preparation.replace_blanktomissing(cube,BLANK)

    ismis = 0
    if any(ismissing,cube) 
        ismis = 1
        cube,missingplaces1D,missingplaces2D  = Data_preparation.pca_prep(cube,DATADIMENSION)
        cube                                  = convert(Array{Float64},cube)
        DATADIMENSION_NOMISSING               = Data_preparation.read_dim(cube)
    else
        cube                                 = reshape(cube,DATADIMENSION[1]*DATADIMENSION[2],DATADIMENSION[3])
        cube                                 = convert(Array{Float64},cube)
        DATADIMENSION_NOMISSING              = (DATADIMENSION[1]*DATADIMENSION[2],DATADIMENSION[3])
    end

    # Important to transform the data into 2D for PCA calculations ; huge gain in time
    #cube = reshape(cube,DATADIMENSION_NOMISSING[1]*DATADIMENSION_NOMISSING[2],DATADIMENSION_NOMISSING[3])

    # Perform the first PCA analysis (same notation as in the MultivariateStats doc)
    println("Perform PCA")
    M, Yt, VARPERCENT,cubereconstructed = Functionforpca.pca(cube,NBPC)

    s = open("$(PATHTOSAVE)/Data/Yt_$(NBPC)PC.bin", "w+")
    write(s,Yt)
    close(s)
    Yt = 0.0
    Ytpath = "$(PATHTOSAVE)/Data/Yt_$(NBPC)PC.bin"

    # Cleaning memory
    cube = 0.0 
    Yt   = 0.0
    head = 0.0
    GC.gc()

    if ismis == 1
        cubereconstructed = Data_preparation.addblank(cubereconstructed,missingplaces2D,BLANK,DATADIMENSION)
    end
    cubereconstructed = reshape(cubereconstructed,DATADIMENSION)

    println("Saving Fits")
    Data_preparation.write_pca_fits(NBPC,VARPERCENT[NBPC],"$(FITSPATH)/$(FILENAME)","RECONSTRUCTED_$(NBPC)PC","$PATHTOSAVE/Data/",cubereconstructed,DATADIMENSION,BLANK)
    println("Data reconstructed from PCA saved in $(PATHTOSAVE)/Data/RECONSTRUCTED_$(NBPC)PC_NumberOfFilesWithTheSameNameAsPrefixe.fits as a fits.")

    cubereconstructed = 0
    GC.gc()

    Graphic.pratio(M,true,NBPC,"Variance reproduction with the number of PC used")
    Plots.savefig("$(PATHTOSAVE)/Convergence_Criteria/TESTpratio_$(NBPC)PC.pdf")
end   #pca



"""
    swo(VARFILEPATH)

Use a SWO (Spectral Window Optimisation) process on a cube. A '.txt' file should be used accordingly as an input (see models inside folders '/varfiles/swo.txt').

Use this function in a julia terminal with :
    julia> Unveil.swo(VARFILEPATH)
"""
function swo(VARFILEPATH)   
    FITSPATH,FILENAME,PATHTOSAVE,UNITVELOCITY,RANGE,BLANK,NOISECANTXT,EXAMPLES = read_var_files(VARFILEPATH)
    NOISECAN = [parse(Int, ss) for ss in split(NOISECANTXT,",")]
    
    
    # Read the fits from the path. Return the data, the VelocityVector, the dimension, the velocity_increment, and the header.
    cube,VELOCITYVECTOR,DATADIMENSION,VELOCITYINCREMENT,HEAD = Data_preparation.read_fits_ppv("$(FITSPATH)/$(FILENAME)",UNITVELOCITY ; check=false)
    
    # Prepare directories where plots and data will be saved.
    Data_preparation.directory_prep(PATHTOSAVE)
    
    # Replace any NaN value into a missing value and deleted them (can't do PCA on missing values with that package).
    cube = Data_preparation.replace_nantomissing(cube)
    ismis = 0
    if any(ismissing,cube) 
        ismis = 1
        cube,missingplaces1D,missingplaces2D  = Data_preparation.pca_prep(cube,DATADIMENSION)
        cube                                  = convert(Array{Float64},cube)
        DATADIMENSION_NOMISSING               = Data_preparation.read_dim(cube)
    else
        cube                                 = reshape(cube,DATADIMENSION[1]*DATADIMENSION[2],DATADIMENSION[3])
        cube                                 = convert(Array{Float64},cube)
        DATADIMENSION_NOMISSING              = (DATADIMENSION[1]*DATADIMENSION[2],DATADIMENSION[3])
    end
    
    
    maskinterv,sigmaT = Spectralwindowopti.optiwind(cube,DATADIMENSION_NOMISSING,VELOCITYVECTOR,NOISECAN,BLANK,RANGE)
    
    if EXAMPLES=="YES"
        Graphic.checkwindowopti(cube,maskinterv,VELOCITYVECTOR,4,4)
        savefig("$(PATHTOSAVE)/checkswo1.pdf")
        Graphic.checkwindowopti(cube,maskinterv,VELOCITYVECTOR,4,4)
        savefig("$(PATHTOSAVE)/checkswo1.pdf")
        Graphic.checkwindowopti(cube,maskinterv,VELOCITYVECTOR,4,4)
        savefig("$(PATHTOSAVE)/checkswo1.pdf")
    end

    if ismis == 1
        maskinterv = Data_preparation.addblank(maskinterv,missingplaces2D,BLANK,DATADIMENSION)
    end
    maskinterv = reshape(maskinterv,DATADIMENSION)

    Data_preparation.write_fits("$(FITSPATH)/$FILENAME","RECONSTRUCTED_$(RANGE)SWO","$PATHTOSAVE/Data/",maskinterv,DATADIMENSION,BLANK,more=["RANGE","$range"])
    println("Data reconstructed from SWO method saved in $(PATHTOSAVE)/Data/RECONSTRUCTED_$(RANGE)SWO_NumberOfFilesWithTheSameNameAsPrefixe.fits as a fits.")



end  #swo






"""
WORK IN PROGRESS
Combine multiple CV maps into one cube. 
Use this script in a julia terminal with :
    julia>Unveil.combinecv()
"""
function combinecv(;VARFILEPATH = "varfiles/multipca.txt")
    cvmap,header,dimens = Data_preparation.read_fits_pp("$FITSPATH/Data/cv_30PC.fits")
    FITSPATH = "/home/delcamps/Data/Pipe_nebula/TestPCA/"
    cubeout = Array{Float64}(undef,size(cvmap)[1],size(cvmap)[2],4)
    count = 0
    for ix in (30,40,50,60)
        global count +=1
        cvmap,header,dimens = Data_preparation.read_fits_pp("$FITSPATH/Data/cv_$(ix)PC.fits")
        cubeout[:,:,count]=cvmap
    end

    Data_preparation.write_fits("$(FITSPATH)/Data/cv_30PC.fits","cv_comb","$(PATHTOSAVE)/Data/",cubeout,(size(cubeout)[1],size(cubeout)[2],size(cubeout)[3]),BLANK,finished=true,more=["NBPC",("30,40,50,60")])

end#function combinecv



end # module Unveil
