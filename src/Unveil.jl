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
include("Functionforpca.jl")   # Calculations of PCA
include("Spectralwindowopti.jl")
include("Graphic.jl")   # Plotting 


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
    pca(VARFILEPATH)

Use a PCA (Principal Component Analysis) process on a cube and with N PCs given as input. A '.txt' file should be used accordingly as an input (see models inside folders '/varfiles/pca.txt').

Use this function in a julia terminal with :
    julia> Unveil.pca(VARFILEPATH)
"""
function pca(;VARFILEPATH = "varfiles/pca.txt")

    println("")
    println(" Path to the variable file ? (txt file containing all the informations relevant to read and work on the data)") 
    #VARFILEPATH = "varfiles/pca.txt"
    GC.gc()
    
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
function swo(;VARFILEPATH="varfiles/swo.txt")
    println("")
    println(" Path to the variable file ? (txt file containing all the informations relevant to read and work on the data)") 
    #VARFILEPATH = "/varfiles/cvoptiwind.txt"
    
    FITSPATH,FILENAME,PATHTOSAVE,UNITVELOCITY,RANGE,BLANK,NOISECANTXT = read_var_files(VARFILEPATH)
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
    
    
    maskinterv,sigmaT = Functionforwindowopti.optiwind(cube,DATADIMENSION_NOMISSING,VELOCITYVECTOR,NOISECAN,BLANK,RANGE)
    println("Saving Fits")
    Data_preparation.write_fits("$(FITSPATH)/$FILENAME","RECONSTRUCTED_$(RANGE)SWO","$PATHTOSAVE/Data/",maskinterv,DATADIMENSION,BLANK,more=["RANGE","$range"])
    println("Data reconstructed from PCA saved in $(PATHTOSAVE)/Data/RECONSTRUCTED_$(RANGE)SWO_NumberOfFilesWithTheSameNameAsPrefixe.fits as a fits.")
end  #swo




"""
    convpca(VARFILEPATH)

Produce calculations to find the PCA convergence criteria based on the matrix projection from PCA. No PCA reconstructed cube will be saved. A '.txt' file should be used accordingly as an input (see models inside folders '/varfiles/convpca.txt').

Use this script in a julia terminal with :
    julia>Unveil.convpca(VARFILEPATH)
"""
function convpca(;VARFILEPATH = "varfiles/convpca.txt")

    println("")
    println(" Path to the variable file ? (txt file containing all the informations relevant to read and work on the data)") 
    #VARFILEPATH = "../varfiles/pcaconvproj.txt"

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

    # Initialize the cube containing all CV maps (one per number of PCs used for the reconstruction)
    #cvcube = Array{Float64}(undef,(DATADIMENSION_NOMISSING[1],HIGHESTPC))


    #First PCA
    println("Perform PCA")
    M, Yt, VARPERCENT,cubereconstructed = Functionforpca.pca(cube,HIGHESTPC)

    proj = projection(M)
    mom1,mom2,mom3,mom4 = Data_analysis.fourmoments(proj,dim=2)
    metric = Data_analysis.calcmetric(mom1,mom2,mom3,mom4,abs(VELOCITYINCREMENT))
    plot(metric)

    if ismis == 1
        proj = Data_preparation.addblank(proj,missingplaces2D[:,1:HIGHESTPC],BLANK,(DATADIMENSION[1],DATADIMENSION[2],HIGHESTPC))
    end
    proj = reshape(proj,(DATADIMENSION[1],DATADIMENSION[2],HIGHESTPC))
    Data_preparation.write_fits("$(FITSPATH)/$FILENAME","projectionmatrix","$(PATHTOSAVE)/Data/",proj,(DATADIMENSION_NOMISSING[1],DATADIMENSION_NOMISSING[2],HIGHESTPC),BLANK,finished=true,overwrite=true)
    xvector = range(1,HIGHESTPC)#[1:HIGHESTPC]
    Graphic.distribcv_multipc(mom1[2:end],mom2[2:end],mom3[2:end],mom4[2:end],metric[2:end],xvector[2:end])
    savefig("$(PATHTOSAVE)/projectionmetric.pdf")

    Data_preparation.write_dat([metric mom1 mom2 mom3 mom4 xvector],"$PATHTOSAVE","metric",overwrite=true,more=["POLARIS_GRAND_CHAMP","Metric  Mom1   Mom2   Mom3   Mom4   PCs"])
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
function convswo(;VARFILEPATH = "varfiles/convswo.txt")
    println("")
    println(" Path to the variable file ? (txt file containing all the informations relevant to read and work on the data)") 
    #VARFILEPATH = "../varfiles/convoptiwind.txt"

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
        cube                                 = reshape(cube,DATADIMENSION[1]*DATADIMENSION[2],DATADIMENSION[3])
        cube                                 = convert(Array{Float64},cube)
        DATADIMENSION_NOMISSING              = (DATADIMENSION[1]*DATADIMENSION[2],DATADIMENSION[3])
    end


    muc,mus,sigc,sigs,gamc,gams,kapc,kaps,metc,mets,SIGMAT=Functionforwindowopti.convoptiwind(cube,DATADIMENSION_NOMISSING,DATADIMENSION,VELOCITYVECTOR,NOISECAN,BLANK,RANGE[1],PATHTOSAVE,missingplaces2D,ismiss=ismis)

    Graphic.distribcv_multiow(muc[1:end],sigc[1:end],gamc[1:end],kapc[1:end],metc[1:end],RANGE[1][1:end],"Parameters of the distribution of the differences \n of Tdv between source cube and window optimized cube",SIGMAT=SIGMAT)
    savefig("$PATHTOSAVE/convoptiwind_cubedif.pdf")

    Graphic.distribcv_multiow(mus[1:end],sigs[1:end],gams[1:end],kaps[1:end],mets[1:end],RANGE[1][1:end],"Parameters of distribution of the percentage of flagged velocity canals",SIGMAT=0)
    savefig("$PATHTOSAVE/convoptiwind_sizenotzero.pdf")

    maskinterv,sigmaT=Functionforwindowopti.optiwind(cube,DATADIMENSION_NOMISSING,VELOCITYVECTOR,NOISECAN,BLANK,4)
    Graphic.checkwindowopti(cube,maskinterv,VELOCITYVECTOR,5,5)
    savefig("$(PATHTOSAVE)/checkwo1.pdf")
    maskinterv,sigmaT=Functionforwindowopti.optiwind(cube,DATADIMENSION_NOMISSING,VELOCITYVECTOR,NOISECAN,BLANK,15)

    Graphic.checkwindowopti(cube,maskinterv,VELOCITYVECTOR,5,5)
    savefig("$(PATHTOSAVE)/checkwo2.pdf")
    maskinterv,sigmaT=Functionforwindowopti.optiwind(cube,DATADIMENSION_NOMISSING,VELOCITYVECTOR,NOISECAN,BLANK,25)

    Graphic.checkwindowopti(cube,maskinterv,VELOCITYVECTOR,5,5)
    savefig("$(PATHTOSAVE)/checkwo3.pdf")

    error()


    if ismis == 1
        maskinterv = Data_preparation.addblank(maskinterv,missingplaces2D,BLANK,DATADIMENSION)
    end

    maskinterv = reshape(maskinterv,(DATADIMENSION[1],DATADIMENSION[2],DATADIMENSION[3]))
    Data_preparation.write_fits("$(FITSPATH)/$(FILENAME)","mask_nopca_optiwind","$(PATHTOSAVE)/Data/",maskinterv,(DATADIMENSION_NOMISSING[1],DATADIMENSION_NOMISSING[2],DATADIMENSION[3]),BLANK,finished=true,overwrite=true)

end #convswo



"""
    cvcvi(VARFILEPATH)

Calculate the CV and CVI from a cube given in input. A '.txt' file should be used accordingly as an input (see models inside folders '/varfiles/cvcvi.txt').

Use this script in a julia terminal with :
    julia>Unveil.cvcvi(VARFILEPATH)
"""
function cvcvi(;VARFILEPATH = "varfiles/cvcvi.txt")

    println("")
    println(" Path to the variable file ? (txt file containing all the informations relevant to read and work on the data)") 
    #VARFILEPATH = "../varfiles/cvi.txt"


    FITSPATH,PATHTOSAVE,UNITVELOCITY,NBPC,BLANK,LAG,NANGLE,DIFFTYPE = read_var_files(VARFILEPATH)
    LAG = [parse(Int, ss) for ss in split(LAG,",")]
    (NBPC == 0) && (NBPC="raw")


    # Read the fits from the path. Return the data, the VelocityVector, the dimension, the velocity_increment, and the header.
    cube,VELOCITYVECTOR,DATADIMENSION,VELOCITYINCREMENT,HEAD = Data_preparation.read_fits_ppv(FITSPATH,UNITVELOCITY ; check=false)

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



    println("------ CVI CALCULATION ------")
    if DIFFTYPE=="relative"
        cvimap_averaged = Functionforcvi.construct_cvimap(cvmap,LAG,NANGLE,(DATADIMENSION[1],DATADIMENSION[2]),diff="relative")
    elseif DIFFTYPE=="abs" 
        cvimap_averaged = Functionforcvi.construct_cvimap(cvmap,LAG,NANGLE,(DATADIMENSION[1],DATADIMENSION[2]),diff="absolute")
    else
        error("Not good argument in DIFFTYPE (should be abs or relative)")
    end


    Data_preparation.write_fits("$(FITSPATH)","cv_$(NBPC)PC","$(PATHTOSAVE)/Data/",cvmap,(DATADIMENSION_NOMISSING[1],DATADIMENSION_NOMISSING[2]),BLANK,finished=true)
    #cvmap = 0.0

    GC.gc()
    println("CV map saved in the $(PATHTOSAVE)/Data/cv_$(NBPC)_NumberOfFilesWithTheSameNameAsPrefixe.fits as a fits.")



    if ismis == 1
        cvimap_averaged = Data_preparation.addblank(cvimap_averaged,missingplaces2D[:,1],BLANK,(DATADIMENSION[1],DATADIMENSION[2]))
        cvimap_averaged = Data_preparation.replace_nantoblank(cvimap_averaged,BLANK)
        cvimap_averaged = Data_preparation.replace_missingtoblank(cvimap_averaged,BLANK)
    end

    cvimap_averaged = reshape(cvimap_averaged,DATADIMENSION[1],DATADIMENSION[2],size(LAG)[1])
    cvimap_averaged = Data_preparation.replace_missingtoblank(cvimap_averaged,BLANK)
    cvimap_averaged = Data_preparation.blank_equal(cvimap_averaged,0.0,BLANK)
    cvimap_averaged = convert(Array{Float64},cvimap_averaged)

    Data_preparation.write_fits("$(FITSPATH)","cvi$(DIFFTYPE)_multlag_$(NBPC)PC","$(PATHTOSAVE)/Data/",cvimap_averaged,(DATADIMENSION_NOMISSING[1],DATADIMENSION_NOMISSING[2],size(LAG)[1]),BLANK,finished=true,cvi=true,lags=LAG)
    println("CVI map reconstructed from PCA saved in the $(PATHTOSAVE)/Data/cvi$(DIFFTYPE)_multlag_$(NBPC)PC_NumberOfFilesWithTheSameNameAsPrefixe.fits as a fits.")


    ## ADDED THE 13/04/22
    ## NEED THE CVI VALUES FOR EACH ANGLE AT A GIVEN LAG
    if DIFFTYPE=="relative"
        cviallangle = Functionforcvi.cv_increment(cvmap,LAG,NANGLE, diff="relative",periodic=false) 
    elseif DIFFTYPE=="abs" 
        cviallangle = Functionforcvi.cv_increment(cvmap,LAG,NANGLE, diff="abs",periodic=false) 
    else
        error("Not good argument in DIFFTYPE (should be abs or relative)")
    end
    ## DO NOT NEED TO RECONSTRUCT THE MAP WITH THE BLANK VALUES. ONLY USED TO COMPUTE PDF, SO ONLY Statistics
    cviallangle = Data_preparation.replace_nantoblank(cviallangle,BLANK)
    cviallangle = Data_preparation.replace_missingtoblank(cviallangle,BLANK)
    cviallangle = Data_preparation.blank_equal(cviallangle,0.0,BLANK)
    cviallangle = convert(Array{Float64},cviallangle)
    Data_preparation.write_fits("$(FITSPATH)","cvi$(DIFFTYPE)_multangle_multlag_$(NBPC)PC","$(PATHTOSAVE)/Data/",cviallangle,(DATADIMENSION_NOMISSING[1]*DATADIMENSION_NOMISSING[2],NANGLE,size(LAG)[1]),BLANK,finished=true)
    println("CVI map with all angles values reconstructed from PCA saved in the $(PATHTOSAVE)/Data/cvi$(DIFFTYPE)_multangle_multlag_$(NBPC)PC_NumberOfFilesWithTheSameNameAsPrefixe.fits as a fits.")

end #function cvcvi




"""
    cvi(VARFILEPATH)

Calculate the CV and CVI from a cube given in input. A '.txt' file should be used accordingly as an input (see models inside folders '/varfiles/cvi.txt').

Use this script in a julia terminal with :
    julia>Unveil.cvi(VARFILEPATH)
"""
function cvi(;VARFILEPATH = "varfiles/cvi.txt")



    println("")
    println(" Path to the variable file ? (txt file containing all the informations relevant to read and work on the data)") 
    #VARFILEPATH = "../varfiles/cvi.txt"


    FITSPATH,PATHTOSAVE,UNITVELOCITY,NBPC,BLANK,LAG,NANGLE,DIFFTYPE = read_var_files(VARFILEPATH)
    LAG = [parse(Int, ss) for ss in split(LAG,",")]
    (NBPC == 0) && (NBPC="raw")


    # Read the fits from the path. Return the data, the VelocityVector, the dimension, the velocity_increment, and the header.
    cvmap,HEAD,DATADIMENSION = Data_preparation.read_fits_pp(FITSPATH)


    println("------ CVI CALCULATION ------")
    if DIFFTYPE=="relative"
        cvimap_averaged = Functionforcvi.construct_cvimap(cvmap,LAG,NANGLE,(DATADIMENSION[1],DATADIMENSION[2]),diff="relative")
    elseif DIFFTYPE=="abs" 
        cvimap_averaged = Functionforcvi.construct_cvimap(cvmap,LAG,NANGLE,(DATADIMENSION[1],DATADIMENSION[2]),diff="absolute")
    else
        error("Not good argument in DIFFTYPE (should be abs or relative)")
    end

    #if ismis == 1
    #    cvimap_averaged = Data_preparation.addblank(cvimap_averaged,missingplaces2D[:,1],BLANK,(DATADIMENSION[1],DATADIMENSION[2]))
    #    cvimap_averaged = Data_preparation.replace_nantoblank(cvimap_averaged,BLANK)
    #    cvimap_averaged = Data_preparation.replace_missingtoblank(cvimap_averaged,BLANK)
    #end

    cvimap_averaged = reshape(cvimap_averaged,DATADIMENSION[1],DATADIMENSION[2],size(LAG)[1])
    cvimap_averaged = Data_preparation.replace_missingtoblank(cvimap_averaged,BLANK)
    cvimap_averaged = Data_preparation.blank_equal(cvimap_averaged,0.0,BLANK)
    cvimap_averaged = convert(Array{Float64},cvimap_averaged)

    Data_preparation.write_fits("$(FITSPATH)","cvi$(DIFFTYPE)_multlag","$(PATHTOSAVE)/Data/",cvimap_averaged,(DATADIMENSION[1],DATADIMENSION[2],size(LAG)[1]),BLANK,finished=true,cvi=true,lags=LAG)
    println("CVI map reconstructed from PCA saved in the $(PATHTOSAVE)/Data/cvi$(DIFFTYPE)_multlag_$(NBPC)PC_NumberOfFilesWithTheSameNameAsPrefixe.fits as a fits.")


    ## ADDED THE 13/04/22
    ## NEED THE CVI VALUES FOR EACH ANGLE AT A GIVEN LAG
    if DIFFTYPE=="relative"
        cviallangle = Functionforcvi.cv_increment(cvmap,LAG,NANGLE, diff="relative",periodic=false) 
    elseif DIFFTYPE=="abs" 
        cviallangle = Functionforcvi.cv_increment(cvmap,LAG,NANGLE, diff="abs",periodic=false) 
    else
        error("Not good argument in DIFFTYPE (should be abs or relative)")
    end
    ## DO NOT NEED TO RECONSTRUCT THE MAP WITH THE BLANK VALUES. ONLY USED TO COMPUTE PDF, SO ONLY Statistics
    cviallangle = Data_preparation.replace_nantoblank(cviallangle,BLANK)
    cviallangle = Data_preparation.replace_missingtoblank(cviallangle,BLANK)
    cviallangle = Data_preparation.blank_equal(cviallangle,0.0,BLANK)
    cviallangle = convert(Array{Float64},cviallangle)
    Data_preparation.write_fits("$(FITSPATH)","cvi$(DIFFTYPE)_multangle_multlag","$(PATHTOSAVE)/Data",cviallangle,(DATADIMENSION[1]*DATADIMENSION[2],NANGLE,size(LAG)[1]),BLANK,finished=true)
    println("CVI map with all angles values reconstructed from PCA saved in the $(PATHTOSAVE)/Data/cvi$(DIFFTYPE)_multangle_multlag_$(NBPC)PC_NumberOfFilesWithTheSameNameAsPrefixe.fits as a fits.")


end #function cvi



"""
    multipca(VARFILEPATH)

Use multiple PCA processes on a cube and with multiple numbers of PCs given as input. A '.txt' file should be used accordingly as an input (see models inside folders '/varfiles/multipca.txt').

Use this script in a julia terminal with :
    julia>Unveil.multipca(VARFILEPATH)
"""
function multipca(;VARFILEPATH = "varfiles/multipca.txt")

    ###################################################################

    ###################################################################

    println("")
    println(" Path to the variable file ? (txt file containing all the informations relevant to read and work on the data)") 
    #VARFILEPATH = "../varfiles/multipca.txt"
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


    # Important to transform the data into 2D for PCA calculations ; huge gain in time
    #cube = reshape(cube,DATADIMENSION_NOMISSING[1]*DATADIMENSION_NOMISSING[2],DATADIMENSION_NOMISSING[3])


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
