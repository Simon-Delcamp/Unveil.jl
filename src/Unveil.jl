
module Unveil


include("Dataprep.jl") # Read and write fits
include("Analysis.jl") 
include("PCA.jl")   # Functions for method PCA
include("SWO.jl")  # Functions for method SWO
include("Graphic.jl")   # Plotting 
include("CVI.jl")   # Functions for CVI computations
include("Structure_functions.jl")

using .Dataprep
using .SWO
using .Graphic
using .Analysis
using .Structure_functions
using .CVI
using .PCA


import Plots 

export pca
export swo
export convpca
export cvcvi
export cvi 
export multipca
export combinecv
export prodvarfile
export structure_functions
export compmethod_stcfct

p = Plots.plot(rand(2,2))
display(p)




"""
    convpca(VARFILEPATH)

Compute the PCA convergence criteria based on the matrix projection from PCA. No PCA reconstructed cube will be saved. A '.txt' file should be used accordingly as an input (see models inside folders '/varfiles/convpca.txt').

OUTPUTS : A plots with every moments of the projection matrix + the metric (see the doc). Also add a .dat file with moments, metric and number of PC used for each.

Use this script in a julia terminal with :
    julia>Unveil.convpca("VARFILEPATH")
"""
function convpca(VARFILEPATH)
    FITSPATH,FILENAME,PATHTOSAVE,SAVENAME,UNITVELOCITY,HIGHESTPC,BLANK,OVERWRITE = Dataprep.read_var_files(VARFILEPATH)

    # Read the fits from the path. Return the data, the VelocityVector, the dimension, the velocity_increment, and the header.
    cube,VELOCITYVECTOR,DATADIMENSION,VELOCITYINCREMENT,HEAD = Dataprep.read_fits_ppv("$(FITSPATH)/$(FILENAME)",UNITVELOCITY ; check=false)

    # Prepare directories where plots and data will be saved.
    Dataprep.directory_prep(PATHTOSAVE)


    # Replace any NaN value into a missing value and deleted them (can't do PCA on missing values with that package).
    cube = Dataprep.replace_nantomissing(cube)
    ismis = 0
    if any(ismissing,cube) 
        ismis = 1
        cube,missingplaces1D,missingplaces2D  = Dataprep.pca_prep(cube,DATADIMENSION)
        cube                                  = convert(Array{Float64},cube)
        DATADIMENSION_NOMISSING               = Dataprep.read_dim(cube)
    else
        cube                                 = reshape(cube,DATADIMENSION[1]*DATADIMENSION[2],DATADIMENSION[3])
        cube                                 = convert(Array{Float64},cube)
        DATADIMENSION_NOMISSING              = (DATADIMENSION[1]*DATADIMENSION[2],DATADIMENSION[3])
    end

    #First PCA
    println("Perform PCA")
    M, Yt, VARPERCENT,cubereconstructed = PCA.pca(cube,HIGHESTPC)

    # Projection matrix
    proj = PCA.proj(M)
    mom1,mom2,mom3,mom4 = Analysis.fourmoments(proj,dim=2)
    metric = Analysis.metricPCA(mom1,mom2,mom3,mom4,abs(VELOCITYINCREMENT))

    if ismis == 1
        proj = Dataprep.addblank(proj,missingplaces2D[:,1:HIGHESTPC],BLANK,(DATADIMENSION[1],DATADIMENSION[2],HIGHESTPC))
    end
    proj = reshape(proj,(DATADIMENSION[1],DATADIMENSION[2],HIGHESTPC))
    Dataprep.write_fits("$(FITSPATH)/$FILENAME","$(SAVENAME)_projectionmatrix","$(PATHTOSAVE)/Data/",proj,(DATADIMENSION_NOMISSING[1],DATADIMENSION_NOMISSING[2],HIGHESTPC),BLANK,finished=true,overwrite=OVERWRITE)
    xvector = range(1,HIGHESTPC)#[1:HIGHESTPC]

    Graphic.distribcv_multipc(mom1[2:end-1],mom2[2:end-1],mom3[2:end-1],mom4[2:end-1],metric[2:end-1],xvector[2:end-1])
    Plots.savefig("$(PATHTOSAVE)/Figures/pca_$(SAVENAME)_metric.pdf")

    Dataprep.write_dat([metric mom1 mom2 mom3 mom4 xvector],"$PATHTOSAVE","$(SAVENAME)_metricPCA",overwrite=true,more=["$FILENAME","Metric  Mom1   Mom2   Mom3   Mom4   PCs"])
    minimetr = minimum(metric)
    minipc   = xvector[findall(x->x==minimetr,metric)]
    println("Metric minimum=$(minimetr)")
    println("#PC of metric minimum=$(minipc)")
end #convpca




"""
    cv(VARFILEPATH)

Calculate the CV from a cube given in input. A '.txt' file should be used accordingly as an input (see models inside folders '/varfiles/cv.txt').

Use this script in a julia terminal with :
    julia>Unveil.cv(VARFILEPATH)
"""
function cv(VARFILEPATH)
    FITSPATH,FITSNAME,PATHTOSAVE,SAVENAME,UNITVELOCITY,THRESHOLD,NOISECANTXT,BLANK,OVERWRITE = Dataprep.read_var_files(VARFILEPATH)
    NOISECAN = [parse(Int, ss) for ss in split(NOISECANTXT,",")]


    # Read the fits from the path. Return the data, the VelocityVector, the dimension, the velocity_increment, and the header.
    cube,VELOCITYVECTOR,DATADIMENSION,VELOCITYINCREMENT,HEAD = Dataprep.read_fits_ppv("$FITSPATH/$FITSNAME",UNITVELOCITY ; check=false)
    if haskey(HEAD,"METHOD")==1 
        if HEAD["METHOD"]=="PCA"
            METH = HEAD["NBPC"]
            METH = "$(METH)PC"
        elseif HEAD["METHOD"]=="SWO"
            METH = "SWO"
        end
    else
        METH = "raw"
    end


    # Prepare directories where plots and data will be saved.
    Dataprep.directory_prep(PATHTOSAVE)

    # If no threshold, change it to the blank value and SIGMAT to 1, because function moment_one_field will blank every values lower than SIGMAT*THRESHOLD
    if THRESHOLD==0
        THRESHOLD = BLANK
        SIGMAT = 1
    else
        SIGMAT = Analysis.rms_cube(cube,NOISECAN)[2]
    end

    # Replace any NaN value into a missing value and deleted them (can't do PCA on missing values with that package).
    cube = Dataprep.replace_nantomissing(cube)

    ismis = 0
    if any(ismissing,cube) 
        ismis = 1
        cube,missingplaces1D,missingplaces2D  = Dataprep.pca_prep(cube,DATADIMENSION)
        cube                                  = convert(Array{Float64},cube)
        DATADIMENSION_NOMISSING               = Dataprep.read_dim(cube)
    else
        cube                                 = reshape(cube,DATADIMENSION[1]*DATADIMENSION[2],DATADIMENSION[3])
        cube                                 = convert(Array{Float64},cube)
        DATADIMENSION_NOMISSING              = (DATADIMENSION[1]*DATADIMENSION[2],DATADIMENSION[3])
    end

    println("------ CV CALCULATION ------")
    cvmap = CVI.moment_one_field(cube,SIGMAT,THRESHOLD,VELOCITYVECTOR,BLANK) # Calculate the first velocity moment order on data reconstructed
    cube  = 0.0 
    GC.gc()

    if ismis == 1
        cvmap = Dataprep.addblank(cvmap,missingplaces2D[:,1],BLANK,(DATADIMENSION[1],DATADIMENSION[2]))
    end
    cvmap = reshape(cvmap,(DATADIMENSION[1],DATADIMENSION[2]))


    Dataprep.write_fits("$(FITSPATH)/$FITSNAME","CV_$(SAVENAME)_$(METH)","$(PATHTOSAVE)/Data/",cvmap,(DATADIMENSION_NOMISSING[1],DATADIMENSION_NOMISSING[2]),BLANK,finished=true,overwrite=OVERWRITE)
    #cvmap = 0.0

    println("CV map saved in $(PATHTOSAVE)/Data/CV_$(SAVENAME)_$(METH)_NumberOfFilesWithTheSameNameAsPrefixe.fits as a fits.")


end #cv


"""
    cvcvi(VARFILEPATH)

Calculate the CV and CVI from a cube given in input. A '.txt' file should be used accordingly as an input (see models inside folders '/varfiles/cvcvi.txt').

Use this script in a julia terminal with :
    julia>Unveil.cvcvi(VARFILEPATH)
"""
function cvcvi(VARFILEPATH)
    FITSPATH,FITSNAME,PATHTOSAVE,SAVENAME,THRESHOLD,NOISECANTXT,UNITVELOCITY,REMOVE,BLANK,LAG,DIFFTYPE,OVERWRITE = Dataprep.read_var_files(VARFILEPATH)
    NOISECAN = [parse(Int, ss) for ss in split(NOISECANTXT,",")]
    LAG = [parse(Int, ss) for ss in split(LAG,",")]


    # Read the fits from the path. Return the data, the VelocityVector, the dimension, the velocity_increment, and the header.
    cube,VELOCITYVECTOR,DATADIMENSION,VELOCITYINCREMENT,HEAD = Dataprep.read_fits_ppv("$FITSPATH/$FITSNAME",UNITVELOCITY ; check=false)
    if haskey(HEAD,"METHOD")==1 
        if HEAD["METHOD"]=="PCA"
            METH = HEAD["NBPC"]
            METH = "$(METH)PC"
        elseif HEAD["METHOD"]=="SWO"

            METH = "SWO"
        end
    else
        METH = "raw"
    end


    # Prepare directories where plots and data will be saved.
    Dataprep.directory_prep(PATHTOSAVE)

    # If no threshold, change it to the blank value and SIGMAT to 1, because function moment_one_field will blank every values lower than SIGMAT*THRESHOLD
    if THRESHOLD==0
        THRESHOLD = BLANK
        SIGMAT = 1
        
    else
        SIGMAT = Analysis.rms_cube(cube,NOISECAN)[2]
    end

    # Replace any NaN value into a missing value and deleted them (can't do PCA on missing values with that package).

    SIGMAMAP= Analysis.rms_cube(cube,NOISECAN)[1]
    cube = Dataprep.replace_nantomissing(cube)
    cube = Dataprep.replace_blanktomissing(cube,BLANK)
    
    if REMOVE==true
        cube = Dataprep.replace_nosignal(cube,DATADIMENSION,VELOCITYVECTOR,BLANK,SIGMAMAP)
    end

    ismis = 0
    if any(ismissing,cube) 
        ismis = 1
        cube,missingplaces1D,missingplaces2D  = Dataprep.pca_prep(cube,DATADIMENSION)
        cube                                  = convert(Array{Float64},cube)
        DATADIMENSION_NOMISSING               = Dataprep.read_dim(cube)
    else
        cube                                 = reshape(cube,DATADIMENSION[1]*DATADIMENSION[2],DATADIMENSION[3])
        cube                                 = convert(Array{Float64},cube)
        DATADIMENSION_NOMISSING              = (DATADIMENSION[1]*DATADIMENSION[2],DATADIMENSION[3])
    end


    println("------ CV CALCULATION ------")
    cvmap = CVI.moment_one_field(cube,SIGMAT,THRESHOLD,VELOCITYVECTOR,BLANK) # Calculate the first velocity moment order on data reconstructed
    cube  = 0.0 
    GC.gc()

    if ismis == 1
        cvmap = Dataprep.addblank(cvmap,missingplaces2D[:,1],BLANK,(DATADIMENSION[1],DATADIMENSION[2]))
    end
    cvmap = reshape(cvmap,(DATADIMENSION[1],DATADIMENSION[2]))


    Dataprep.write_fits("$(FITSPATH)/$FITSNAME","CV_$(SAVENAME)_$(METH)","$(PATHTOSAVE)/Data/",cvmap,(DATADIMENSION_NOMISSING[1],DATADIMENSION_NOMISSING[2]),BLANK,finished=true,overwrite=OVERWRITE)
    #cvmap = 0.0

    GC.gc()
    println("CV map saved in $(PATHTOSAVE)/Data/CV_$(SAVENAME)_$(METH)_NumberOfFilesWithTheSameNameAsPrefixe.fits as a fits.")

    cvmap = Dataprep.replace_nantomissing(cvmap)
    cvmap = Dataprep.replace_blanktomissing(cvmap,BLANK)
    println("------ CVI CALCULATION ------")
    if DIFFTYPE=="relative" 
        cviallangle,cvimap_averaged,NANGLE = CVI.construct_cvimap(cvmap,LAG,(DATADIMENSION[1],DATADIMENSION[2]),diff="relative")
        DIFFTYPE = "rel"

    elseif DIFFTYPE=="abs" 
        cviallangle,cvimap_averaged,NANGLE = CVI.construct_cvimap(cvmap,LAG,(DATADIMENSION[1],DATADIMENSION[2]),diff="absolute")
        DIFFTYPE = "abs"

    else
        error("Not good argument in DIFFTYPE (should be abs or relative)")
    end


    cvimap_averaged = reshape(cvimap_averaged,DATADIMENSION[1],DATADIMENSION[2],size(LAG)[1])
    for lag=1:size(LAG)[1]
        cvimap_averaged[:,1:LAG[lag],lag] .= missing
        cvimap_averaged[1:LAG[lag],:,lag] .= missing
        cvimap_averaged[(DATADIMENSION[1]-LAG[lag])+1:DATADIMENSION[1],:,lag] .= missing
        cvimap_averaged[:,(DATADIMENSION[2]-LAG[lag])+1:DATADIMENSION[2],lag] .= missing
    end
    cvimap_averaged = Dataprep.replace_missingtoblank(cvimap_averaged,BLANK)
    cvimap_averaged = Dataprep.blank_equal(cvimap_averaged,0.0,BLANK)
    cvimap_averaged = convert(Array{Float64},cvimap_averaged)

    ## DO NOT NEED TO RECONSTRUCT THE MAP WITH THE BLANK VALUES. ONLY USED TO COMPUTE PDF, SO ONLY Statistics
    cviallangle = reshape(cviallangle,DATADIMENSION[1],DATADIMENSION[2],maximum(NANGLE),size(LAG)[1])
    for lag=1:size(LAG)[1]
        cviallangle[:,1:LAG[lag],:,lag] .= missing
        cviallangle[1:LAG[lag],:,:,lag] .= missing
        cviallangle[(DATADIMENSION[1]-LAG[lag])+1:DATADIMENSION[1],:,:,lag] .= missing
        cviallangle[:,(DATADIMENSION[2]-LAG[lag])+1:DATADIMENSION[2],:,lag] .= missing
    end
    cviallangle = reshape(cviallangle,DATADIMENSION[1]*DATADIMENSION[2],maximum(NANGLE),size(LAG)[1])

    cviallangle = Dataprep.replace_nantoblank(cviallangle,BLANK)
    cviallangle = Dataprep.replace_missingtoblank(cviallangle,BLANK)
    cviallangle = Dataprep.blank_equal(cviallangle,0.0,BLANK)
    cviallangle = convert(Array{Float64},cviallangle)


    Dataprep.write_fits("$(FITSPATH)/$FITSNAME","CVI$(DIFFTYPE)_$(SAVENAME)_$(METH)","$(PATHTOSAVE)/Data/",cvimap_averaged,(DATADIMENSION[1],DATADIMENSION[2],size(LAG)[1]),BLANK,finished=true,overwrite=OVERWRITE,more=["LAG",LAG])
    println("CVI map saved in $(PATHTOSAVE)/Data/CVI$(DIFFTYPE)_$(SAVENAME)_$(METH)_NumberOfFilesWithTheSameNameAsPrefixe.fits as a fits.")

    Dataprep.write_fits("$(FITSPATH)/$FITSNAME","CVI$(DIFFTYPE)_$(SAVENAME)_allangle_$(METH)","$(PATHTOSAVE)/Data",cviallangle,(DATADIMENSION[1]*DATADIMENSION[2],maximum(NANGLE),size(LAG)[1]),BLANK,finished=true,overwrite=OVERWRITE,more=["LAG",LAG])
    println("CVI map with all angles values saved in the $(PATHTOSAVE)/Data/CVI$(DIFFTYPE)_$(SAVENAME)_allangle_$(METH)_NumberOfFilesWithTheSameNameAsPrefixe.fits as a fits.")

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
    FITSPATH,FITSNAME,PATHTOSAVE,SAVENAME,BLANK,LAG,DIFFTYPE,OVERWRITE = Dataprep.read_var_files(VARFILEPATH)
    LAG = [parse(Int, ss) for ss in split(LAG,",")]

    cvmap,HEAD,DATADIMENSION = Dataprep.read_fits_pp("$FITSPATH/$FITSNAME")
    if haskey(HEAD,"METHOD")==1 
        if HEAD["METHOD"]=="PCA"
            METH = HEAD["NBPC"]
            METH = "$(METH)PC"
        elseif HEAD["METHOD"]=="SWO"
            METH = "SWO"
        end
    else
        METH = "raw"
    end

    # Prepare directories where plots and data will be saved.
    Dataprep.directory_prep(PATHTOSAVE)

    cvmap = Dataprep.replace_nantomissing(cvmap)   
    cvmap = Dataprep.replace_blanktomissing(cvmap,BLANK)


    println("------ CVI CALCULATION ------")
    if DIFFTYPE=="relative" 
        cviallangle,cvimap_averaged,NANGLE = CVI.construct_cvimap(cvmap,LAG,(DATADIMENSION[1],DATADIMENSION[2]),diff="relative",BLANK=BLANK)
        DIFFTYPE = "rel"

    elseif DIFFTYPE=="abs" 
        cviallangle,cvimap_averaged,NANGLE = CVI.construct_cvimap(cvmap,LAG,(DATADIMENSION[1],DATADIMENSION[2]),diff="absolute",BLANK=BLANK)
        DIFFTYPE = "abs"

    else
        error("Not good argument in DIFFTYPE (should be abs or relative)")
    end


    cvimap_averaged = reshape(cvimap_averaged,DATADIMENSION[1],DATADIMENSION[2],size(LAG)[1])
    for lag=1:size(LAG)[1]
        cvimap_averaged[:,1:LAG[lag],lag] .= missing
        cvimap_averaged[1:LAG[lag],:,lag] .= missing
        cvimap_averaged[(DATADIMENSION[1]-LAG[lag])+1:DATADIMENSION[1],:,lag] .= missing
        cvimap_averaged[:,(DATADIMENSION[2]-LAG[lag])+1:DATADIMENSION[2],lag] .= missing
    end

    cvimap_averaged = Dataprep.replace_missingtoblank(cvimap_averaged,BLANK)
    cvimap_averaged = Dataprep.blank_equal(cvimap_averaged,0.0,BLANK)
    cvimap_averaged = convert(Array{Float64},cvimap_averaged)

    ## DO NOT NEED TO RECONSTRUCT THE MAP WITH THE BLANK VALUES. ONLY USED TO COMPUTE PDF, SO ONLY Statistics
    cviallangle = reshape(cviallangle,DATADIMENSION[1],DATADIMENSION[2],maximum(NANGLE),size(LAG)[1])
    for lag=1:size(LAG)[1]
        cviallangle[:,1:LAG[lag],:,lag] .= missing
        cviallangle[1:LAG[lag],:,:,lag] .= missing
        cviallangle[(DATADIMENSION[1]-LAG[lag])+1:DATADIMENSION[1],:,:,lag] .= missing
        cviallangle[:,(DATADIMENSION[2]-LAG[lag])+1:DATADIMENSION[2],:,lag] .= missing
    end
    cviallangle = reshape(cviallangle,DATADIMENSION[1]*DATADIMENSION[2],maximum(NANGLE),size(LAG)[1])

    cviallangle = Dataprep.replace_nantoblank(cviallangle,BLANK)
    cviallangle = Dataprep.replace_missingtoblank(cviallangle,BLANK)
    cviallangle = Dataprep.blank_equal(cviallangle,0.0,BLANK)
    cviallangle = convert(Array{Float64},cviallangle)

    Dataprep.write_fits("$(FITSPATH)/$FITSNAME","CVI$(DIFFTYPE)_$(SAVENAME)_$(METH)","$(PATHTOSAVE)/Data/",cvimap_averaged,(DATADIMENSION[1],DATADIMENSION[2],size(LAG)[1]),BLANK,finished=true,overwrite=OVERWRITE,more=["LAG",LAG])
    println("CVI map saved in $(PATHTOSAVE)/Data/CVI$(DIFFTYPE)_$(SAVENAME)_$(METH)_NumberOfFilesWithTheSameNameAsPrefixe.fits as a fits.")


    Dataprep.write_fits("$(FITSPATH)/$FITSNAME","CVI$(DIFFTYPE)_$(SAVENAME)_allangle_$(METH)","$(PATHTOSAVE)/Data",cviallangle,(DATADIMENSION[1]*DATADIMENSION[2],maximum(NANGLE),size(LAG)[1]),BLANK,finished=true,overwrite=OVERWRITE,more=["LAG",LAG])
    println("CVI map with all angles values saved in the $(PATHTOSAVE)/Data/CVI$(DIFFTYPE)_$(SAVENAME)_allangle_$(METH)_NumberOfFilesWithTheSameNameAsPrefixe.fits as a fits.")

end #function cvi


"""
    multipca(VARFILEPATH)

Use multiple PCA processes on a cube and with multiple numbers of PCs given as input. A '.txt' file should be used accordingly as an input (see models inside folders '/varfiles/multipca.txt').

OUTPUTS : One file per cube reconstructed with one of the number of PC given as an input.

Use this script in a julia terminal with :
    julia>Unveil.multipca(VARFILEPATH)
"""
function multipca(VARFILEPATH)
    FITSPATH,FILENAME,PATHTOSAVE,SAVENAME,UNITVELOCITY,NBPC,BLANK,OVERWRITE = Dataprep.read_var_files(VARFILEPATH)
    NBPC = [parse(Int, ss) for ss in split(NBPC,",")]

    # Read the fits from the path. Return the data, the VelocityVector, the dimension, the velocity_increment, and the header.
    cube,VELOCITYVECTOR,DATADIMENSION,VELOCITYINCREMENT,HEAD = Dataprep.read_fits_ppv("$FITSPATH/$FILENAME",UNITVELOCITY ; check=false)


    # Prepare directories where plots and data will be saved.
    Dataprep.directory_prep(PATHTOSAVE)

    # Replace any NaN value into a missing value and deleted them (can't do PCA on missing values with that package)
    cube = Dataprep.replace_nantomissing(cube)
    cube = Dataprep.replace_blanktomissing(cube,BLANK)

    ismis = 0
    if any(ismissing,cube) 
        ismis = 1
        cube,missingplaces1D,missingplaces2D  = Dataprep.pca_prep(cube,DATADIMENSION)
        cube                                  = convert(Array{Float64},cube)
        DATADIMENSION_NOMISSING               = Dataprep.read_dim(cube)
    else
        cube                                 = reshape(cube,DATADIMENSION[1]*DATADIMENSION[2],DATADIMENSION[3])
        cube                                 = convert(Array{Float64},cube)
        DATADIMENSION_NOMISSING              = (DATADIMENSION[1]*DATADIMENSION[2],DATADIMENSION[3])
    end


    # Perform the first PCA analysis (same notation as in the MultivariateStats doc)
    println("Perform PCA")
    M, Yt, VARPERCENT,cubereconstructed = PCA.pca(cube,DATADIMENSION[2]-1)

    
    s = open("$(PATHTOSAVE)/Data/Yt_$(NBPC)PC.bin", "w+")
    write(s,Yt)
    close(s)
    Ytpath = "$(PATHTOSAVE)/Data/Yt_$(NBPC)PC.bin"

    for ix=1:size(NBPC)[1]
        cubereconstructed = PCA.pca_nomorecalc(M,Yt,1,NBPC[ix],DATADIMENSION)
        if ismis == 1
            cubereconstructed = Dataprep.addblank(cubereconstructed,missingplaces2D,BLANK,DATADIMENSION)
        end 


        cubereconstructed = reshape(cubereconstructed,DATADIMENSION)


        println("Saving Fits")
        Dataprep.write_fits("$FITSPATH/$FILENAME","RECONSTRUCTED_$(SAVENAME)_$(NBPC[ix])PC","$PATHTOSAVE/Data",cubereconstructed,DATADIMENSION,BLANK,more=["NBPC",NBPC[ix],"VARPERCENT",VARPERCENT[NBPC[ix]]],overwrite=OVERWRITE)
        println("Data reconstructed from PCA saved in $(PATHTOSAVE)/Data/RECONSTRUCTED_$(SAVENAME)_$(NBPC[ix])PC_NumberOfFilesWithTheSameNameAsPrefixe.fits as a fits.")

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
    FITSPATH,FILENAME,PATHTOSAVE,SAVENAME,UNITVELOCITY,NBPC,BLANK,NOISECANTXT,OVERWRITE = Dataprep.read_var_files(VARFILEPATH)
    (NBPC == 0) && (NBPC="raw")
    NOISECAN = [parse(Int, ss) for ss in split(NOISECANTXT,",")]

    # Read the fits from the path. Return the data, the VelocityVector, the dimension, the velocity_increment, and the header.
    cube,VELOCITYVECTOR,DATADIMENSION,VELOCITYINCREMENT,HEAD = Dataprep.read_fits_ppv("$(FITSPATH)/$(FILENAME)",UNITVELOCITY ; check=false)

    # Prepare directories where plots and data will be saved.
    Dataprep.directory_prep(PATHTOSAVE)

    # Replace any NaN value into a missing value and deleted them (can't do PCA on missing values with that package)
    SIGMAMAP= Analysis.rms_cube(cube,NOISECAN)[1]
    cube = Dataprep.replace_nantomissing(cube)
    cube = Dataprep.replace_blanktomissing(cube,BLANK)
    
    cube = Dataprep.replace_nosignal(cube,DATADIMENSION,VELOCITYVECTOR,BLANK,SIGMAMAP)


    ismis = 0
    if any(ismissing,cube) 
        ismis = 1
        cube,missingplaces1D,missingplaces2D  = Dataprep.pca_prep(cube,DATADIMENSION)
        cube                                  = convert(Array{Float64},cube)
        DATADIMENSION_NOMISSING               = Dataprep.read_dim(cube)
    else
        cube                                 = reshape(cube,DATADIMENSION[1]*DATADIMENSION[2],DATADIMENSION[3])
        cube                                 = convert(Array{Float64},cube)
        DATADIMENSION_NOMISSING              = (DATADIMENSION[1]*DATADIMENSION[2],DATADIMENSION[3])
    end

    # Important to transform the data into 2D for PCA calculations ; huge gain in time
    #cube = reshape(cube,DATADIMENSION_NOMISSING[1]*DATADIMENSION_NOMISSING[2],DATADIMENSION_NOMISSING[3])

    # Perform the first PCA analysis (same notation as in the MultivariateStats doc)
    println("Perform PCA")
    M, Yt, VARPERCENT,cubereconstructed = PCA.pca(cube,NBPC)
    println(VARPERCENT[NBPC])

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
        cubereconstructed = Dataprep.addblank(cubereconstructed,missingplaces2D,BLANK,DATADIMENSION)
        projec            = Dataprep.addblank(PCA.proj(M),missingplaces2D[:,1],BLANK,DATADIMENSION) 
    end
    cubereconstructed = reshape(cubereconstructed,DATADIMENSION)
    projec            = reshape(projec,DATADIMENSION)

    println("Saving Fits")
    Dataprep.write_fits("$(FITSPATH)/$(FILENAME)","RECONSTRUCTED_$(SAVENAME)_$(NBPC)PC","$PATHTOSAVE/Data/",cubereconstructed,DATADIMENSION,BLANK,overwrite=OVERWRITE,more=["NBPC",NBPC,"VARPERC",VARPERCENT[NBPC]*100,"METHOD","PCA"])
    println("Data reconstructed from PCA saved in $(PATHTOSAVE)/Data/RECONSTRUCTED_$(SAVENAME)_$(NBPC)PC_NumberOfFilesWithTheSameNameAsPrefixe.fits as a fits.")

    Dataprep.write_fits("$(FITSPATH)/$(FILENAME)","PROJMATRIX_$(SAVENAME)_$(NBPC)PC","$PATHTOSAVE/Data/",projec,DATADIMENSION,BLANK,overwrite=OVERWRITE,more=["NBPC",NBPC,"VARPERC",VARPERCENT[NBPC]*100,"METHOD","PCA"])
    println("Projection matrix from PCA saved in $(PATHTOSAVE)/Data/PROJMATRIX_$(SAVENAME)_$(NBPC)PC_NumberOfFilesWithTheSameNameAsPrefixe.fits as a fits.")
    cubereconstructed = 0
    GC.gc()

    #Graphic.pratio(M,true,NBPC,"Variance reproduction with the number of PC used")
    #Plots.savefig("$(PATHTOSAVE)/Convergence_Criteria/TESTpratio_$(NBPC)PC.pdf")
end   #pca



"""
    swo(VARFILEPATH)

Use a SWO (Spectral Window Optimisation) process on a cube. A '.txt' file should be used accordingly as an input (see models inside folders '/varfiles/swo.txt').

Use this function in a julia terminal with :
    julia> Unveil.swo(VARFILEPATH)

"""
function swo(VARFILEPATH)   
    FITSPATH,FILENAME,PATHTOSAVE,SAVENAME,UNITVELOCITY,BLANK,NOISECANTXT,EXAMPLES,OVERWRITE = Dataprep.read_var_files(VARFILEPATH)
    NOISECAN = [parse(Int, ss) for ss in split(NOISECANTXT,",")]

    # Read the fits from the path. Return the data, the VelocityVector, the dimension, the velocity_increment, and the header.
    cube,VELOCITYVECTOR,DATADIMENSION,VELOCITYINCREMENT,HEAD = Dataprep.read_fits_ppv("$(FITSPATH)/$(FILENAME)",UNITVELOCITY ; check=false)
    
    # Prepare directories where plots and data will be saved.
    Dataprep.directory_prep(PATHTOSAVE)
    
    # Replace any NaN value into a missing value and deleted them (can't do PCA on missing values with that package).
    cube = Dataprep.replace_nantomissing(cube)
    ismis = 0
    if any(ismissing,cube) 
        ismis = 1
        cube,missingplaces1D,missingplaces2D  = Dataprep.pca_prep(cube,DATADIMENSION)
        cube                                  = convert(Array{Float64},cube)
        DATADIMENSION_NOMISSING               = Dataprep.read_dim(cube)
    else
        cube                                 = reshape(cube,DATADIMENSION[1]*DATADIMENSION[2],DATADIMENSION[3])
        cube                                 = convert(Array{Float64},cube)
        DATADIMENSION_NOMISSING              = (DATADIMENSION[1]*DATADIMENSION[2],DATADIMENSION[3])
    end
    
    maskinterv,mask = SWO.bestsnr(cube,DATADIMENSION_NOMISSING,VELOCITYVECTOR,NOISECAN)

    if EXAMPLES=="YES"
        Graphic.checkwindowopti(cube,maskinterv,mask,VELOCITYVECTOR,6,6)
        Plots.savefig("$(PATHTOSAVE)/Figures/checkswo1.pdf")

        Graphic.checkwindowopti(cube,maskinterv,mask,VELOCITYVECTOR,6,6)
        Plots.savefig("$(PATHTOSAVE)/Figures/checkswo2.pdf")

        Graphic.checkwindowopti(cube,maskinterv,mask,VELOCITYVECTOR,6,6)
        Plots.savefig("$(PATHTOSAVE)/Figures/checkswo3.pdf")
    end

    if ismis == 1
        maskinterv = Dataprep.addblank(maskinterv,missingplaces2D,BLANK,DATADIMENSION)
    end
    maskinterv = reshape(maskinterv,DATADIMENSION)

    maskinterv = Dataprep.blank_equal(maskinterv,BLANK,0)

    Dataprep.write_fits("$(FITSPATH)/$FILENAME","RECONSTRUCTED_$(SAVENAME)_SWO","$PATHTOSAVE/Data/",maskinterv,DATADIMENSION,BLANK,overwrite=OVERWRITE,more=["METHOD","SWO"])
    println("Data reconstructed from SWO method saved in $(PATHTOSAVE)/Data/RECONSTRUCTED_$(SAVENAME)_SWO_NumberOfFilesWithTheSameNameAsPrefixe.fits as a fits.")
end




"""
    structure_functions(VARFILEPATH)

Compute the structure functions ``S_p(l)`` of a cvi cube (following definition in Kritsuk+2007). A '.txt' file should be used accordingly as an input (see models inside folders '/varfiles/structure_functions.txt'). The cube given as input should be a CVI (use function **Unveil.cvi** if needed). Prefer a cube with rotations of every lags than azimutal average, like that : (Pixel positions,angles,lag). 

OUTPUTS : One figure with ``S_p(l)`` vs ``S_3(l)``, same with fit, and exponants of ``S_p(l)`` function of p. Also a *.dat* file with values of the exponants.
Use this function in a julia terminal with :
    julia> Unveil.structure_functions(VARFILEPATH)
"""
function structure_functions(VARFILEPATH)
    FITSPATH,FILENAME,PATHTOSAVE,SAVENAME,ORDERSTXT,BLANK,OVERWRITE = Dataprep.read_var_files(VARFILEPATH)
    ORDERS = [parse(Int, ss) for ss in split(ORDERSTXT,",")]


    # Read the fits from the path. Return the data, the VelocityVector, the dimension, the velocity_increment, and the header.
    cvicube,DATADIMENSION,HEAD = Dataprep.read_fits_cvi("$(FITSPATH)/$(FILENAME)" ; check=false)
    if haskey(HEAD,"METHOD")==1 
        if HEAD["METHOD"]=="PCA"
            METH = HEAD["NBPC"]
            METH = "$(METH)PC"
            METHV = HEAD["NBPC"]
        elseif HEAD["METHOD"]=="SWO"
            METH = "SWO"
            METHV = -1
        end
    else
        METHV = 0
        METH = "raw"
    end

    LAG = [parse(Int,ss) for ss in split(HEAD["LAG"][2:end-1],",")]

    # Prepare directories where plots and data will be saved.
    Dataprep.directory_prep(PATHTOSAVE)

    cvicube = Dataprep.replace_nantomissing(cvicube)
    cvicube = Dataprep.replace_blanktomissing(cvicube,BLANK)
    #cvicube = Dataprep.replace_blanktomissing(cvicube,0)


    sct = Structure_functions.fct_sct(cvicube,LAG,ORDERS)  # order,lag
    Graphic.StcFct(sct,sct[3,:],ORDERS,"$SAVENAME")
    if OVERWRITE==true
        Plots.savefig("$(PATHTOSAVE)/Figures/stc_fct_$(SAVENAME)_$(METH).pdf")
    elseif (OVERWRITE==false && isfile("$(PATHTOSAVE)/Figures/stc_fct_$(SAVENAME)_$(METH).pdf")==true)
        println("THE GIVEN FILE NAME ALREADY EXIST. AN INDICE WILL BE ADDED AT THE END OF THE GIVEN NAME, EQUAL TO THE NUMBER OF FILES WITH THE SAME NAME +1 ")
        count = 1
        for ix=1:size((findall.("stc_fct_$(SAVENAME)_$(METH).pdf",readdir("$(PATHTOSAVE)/Figures/"))))[1]
            if size(findall("stc_fct_$(SAVENAME)_$(METH).pdf",readdir("$(PATHTOSAVE)/Figures/")[ix]))[1]!=0
                count += 1
            end
        end
        newname = "stc_fct_$(SAVENAME)_$(METH)"
        Plots.savefig("$(PATHTOSAVE)/Figures/$(newname).pdf")
    end

    println(" ")
    println("On which Lag to fit ? Give the indices in the array Lag of the var file (first indice=1). Size of the array : $(size(LAG)[1])")
    println("First indice : ")
    CANALINF   = parse(Int64,readline())
    println("Second indice : ")
    CANALSUP   = parse(Int64,readline())
    CANALTOFIT = CANALINF:CANALSUP

    zeta = Structure_functions.xhi_fct_p(ORDERS[:],sct[:,CANALTOFIT])
    Graphic.StcFctWithFit(sct,sct[3,:],ORDERS,zeta,CANALTOFIT,"$SAVENAME")
    if OVERWRITE==true
        Plots.savefig("$(PATHTOSAVE)/Figures/stc_fctfit_$(SAVENAME)_$(METH).pdf")
    elseif (OVERWRITE==false && isfile("$(PATHTOSAVE)/Figures/stc_fctfit_$(SAVENAME)_$(METH).pdf")==true)
        println("THE GIVEN FILE NAME ALREADY EXIST. AN INDICE WILL BE ADDED AT THE END OF THE GIVEN NAME, EQUAL TO THE NUMBER OF FILES WITH THE SAME NAME +1 ")
        count = 1
        for ix=1:size((findall.("stc_fctfit_$(SAVENAME)_$(METH).pdf",readdir("$(PATHTOSAVE)/Figures/"))))[1]
            if size(findall("stc_fctfit_$(SAVENAME)_$(METH).pdf",readdir("$(PATHTOSAVE)/Figures/")[ix]))[1]!=0
                count += 1
            end
        end
        newname = "stc_fctfit_$(SAVENAME)_$(METH)"
        Plots.savefig("$(PATHTOSAVE)/Figures/$(newname).pdf")
    else
        Plots.savefig("$(PATHTOSAVE)/Figures/stc_fctfit_$(SAVENAME)_$(METH).pdf")
    end

    Graphic.StcFctExponent(zeta,zeta[3,1],ORDERS,[0,ORDERS[end]+1],[0,1.8],"Using $(METH)","$SAVENAME")
    if OVERWRITE==true
        Plots.savefig("$(PATHTOSAVE)/Figures/zeta_$(SAVENAME)_$(METH).pdf")
    elseif (OVERWRITE==false && isfile("$(PATHTOSAVE)/Figures/zeta_$(SAVENAME)_$(METH).pdf")==true) #Not sure this is working
        println("THE GIVEN FILE NAME ALREADY EXIST. AN INDICE WILL BE ADDED AT THE END OF THE GIVEN NAME, EQUAL TO THE NUMBER OF FILES WITH THE SAME NAME +1 ")
        count = 1
        for ix=1:size((findall.("zeta_$(SAVENAME)_$(METH).pdf",readdir("$(PATHTOSAVE)/Figures/"))))[1]
            if size(findall("zeta_$(SAVENAME)_$(METH).pdf",readdir("$(PATHTOSAVE)/Figures/")[ix]))[1]!=0
                count += 1
            end
        end
        newname = "zeta_$(SAVENAME)_$(METH)_$(count)"
        Plots.savefig("$(PATHTOSAVE)/Figures/$(newname).pdf")
    else 
        Plots.savefig("$(PATHTOSAVE)/Figures/zeta_$(SAVENAME)_$(METH).pdf")
    end


    Dataprep.write_dat(cat([METHV 0 0 0],cat(zeta,ORDERS,dims=2),dims=1),"$(PATHTOSAVE)/Data/","$(SAVENAME)_stcfct_$(METH)", more=["METHOD $(METH) ; FILE : $(SAVENAME). ROW are results for differents orders which are given at the last column. First column is the exponant, second column is the factor A : Sp(l)=A*l^B. Also, at first row and first column is the method used : <0 for SWO, 0 for raw cube, >0 for PCA. The absolute value gives the parameter (RANGE for SWO, and PC for PCA)" ], overwrite=OVERWRITE)


end #function structure_function





"""
Comparing PCA and SWO method by printing exponant of the structure functions with the order on the same figure
"""
function compmethod_stcfct(VARFILEPATH)
    DATPATH,PCNAME,SWONAME,NFNAME,RAWNAME,PATHTOSAVE,SAVENAME,OVERWRITE = Dataprep.read_var_files(VARFILEPATH)

    # Prepare directories where plots and data will be saved.
    Dataprep.directory_prep(PATHTOSAVE)

    pcdat  = Dataprep.read_dat("$DATPATH/$PCNAME")
    swdat  = Dataprep.read_dat("$DATPATH/$SWONAME")

    NBPC = floor(Int,pcdat[1,1])   # FIRST ROW USED TO KNOW HOW MANY PCs WHERE USED FOR THE RECONSTRUCTION DURING PCA METHOD
    RANGE = floor(Int,swdat[1,1])  # FIRST ROW USED TO KNOW THE SIZE OF THE RANGE USED FOR THE RECONSTRUCTION DURING SWO METHOD
    Graphic.StcFctExponent(pcdat[2:end,:],pcdat[4,1],pcdat[2:end,4],[0,pcdat[:,4][end]+1],[0,2],"Using $(NBPC)PC","$SAVENAME",markers=:rect)
    Graphic.StcFctExponent(swdat[2:end,:],swdat[4,1],swdat[2:end,4],[0,swdat[:,4][end]+1],[0,2],"Using SWO","$SAVENAME",add=true,markers=:diamond)

    if length(NFNAME)!=0
        nfdat = Dataprep.read_dat("$DATPATH/$NFNAME")
        Graphic.StcFctExponent(nfdat[2:end,:],nfdat[4,1],nfdat[2:end,4],[0,nfdat[:,4][end]+1],[0,2],"Using NF","$SAVENAME",add=true)
    end

    if length(RAWNAME)!=0
        raw = Dataprep.read_dat("$DATPATH/$RAWNAME")
        Graphic.StcFctExponent(raw[2:end,:],raw[4,1],raw[2:end,4],[0,raw[:,4][end]+1],[0,2],"Using RAW","$SAVENAME",add=true)
    end

    if OVERWRITE==true

        Plots.savefig("$(PATHTOSAVE)/Figures/zetacomp_$(SAVENAME)_$(NBPC)PC_SWO.pdf")
    elseif (OVERWRITE==false && isfile("$(PATHTOSAVE)/Figures/zetacomp_$(SAVENAME)_$(NBPC)PC_SWO.pdf")==true)
        println("THE GIVEN FILE NAME ALREADY EXIST. AN INDICE WILL BE ADDED AT THE END OF THE GIVEN NAME, EQUAL TO THE NUMBER OF FILES WITH THE SAME NAME +1 ")
        count = 1
        for ix=1:size((findall.("zetacomp_$(SAVENAME)_$(NBPC)PC_SWO.pdf",readdir("$(PATHTOSAVE)/Figures/"))))[1]
            if size(findall("zetacomp_$(SAVENAME)_$(NBPC)PC_SWO.pdf",readdir("$(PATHTOSAVE)/Figures/")[ix]))[1]!=0
                count += 1
            end
        end
        newname = "zetacomp_$(SAVENAME)_$(NBPC)PC_SWO_$(count)"
        Plots.savefig("$(PATHTOSAVE)/Figures/$(newname).pdf")
    else 
        Plots.savefig("$(PATHTOSAVE)/Figures/zetacomp_$(SAVENAME)_$(NBPC)PC_SWO.pdf")

    end
        
end




function prodallvarfile(;PATH=".")
    Dataprep.prodvarfile(PATH="$PATH")
end



"""
WORK IN PROGRESS
Combine multiple CV maps into one cube. 
Use this script in a julia terminal with :
    julia>Unveil.combinecv()
"""
function combinecv(;VARFILEPATH = "varfiles/multipca.txt")
    cvmap,header,dimens = Dataprep.read_fits_pp("$FITSPATH/Data/cv_30PC.fits")
    FITSPATH = "/home/delcamps/Data/Pipe_nebula/TestPCA/"
    cubeout = Array{Float64}(undef,size(cvmap)[1],size(cvmap)[2],4)
    count = 0
    for ix in (30,40,50,60)
        global count +=1
        cvmap,header,dimens = Dataprep.read_fits_pp("$FITSPATH/Data/cv_$(ix)PC.fits")
        cubeout[:,:,count]=cvmap
    end

    Dataprep.write_fits("$(FITSPATH)/Data/cv_30PC.fits","cv_comb","$(PATHTOSAVE)/Data/",cubeout,(size(cubeout)[1],size(cubeout)[2],size(cubeout)[3]),BLANK,finished=true,more=["NBPC",("30,40,50,60")])

end#function combinecv





end # module Unveil
