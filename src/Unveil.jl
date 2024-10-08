
module Unveil


include("Dataprep.jl")              # Read and write fits
include("Analysis.jl")              # PCA metric, power_spectra, rms...
include("PCA.jl")                   # Functions for method PCA
include("SWO.jl")                   # Functions for method SWO
include("Graphic.jl")               # Plotting 
include("CVI.jl")                   # Functions for CVI computations
include("Structure_functions.jl")   # Functions for Structure functions computations

using .Dataprep
using .SWO
using .Graphic
using .Analysis
using .Structure_functions
using .CVI
using .PCA
using Plots
using ProgressBars
using StatsBase
using Makie
using DelimitedFiles

export checkPCASWO
export pca
export swo
export convpca
export cvcvi
export cvi 
export multipca
export combinecv
export fitspl
export prodvarfile
export structure_functions
export compmethod_stcfct




"""
    checkPCASWO(VARFILEPATH)

Check PCA and SWO results by plotting randomly taken spectra on the cubes. Plot the noisy original cube superposed by the PCA spectra reconstructed and with the SWO windows. Plot also the residual (original noisy cube less PCA) for the same cubes.

OUTPUTS : 2 plots

Use this script in a julia terminal with :
    julia>Unveil.checkPCASWO("VARFILEPATH")
"""
function checkPCASWO(VARFILEPATH)
    FITSPATH, PCANAME, SWONAME, ORINAME, PATHTOSAVE, SAVENAME, NCOL, NROW, SPECFILE, UNITVELOCITY, NOISECANTXT,OVERWRITE = Dataprep.read_var_files(VARFILEPATH)
    NOISECAN = [parse(Int, ss) for ss in split(NOISECANTXT,",")]


    # Read the fits from the path. Return the data, the VelocityVector, the dimension, the velocity_increment, and the header.
    cubepca,VELOCITYVECTOR,DATADIMENSION,VELOCITYINCREMENT,HEAD = Dataprep.read_fits_ppv("$(FITSPATH)/$(PCANAME)",UNITVELOCITY ; check=false)

    cubeswo,VELOCITYVECTOR,DATADIMENSION,VELOCITYINCREMENT,HEAD = Dataprep.read_fits_ppv("$(FITSPATH)/$(SWONAME)",UNITVELOCITY ; check=false)

    cubeori,VELOCITYVECTOR,DATADIMENSION,VELOCITYINCREMENT,HEAD = Dataprep.read_fits_ppv("$(FITSPATH)/$(ORINAME)",UNITVELOCITY ; check=false)

    # Prepare directories where plots and data will be saved.
    Dataprep.directory_prep(PATHTOSAVE)

    # Replace any NaN value into a missing value and deleted them (can't do PCA on missing values with that package).
    cubepca = Dataprep.replace_nantomissing(cubepca)
    cubeswo = Dataprep.replace_nantomissing(cubeswo)
    cubeori = Dataprep.replace_nantomissing(cubeori)

    cubepca = Dataprep.replace_missingtoblank(cubepca,-1000)
    cubeswo = Dataprep.replace_missingtoblank(cubeswo,-1000)
    cubeori = Dataprep.replace_missingtoblank(cubeori,-1000)
    f = Graphic.spectrePCASWO(cubepca,cubeswo,cubeori,VELOCITYVECTOR,NCOL,NROW,file="$SPECFILE")
    newname = "$(SAVENAME)_CHECKPCASWO"
    if (OVERWRITE==false && isfile("$(PATHTOSAVE)/Figures/$(SAVENAME)_CHECKPCASWO.pdf")==true)
        println("THE GIVEN FILE NAME ALREADY EXIST. AN INDICE WILL BE ADDED AT THE END OF THE GIVEN NAME, EQUAL TO THE NUMBER OF FILES WITH THE SAME NAME +1 ")
        count = 0
        for ix=1:size((findall.("$(SAVENAME)_CHECKPCASWO",readdir("$(PATHTOSAVE)/Figures/"))))[1]
            if size(findall("$(SAVENAME)_CHECKPCASWO",readdir("$(PATHTOSAVE)/Figures/")[ix]))[1]!=0
                count += 1
            end
        end
        newname = "$(newname)_$(count)"
    end #if
    Makie.save("$PATHTOSAVE/Figures/$(SAVENAME)_CHECKPCASWO.pdf",f)



    f = Graphic.residusPCA(cubepca,cubeori,VELOCITYVECTOR,NCOL,NROW,NOISECAN[end],file="$SPECFILE")
    newname = "$(SAVENAME)_RESIDUSPCA"
    if (OVERWRITE==false && isfile("$(PATHTOSAVE)/Figures/$(SAVENAME)_RESIDUSPCA.pdf")==true)
        println("THE GIVEN FILE NAME ALREADY EXIST. AN INDICE WILL BE ADDED AT THE END OF THE GIVEN NAME, EQUAL TO THE NUMBER OF FILES WITH THE SAME NAME +1 ")
        count = 0
        for ix=1:size((findall.("$(SAVENAME)_RESIDUSPCA",readdir("$(PATHTOSAVE)/Figures/"))))[1]
            if size(findall("$(SAVENAME)_RESIDUSPCA",readdir("$(PATHTOSAVE)/Figures/")[ix]))[1]!=0
                count += 1
            end
        end
        newname = "$(newname)_$(count)"
    end #if
    Makie.save("$PATHTOSAVE/Figures/$(SAVENAME)_RESIDUSPCA.pdf",f)

end #checkPCASWO





"""
    compmethod_stcfct(VARFILEPATH)

Comparing PCA and SWO methods by plotting their exponents of structures functions on the same figure. Can also add the exponents from a noise free cube and from a raw cube (e.g. without any method used). A '.txt' file should be used accordingly as an input : use the function 'Unveil.prodallvarfile' to produce it.

OUTPUTS : A plots with exponents of the structures functions of multiple methods.

Use this script in a julia terminal with :
    julia>Unveil.compmethod_stcfct("VARFILEPATH")
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

# Coucou, c'est moi


"""
    convpca(VARFILEPATH)

Compute the PCA convergence criteria based on the matrix projection from PCA. No PCA reconstructed cube will be saved. A '.txt' file should be used accordingly as an input : use the function 'Unveil.prodallvarfile' to produce it. 

OUTPUTS : A plots with every moments of the projection matrix + a plot of the convergence criteria (see the doc). Also add a .dat file with moments, metric and number of PC used for each.

Use this script in a julia terminal with :
    julia>Unveil.convpca("VARFILEPATH")

"""
function convpca(VARFILEPATH; plot=true, onlyplot="")
    FITSPATH,FILENAME,PATHTOSAVE,SAVENAME,NOISECANTXT,UNITVELOCITY,HIGHESTPC,BLANK,OVERWRITE = Dataprep.read_var_files(VARFILEPATH)
    NOISECAN = [parse(Int, ss) for ss in split(NOISECANTXT,",")]

    if onlyplot==""

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

        SIGMAT = Analysis.rms_cube(cube,NOISECAN)[2]

        #First PCA
        println("Perform PCA")
        M, Yt, VARPERCENT,cubereconstructed = PCA.pca(cube,HIGHESTPC)

        # Projection matrix
        proj = PCA.proj(M)
        mom1,mom2,mom3,mom4 = Analysis.fourmoments(proj,dim=2)


        if ismis == 1
            proj = Dataprep.addblank(proj,missingplaces2D[:,1:HIGHESTPC],BLANK,(DATADIMENSION[1],DATADIMENSION[2],HIGHESTPC))
        end
        proj = reshape(proj,(DATADIMENSION[1],DATADIMENSION[2],HIGHESTPC))
        Dataprep.write_fits("$(FITSPATH)/$FILENAME","$(SAVENAME)_projectionmatrix","$(PATHTOSAVE)/Data/",proj,(DATADIMENSION_NOMISSING[1],DATADIMENSION_NOMISSING[2],HIGHESTPC),BLANK,finished=true,overwrite=OVERWRITE)
        xvector = range(1,HIGHESTPC)#[1:HIGHESTPC]




        cubereconstructed = 0
        cube = 0
        missingplaces1D = 0.0 
        missingplaces2D = 0
        M = 0
        Yt = 0
        GC.gc()   # CLEANING MEMORY, NECESSARY FOR LARGE DATASET

        #BL#Graphic.distribmom_multipc(mom1[2:end-1],mom2[2:end-1],mom3[2:end-1],mom4[2:end-1],xvector[2:end-1])
        newname = "$(SAVENAME)_mom"
        if (OVERWRITE==false && isfile("$(PATHTOSAVE)/Figures/$(SAVENAME)_mom.pdf")==true)
            println("THE GIVEN FILE NAME ALREADY EXIST. AN INDICE WILL BE ADDED AT THE END OF THE GIVEN NAME, EQUAL TO THE NUMBER OF FILES WITH THE SAME NAME +1 ")
            count = 0
            for ix=1:size((findall.("$(SAVENAME)_mom",readdir("$(PATHTOSAVE)/Figures/"))))[1]
                if size(findall("$(SAVENAME)_mom",readdir("$(PATHTOSAVE)/Figures/")[ix]))[1]!=0
                    count += 1
                end
            end
            newname = "$(newname)_$(count)"
        end #if
        #BL#Plots.savefig("$(PATHTOSAVE)/Figures/$(newname).pdf")


        #if PCOPT==0
        #    println(" ")
        #    println("From which PC number moments are converging ?")
        #    println("First indice : ")
        #    PCOPT   = parse(Int64,readline())
        #end

        metric = Analysis.metricPCA(mom1,mom2,mom3,mom4,abs(VELOCITYINCREMENT))#,SIGMAT)#abs(VELOCITYINCREMENT))
        #metric = Analysis.metricPCA(mom1,mom2,mom3,mom4)#,SIGMAT)#abs(VELOCITYINCREMENT))
        #metric = Analysis.metricPCA(mom1,mom2,mom3,mom4,abs(VELOCITYINCREMENT))
        #println(metric)
        println("Metric calculated")
        #BL#Graphic.distribcv_multipc(mom1[2:end-1],mom2[2:end-1],mom3[2:end-1],mom4[2:end-1],metric[2:end-1],xvector[2:end-1])
        newname = "$(SAVENAME)_metric"
        if (OVERWRITE==false && isfile("$(PATHTOSAVE)/Figures/$(SAVENAME)_metric.pdf")==true)
            println("THE GIVEN FILE NAME ALREADY EXIST. AN INDICE WILL BE ADDED AT THE END OF THE GIVEN NAME, EQUAL TO THE NUMBER OF FILES WITH THE SAME NAME +1 ")
            count = 0
            for ix=1:size((findall.("$(SAVENAME)_metric",readdir("$(PATHTOSAVE)/Figures/"))))[1]
                if size(findall("$(SAVENAME)_metric",readdir("$(PATHTOSAVE)/Figures/")[ix]))[1]!=0
                    count += 1
                end
            end
            newname = "$(newname)_$(count)"

        end #if

        Dataprep.write_dat([metric mom1 mom2 mom3 mom4.-3 xvector],"$PATHTOSAVE/Data/","$(SAVENAME)_metricPCA",overwrite=OVERWRITE,more=["$FILENAME","Metric  Mom1   Mom2   Mom3   Mom4   PCs"])
        if plot==true
            f=Graphic.convpca("$PATHTOSAVE/Data/$(SAVENAME)_metricPCA.dat")
            newname = "$(SAVENAME)_metric"
            if (OVERWRITE==false && isfile("$(PATHTOSAVE)/Figures/$(SAVENAME)_metric.pdf")==true)
                println("THE GIVEN FILE NAME ALREADY EXIST. AN INDICE WILL BE ADDED AT THE END OF THE GIVEN NAME, EQUAL TO THE NUMBER OF FILES WITH THE SAME NAME +1 ")
                count = 0
                for ix=1:size((findall.("$(SAVENAME)_metric",readdir("$(PATHTOSAVE)/Figures/"))))[1]
                    if size(findall("$(SAVENAME)_metric",readdir("$(PATHTOSAVE)/Figures/")[ix]))[1]!=0
                        count += 1
                    end
                end
                newname = "$(newname)_$(count)"
    
            end #if
            Makie.save("$PATHTOSAVE/Figures/$(SAVENAME)_metric.pdf",f)

            println("First PC with convergence ?")
            FIRST = parse(Int64,readline())
            println("Last PC with convergence ?")
            BURN = parse(Int64,readline())
            BURNMET = BURN
            g=Graphic.convpca("$PATHTOSAVE/Data/$(SAVENAME)_metricPCA.dat",FIRST,BURN,BURNMET=BURNMET)
            newname = "$(SAVENAME)_NOPT"
            if (OVERWRITE==false && isfile("$(PATHTOSAVE)/Figures/$(SAVENAME)_NOPT.pdf")==true)
                println("THE GIVEN FILE NAME ALREADY EXIST. AN INDICE WILL BE ADDED AT THE END OF THE GIVEN NAME, EQUAL TO THE NUMBER OF FILES WITH THE SAME NAME +1 ")
                count = 0
                for ix=1:size((findall.("$(SAVENAME)_NOPT",readdir("$(PATHTOSAVE)/Figures/"))))[1]
                    if size(findall("$(SAVENAME)_NOPT",readdir("$(PATHTOSAVE)/Figures/")[ix]))[1]!=0
                        count += 1
                    end
                end
                newname = "$(newname)_$(count)"
    
            end #if
            Makie.save("$PATHTOSAVE/Figures/$(SAVENAME)_NOPT.pdf",g)

            println("How many rows for Eigenimages plots ?")
            NROW = parse(Int64,readline())
            println("How many columns for Eigenimages plots ?")
            NCOL = parse(Int64,readline())
            g=Graphic.distreigenimage(proj,NROW,NCOL)
            newname = "$(SAVENAME)_distreigenimage"
            if (OVERWRITE==false && isfile("$(PATHTOSAVE)/Figures/$(SAVENAME)_distreigenimage.pdf")==true)
                println("THE GIVEN FILE NAME ALREADY EXIST. AN INDICE WILL BE ADDED AT THE END OF THE GIVEN NAME, EQUAL TO THE NUMBER OF FILES WITH THE SAME NAME +1 ")
                count = 0
                for ix=1:size((findall.("$(SAVENAME)_distreigenimage",readdir("$(PATHTOSAVE)/Figures/"))))[1]
                    if size(findall("$(SAVENAME)_distreigenimage",readdir("$(PATHTOSAVE)/Figures/")[ix]))[1]!=0
                        count += 1
                    end
                end
                newname = "$(newname)_$(count)"
            end #if
            Makie.save("$PATHTOSAVE/Figures/$(SAVENAME)_distreigenimage.pdf",g)

            println("Max of the color scale ?")
            ZLIM = parse(Float64,readline())
            g=Graphic.mapeigenimage(proj,NROW,NCOL,ZLIM)
            newname = "$(SAVENAME)_eigenimage"
            if (OVERWRITE==false && isfile("$(PATHTOSAVE)/Figures/$(SAVENAME)_eigenimage.pdf")==true)
                println("THE GIVEN FILE NAME ALREADY EXIST. AN INDICE WILL BE ADDED AT THE END OF THE GIVEN NAME, EQUAL TO THE NUMBER OF FILES WITH THE SAME NAME +1 ")
                count = 0
                for ix=1:size((findall.("$(SAVENAME)_eigenimage",readdir("$(PATHTOSAVE)/Figures/"))))[1]
                    if size(findall("$(SAVENAME)_eigenimage",readdir("$(PATHTOSAVE)/Figures/")[ix]))[1]!=0
                        count += 1
                    end
                end
                newname = "$(newname)_$(count)"
            end #if
            Makie.save("$PATHTOSAVE/Figures/$(SAVENAME)_eigenimage.pdf",g)

            println("Change color scale ? Y/N")
            CONT = readline()
            while CONT=="Y" || CONT=="y"
                println("Max of the color scale ?")
                ZLIM = parse(Float64,readline())
                g=Graphic.mapeigenimage(proj,NROW,NCOL,ZLIM)
                Makie.save("$PATHTOSAVE/Figures/$(SAVENAME)_eigenimage.pdf",g)
                println("Change color scale ?")
                CONT = readline()
            end
        end 

    else
        file = readdlm("$onlyplot",skipstart=2)
        metric = file[:,1]
        println(size(metric))
        println("First PC with convergence ?")
        FIRST = parse(Int64,readline())
        println("Last PC with convergence ?")
        BURN = parse(Int64,readline())
        BURNMET = BURN
        g=Graphic.convpca("$onlyplot",FIRST,BURN,BURNMET=BURNMET)
        if (OVERWRITE==false && isfile("$(PATHTOSAVE)/Figures/$(SAVENAME)_NOPT.pdf")==true)
            println("THE GIVEN FILE NAME ALREADY EXIST. AN INDICE WILL BE ADDED AT THE END OF THE GIVEN NAME, EQUAL TO THE NUMBER OF FILES WITH THE SAME NAME +1 ")
            count = 0
            for ix=1:size((findall.("$(SAVENAME)_NOPT",readdir("$(PATHTOSAVE)/Figures/"))))[1]
                if size(findall("$(SAVENAME)_NOPT",readdir("$(PATHTOSAVE)/Figures/")[ix]))[1]!=0
                    count += 1
                end
            end
            newname = "$(newname)_$(count)"
        end #if
        Makie.save("$PATHTOSAVE/Figures/$(SAVENAME)_NOPT.pdf",g)

        println("How many rows for Eigenimages plots ?")
            NROW = parse(Int64,readline())
            println("How many columns for Eigenimages plots ?")
            NCOL = parse(Int64,readline())
            g=Graphic.distreigenimage(proj,NROW,NCOL)
            newname = "$(SAVENAME)_distreigenimage"
            if (OVERWRITE==false && isfile("$(PATHTOSAVE)/Figures/$(SAVENAME)_distreigenimage.pdf")==true)
                println("THE GIVEN FILE NAME ALREADY EXIST. AN INDICE WILL BE ADDED AT THE END OF THE GIVEN NAME, EQUAL TO THE NUMBER OF FILES WITH THE SAME NAME +1 ")
                count = 0
                for ix=1:size((findall.("$(SAVENAME)_distreigenimage",readdir("$(PATHTOSAVE)/Figures/"))))[1]
                    if size(findall("$(SAVENAME)_distreigenimage",readdir("$(PATHTOSAVE)/Figures/")[ix]))[1]!=0
                        count += 1
                    end
                end
                newname = "$(newname)_$(count)"
            end #if
            Makie.save("$PATHTOSAVE/Figures/$(SAVENAME)_distreigenimage.pdf",g)

            println("Max of the color scale ?")
            ZLIM = parse(Float64,readline())
            g=Graphic.mapeigenimage(proj,NROW,NCOL,ZLIM)
            newname = "$(SAVENAME)_eigenimage"
            if (OVERWRITE==false && isfile("$(PATHTOSAVE)/Figures/$(SAVENAME)_eigenimage.pdf")==true)
                println("THE GIVEN FILE NAME ALREADY EXIST. AN INDICE WILL BE ADDED AT THE END OF THE GIVEN NAME, EQUAL TO THE NUMBER OF FILES WITH THE SAME NAME +1 ")
                count = 0
                for ix=1:size((findall.("$(SAVENAME)_eigenimage",readdir("$(PATHTOSAVE)/Figures/"))))[1]
                    if size(findall("$(SAVENAME)_eigenimage",readdir("$(PATHTOSAVE)/Figures/")[ix]))[1]!=0
                        count += 1
                    end
                end
                newname = "$(newname)_$(count)"
            end #if
            Makie.save("$PATHTOSAVE/Figures/$(SAVENAME)_eigenimage.pdf",g)

            println("Change color scale ?")
            CONT = readline()
            while CONT=="Y" || CONT=="y"
                println("Max of the color scale ?")
                ZLIM = parse(Float64,readline())
                g=Graphic.mapeigenimage(proj,NROW,NCOL,ZLIM)
                Makie.save("$PATHTOSAVE/Figures/$(SAVENAME)_eigenimage.pdf",g)
                println("Change color scale ?")
                CONT = readline()
            end
    end

end #convpca





"""
    cv(VARFILEPATH)

Calculate the Centroide Velocity (CV) from a cube given in input. A '.txt' file should be used accordingly as an input : use the function 'Unveil.prodallvarfile' to produce it.

OUTPUTS : fits file of the CV map

Use this script in a julia terminal with :
    julia>Unveil.cv(VARFILEPATH)
"""
function cv(VARFILEPATH)
    FITSPATH,FITSNAME,PATHTOSAVE,FITSOURCE,SAVENAME,UNITVELOCITY,THRESHOLD,NOISECANTXT,VSHIFT,BLANK,OVERWRITE = Dataprep.read_var_files(VARFILEPATH)
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


    cubesource = Dataprep.read_fits_ppv("$FITSOURCE",UNITVELOCITY ; check=false)[1]


    # If no threshold, change it to the blank value and SIGMAT to 1, because function moment_one_field will blank every values lower than SIGMAT*THRESHOLD
    if THRESHOLD==0
        THRESHOLD = BLANK
        SIGMAT = 1
    else
        cubesource = Dataprep.replace_nantomissing(cubesource)
        cubesource = Dataprep.replace_blanktomissing(cubesource,BLANK)
        cubesource = Dataprep.pca_prep(cubesource,DATADIMENSION)[1]
        cubesource = convert(Array{Float64},cubesource)

        SIGMAT     = Analysis.rms_cube(cubesource,NOISECAN)[2]
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
    VELOCITYVECTOR = Dataprep.shiftspec(VELOCITYVECTOR,VSHIFT)
    cvmap = CVI.moment_one_field(cube,SIGMAT,THRESHOLD,VELOCITYVECTOR,BLANK) # Calculate the first velocity moment order on data reconstructed
    cube  = 0.0 
    GC.gc()

    if ismis == 1
        cvmap = Dataprep.addblank(cvmap,missingplaces2D[:,1],BLANK,(DATADIMENSION[1],DATADIMENSION[2]))
    end
    cvmap = reshape(cvmap,(DATADIMENSION[1],DATADIMENSION[2]))

    cvmap .= cvmap.+VSHIFT
    VELOCITYVECTOR = Dataprep.shiftspec(VELOCITYVECTOR,-VSHIFT)

    Dataprep.write_fits("$(FITSPATH)/$FITSNAME","CV_$(SAVENAME)_$(METH)","$(PATHTOSAVE)/Data/",cvmap,(DATADIMENSION_NOMISSING[1],DATADIMENSION_NOMISSING[2]),BLANK,finished=true,overwrite=OVERWRITE,more=["THRESH",THRESHOLD])
    #cvmap = 0.0

    println("CV map saved in $(PATHTOSAVE)/Data/CV_$(SAVENAME)_$(METH)_NumberOfFilesWithTheSameNameAsPrefixe.fits as a fits.")


end #cv







"""
    cvcvi(VARFILEPATH)

Calculate the Centroid Velocity (CV) and their increments (CVI) from a cube given in input. A '.txt' file should be used accordingly as an input : use the function 'Unveil.prodallvarfile' to produce it.

OUTPUT : CV map, CVI map with azimutal average, CVI cube with all lags and rotations

Use this script in a julia terminal with :
    julia>Unveil.cvcvi(VARFILEPATH)
"""
function cvcvi(VARFILEPATH)    #  ; thresh=0, add="0")  #<- OPTION used to benchmark intensity threshold easiest    
    FITSPATH,FITSNAME,PATHTOSAVE,FITSOURCE,SAVENAME,THRESHOLD,NOISECANTXT,UNITVELOCITY,REMOVE,VSHIFT,BLANK,LAG,DIFFTYPE,OVERWRITE = Dataprep.read_var_files(VARFILEPATH)
    NOISECAN = [parse(Int, ss) for ss in split(NOISECANTXT,",")]
    if length(LAG)!=1
        LAG = [parse(Int, ss) for ss in split(LAG,",")]
        mult = true
    else
        mult = false
    end 

    
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


    cubesource,VEL,DATADIMSOURCE,INCREM,HEAD = Dataprep.read_fits_ppv("$FITSOURCE",UNITVELOCITY ; check=false)



    #THRESHOLD = thresh
    #cubesource = Dataprep.replace_nantomissing(cubesource)
    #cubesource = Dataprep.replace_blanktomissing(cubesource,BLANK)
    #if any(ismissing,cubesource) 
    #    cubesource = Dataprep.pca_prep(cubesource,DATADIMSOURCE)[1]
    #    cubesource = convert(Array{Float64},cubesource)
    #else
    #    cubesource  = reshape(cubesource,DATADIMSOURCE[1]*DATADIMSOURCE[2],DATADIMSOURCE[3])
    #    cubesource  = convert(Array{Float64},cubesource)
    #end
    #SIGMAT     = Analysis.rms_cube(cubesource,NOISECAN)[2]


    
    if THRESHOLD==0
        THRESHOLD = BLANK
        SIGMAT = 1
        
    else
        cubesource = Dataprep.replace_nantomissing(cubesource)
        cubesource = Dataprep.replace_blanktomissing(cubesource,BLANK)
        if any(ismissing,cubesource) 
            cubesource = Dataprep.pca_prep(cubesource,DATADIMSOURCE)[1]
            cubesource = convert(Array{Float64},cubesource)
        else
            cubesource  = reshape(cubesource,DATADIMSOURCE[1]*DATADIMSOURCE[2],DATADIMSOURCE[3])
            cubesource  = convert(Array{Float64},cubesource)
        end
        SIGMAT     = Analysis.rms_cube(cubesource,NOISECAN)[2]
    end
    
    #THRESHOLD = thresh
    #SIGMAT     = Analysis.rms_cube(cubesource,NOISECAN)[2]

    # Replace any NaN value into a missing value and deleted them (can't do PCA on missing values with that package).


    cube = Dataprep.replace_nantomissing(cube)
    cube = Dataprep.replace_blanktomissing(cube,BLANK)
    
    if REMOVE==true
        SIGMAMAP= Analysis.rms_cube(cubesource,NOISECAN)[1]

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
    VELOCITYVECTOR = Dataprep.shiftspec(VELOCITYVECTOR,VSHIFT)
    cvmap = CVI.moment_one_field(cube,SIGMAT,THRESHOLD,VELOCITYVECTOR,BLANK) # Calculate the first velocity moment order on data reconstructed
    cube  = 0.0 
    GC.gc()

    cvmap .= cvmap.+VSHIFT
    VELOCITYVECTOR = Dataprep.shiftspec(VELOCITYVECTOR,-VSHIFT)

    if ismis == 1
        cvmap = Dataprep.addblank(cvmap,missingplaces2D[:,1],BLANK,(DATADIMENSION[1],DATADIMENSION[2]))
    end
    cvmap = reshape(cvmap,(DATADIMENSION[1],DATADIMENSION[2]))


    #Dataprep.write_fits("$(FITSPATH)/$FITSNAME","CV_$(SAVENAME)_$(METH)","$(PATHTOSAVE)/Data/",cvmap,(DATADIMENSION_NOMISSING[1],DATADIMENSION_NOMISSING[2]),BLANK,finished=true,overwrite=OVERWRITE)
    Dataprep.write_fits("$(FITSPATH)/$FITSNAME","CV_$(SAVENAME)_$(METH)","$(PATHTOSAVE)/Data/",cvmap,(DATADIMENSION_NOMISSING[1],DATADIMENSION_NOMISSING[2]),BLANK,finished=true,overwrite=OVERWRITE,more=["THRESH",THRESHOLD])
    #cvmap = 0.0

    GC.gc()
    println("CV map saved in $(PATHTOSAVE)/Data/CV_$(SAVENAME)_$(METH)_NumberOfFilesWithTheSameNameAsPrefixe.fits as a fits.")

    cvmap = Dataprep.replace_nantomissing(cvmap)
    cvmap = Dataprep.replace_blanktomissing(cvmap,BLANK)

    #powerspec,karr = Analysis.power_spectra(cvmap,DATADIMENSION[1])
    #Graphic.energyspec(powerspec,karr,DATADIMENSION[1],PATHTOSAVE,SAVENAME="$(SAVENAME)")
    cube = 0
    missingplaces1D=0
    missingplaces2D=0

    GC.gc()
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
    
    cvmap = 0
    GC.gc()
    if  mult==true
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

        # REDUCING THE SIZE OF THE ALLANGLE CVI FITS BY REMOVING SOME BLANKING VALUES AND CREATING A 2D ARRAY INSTEAD OF A 3D
        sizee = Array{Float64}(undef,size(LAG)[1])
        for lx=1:size(LAG)[1]
            temp     = size(Dataprep.delete_allnotvalue(cviallangle[:,:,lx],BLANK))[1] |> Int64 #Data without missing value
            sizee[lx]=temp
        end 
        maxi = maximum(sizee) |> Int64
        cviallanglereduced = Array{Float64}(undef,maxi,size(LAG)[1])
        cviallanglereduced .= BLANK
        for lx=1:size(LAG)[1]
            tr = sizee[lx] |> Int64
            cviallanglereduced[1:tr,lx] .= Dataprep.delete_allnotvalue(cviallangle[:,:,lx],BLANK)
        end
        println("Start saving")
        Dataprep.write_fits("$(FITSPATH)/$FITSNAME","CVI$(DIFFTYPE)_$(SAVENAME)_$(METH)","$(PATHTOSAVE)/Data/",cvimap_averaged,(DATADIMENSION[1],DATADIMENSION[2],size(LAG)[1]),BLANK,finished=true,overwrite=OVERWRITE,more=["LAG",LAG,"THRESH",THRESHOLD])
        println("CVI map saved in $(PATHTOSAVE)/Data/CVI$(DIFFTYPE)_$(SAVENAME)_$(METH)_NumberOfFilesWithTheSameNameAsPrefixe.fits as a fits.")
       #Dataprep.write_fits("$(FITSPATH)/$FITSNAME","CVI$(DIFFTYPE)_$(SAVENAME)_allangle_$(METH)","$(PATHTOSAVE)/Data",cviallangle,(DATADIMENSION[1]*DATADIMENSION[2],maximum(NANGLE),size(LAG)[1]),BLANK,finished=true,overwrite=OVERWRITE,more=["LAG",LAG],cvi=true)
        Dataprep.write_fits("$(FITSPATH)/$FITSNAME","CVI$(DIFFTYPE)_$(SAVENAME)_allangle_$(METH)","$(PATHTOSAVE)/Data",cviallanglereduced,(maxi,size(LAG)[1]),BLANK,finished=true,overwrite=OVERWRITE,more=["LAG",LAG,"THRESH",THRESHOLD],cvi=true)
        println("CVI map with all angles values saved in the $(PATHTOSAVE)/Data/CVI$(DIFFTYPE)_$(SAVENAME)_allangle_$(METH)_NumberOfFilesWithTheSameNameAsPrefixe.fits as a fits.")
    
    else
        cvimap_averaged = reshape(cvimap_averaged,DATADIMENSION[1],DATADIMENSION[2])
        cvimap_averaged[:,1:LAG] .= missing
        cvimap_averaged[1:LAG,:] .= missing
        cvimap_averaged[(DATADIMENSION[1]-LAG)+1:DATADIMENSION[1],:] .= missing
        cvimap_averaged[:,(DATADIMENSION[2]-LAG)+1:DATADIMENSION[2]] .= missing
        cvimap_averaged = Dataprep.replace_missingtoblank(cvimap_averaged,BLANK)
        cvimap_averaged = Dataprep.blank_equal(cvimap_averaged,0.0,BLANK)
        cvimap_averaged = convert(Array{Float64},cvimap_averaged)

        cviallangle = reshape(cviallangle,DATADIMENSION[1],DATADIMENSION[2],maximum(NANGLE))
        cviallangle[:,1:LAG,:] .= missing
        cviallangle[1:LAG,:,:] .= missing
        cviallangle[(DATADIMENSION[1]-LAG)+1:DATADIMENSION[1],:,:] .= missing
        cviallangle[:,(DATADIMENSION[2]-LAG)+1:DATADIMENSION[2],:] .= missing
        
        cviallangle = reshape(cviallangle,DATADIMENSION[1]*DATADIMENSION[2],maximum(NANGLE))
        cviallangle = Dataprep.replace_nantoblank(cviallangle,BLANK)
        cviallangle = Dataprep.replace_missingtoblank(cviallangle,BLANK)
        cviallangle = Dataprep.blank_equal(cviallangle,0.0,BLANK)
        cviallangle = convert(Array{Float64},cviallangle)
        #cviallangle = Dataprep.delete_allnotvalue(cviallangle,BLANK)

        # REDUCING THE SIZE OF THE ALLANGLE CVI FITS BY REMOVING SOME BLANKING VALUES AND CREATING A 2D ARRAY INSTEAD OF A 3D
        temp     = size(Dataprep.delete_allnotvalue(cviallangle[:,:],BLANK))[1] |> Int64 #Data without missing value
        sizee=temp
 
        maxi = maximum(sizee) |> Int64
        cviallanglereduced = Array{Float64}(undef,maxi)
        cviallanglereduced .= BLANK
        tr = sizee |> Int64
        cviallanglereduced[1:tr] .= Dataprep.delete_allnotvalue(cviallangle[:,:],BLANK)
        println("Start saving")
        Dataprep.write_fits("$(FITSPATH)/$FITSNAME","CVI$(DIFFTYPE)_$(SAVENAME)_$(METH)","$(PATHTOSAVE)/Data/",cvimap_averaged,(DATADIMENSION[1],DATADIMENSION[2]),BLANK,finished=true,overwrite=OVERWRITE,more=["LAG",LAG,"THRESH",THRESHOLD])
        println("CVI map saved in $(PATHTOSAVE)/Data/CVI$(DIFFTYPE)_$(SAVENAME)_$(METH)_NumberOfFilesWithTheSameNameAsPrefixe.fits as a fits.")

        #Dataprep.write_fits("$(FITSPATH)/$FITSNAME","CVI$(DIFFTYPE)_$(SAVENAME)_$(METH)","$(PATHTOSAVE)/Data/",cvimap_averaged,(DATADIMENSION[1],DATADIMENSION[2]),BLANK,finished=true,overwrite=OVERWRITE,more=["LAG",LAG])
        #println("CVI map saved in $(PATHTOSAVE)/Data/CVI$(DIFFTYPE)_$(SAVENAME)_$(METH)_NumberOfFilesWithTheSameNameAsPrefixe.fits as a fits.")
        Dataprep.write_fits("$(FITSPATH)/$FITSNAME","CVI$(DIFFTYPE)_$(SAVENAME)_allangle_$(METH)","$(PATHTOSAVE)/Data",cviallangle,(DATADIMENSION[1]*DATADIMENSION[2],size(cviallangle)[2]),BLANK,finished=true,overwrite=OVERWRITE,more=["LAG",LAG,"THRESH",THRESHOLD],cvi=true)
        println("CVI map with all angles values saved in the $(PATHTOSAVE)/Data/CVI$(DIFFTYPE)_$(SAVENAME)_allangle_$(METH)_NumberOfFilesWithTheSameNameAsPrefixe.fits as a fits.")
    end

    
end #function cvcvi




"""
    cvi(VARFILEPATH)

Calculate the CVI from a CV map. A '.txt' file should be used accordingly as an input : use the function 'Unveil.prodallvarfile' to produce it.

OUTPUTS : Save the CV fits, 

Use this script in a julia terminal with :
    julia>Unveil.cvi(VARFILEPATH)
"""
function cvi(VARFILEPATH)
    FITSPATH,FITSNAME,PATHTOSAVE,SAVENAME,BLANK,LAG,DIFFTYPE,OVERWRITE = Dataprep.read_var_files(VARFILEPATH)
    if length(LAG)!=1
        LAG = [parse(Int, ss) for ss in split(LAG,",")]
        mult = true
    else
        mult = false
    end 

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
        cviallangle,cvimap_averaged,NANGLE = CVI.construct_cvimap(cvmap,LAG,(DATADIMENSION[1],DATADIMENSION[2]),diff="relative")
        DIFFTYPE = "rel"

    elseif DIFFTYPE=="abs" 
        cviallangle,cvimap_averaged,NANGLE = CVI.construct_cvimap(cvmap,LAG,(DATADIMENSION[1],DATADIMENSION[2]),diff="absolute")
        DIFFTYPE = "abs"

    else
        error("Not good argument in DIFFTYPE (should be abs or relative)")
    end
    cvmap = 0
    GC.gc()
    if  mult==true
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

        # REDUCING THE SIZE OF THE ALLANGLE CVI FITS BY REMOVING SOME BLANKING VALUES AND CREATING A 2D ARRAY INSTEAD OF A 3D
        sizee = Array{Float64}(undef,size(LAG)[1])
        for lx=1:size(LAG)[1]
            temp     = size(Dataprep.delete_allnotvalue(cviallangle[:,:,lx],BLANK))[1] |> Int64 #Data without missing value
            sizee[lx]=temp
        end 
        maxi = maximum(sizee) |> Int64
        cviallanglereduced = Array{Float64}(undef,maxi,size(LAG)[1])
        cviallanglereduced .= BLANK
        for lx=1:size(LAG)[1]
            tr = sizee[lx] |> Int64
            cviallanglereduced[1:tr,lx] .= Dataprep.delete_allnotvalue(cviallangle[:,:,lx],BLANK)
        end
        println("Start saving")
        Dataprep.write_fits("$(FITSPATH)/$FITSNAME","CVI$(DIFFTYPE)_$(SAVENAME)_$(METH)","$(PATHTOSAVE)/Data/",cvimap_averaged,(DATADIMENSION[1],DATADIMENSION[2],size(LAG)[1]),BLANK,finished=true,overwrite=OVERWRITE,more=["LAG",LAG])
        println("CVI map saved in $(PATHTOSAVE)/Data/CVI$(DIFFTYPE)_$(SAVENAME)_$(METH)_NumberOfFilesWithTheSameNameAsPrefixe.fits as a fits.")
       #Dataprep.write_fits("$(FITSPATH)/$FITSNAME","CVI$(DIFFTYPE)_$(SAVENAME)_allangle_$(METH)","$(PATHTOSAVE)/Data",cviallangle,(DATADIMENSION[1]*DATADIMENSION[2],maximum(NANGLE),size(LAG)[1]),BLANK,finished=true,overwrite=OVERWRITE,more=["LAG",LAG],cvi=true)
        Dataprep.write_fits("$(FITSPATH)/$FITSNAME","CVI$(DIFFTYPE)_$(SAVENAME)_allangle_$(METH)","$(PATHTOSAVE)/Data",cviallanglereduced,(maxi,size(LAG)[1]),BLANK,finished=true,overwrite=OVERWRITE,more=["LAG",LAG],cvi=true)
        println("CVI map with all angles values saved in the $(PATHTOSAVE)/Data/CVI$(DIFFTYPE)_$(SAVENAME)_allangle_$(METH)_NumberOfFilesWithTheSameNameAsPrefixe.fits as a fits.")
    
    else
        cvimap_averaged = reshape(cvimap_averaged,DATADIMENSION[1],DATADIMENSION[2])
        cvimap_averaged[:,1:LAG] .= missing
        cvimap_averaged[1:LAG,:] .= missing
        cvimap_averaged[(DATADIMENSION[1]-LAG)+1:DATADIMENSION[1],:] .= missing
        cvimap_averaged[:,(DATADIMENSION[2]-LAG)+1:DATADIMENSION[2]] .= missing
        cvimap_averaged = Dataprep.replace_missingtoblank(cvimap_averaged,BLANK)
        cvimap_averaged = Dataprep.blank_equal(cvimap_averaged,0.0,BLANK)
        cvimap_averaged = convert(Array{Float64},cvimap_averaged)

        cviallangle = reshape(cviallangle,DATADIMENSION[1],DATADIMENSION[2],maximum(NANGLE))
        cviallangle[:,1:LAG,:] .= missing
        cviallangle[1:LAG,:,:] .= missing
        cviallangle[(DATADIMENSION[1]-LAG)+1:DATADIMENSION[1],:,:] .= missing
        cviallangle[:,(DATADIMENSION[2]-LAG)+1:DATADIMENSION[2],:] .= missing
        
        cviallangle = reshape(cviallangle,DATADIMENSION[1]*DATADIMENSION[2],maximum(NANGLE))
        cviallangle = Dataprep.replace_nantoblank(cviallangle,BLANK)
        cviallangle = Dataprep.replace_missingtoblank(cviallangle,BLANK)
        cviallangle = Dataprep.blank_equal(cviallangle,0.0,BLANK)
        cviallangle = convert(Array{Float64},cviallangle)
        #cviallangle = Dataprep.delete_allnotvalue(cviallangle,BLANK)

        # REDUCING THE SIZE OF THE ALLANGLE CVI FITS BY REMOVING SOME BLANKING VALUES AND CREATING A 2D ARRAY INSTEAD OF A 3D
        temp     = size(Dataprep.delete_allnotvalue(cviallangle[:,:],BLANK))[1] |> Int64 #Data without missing value
        sizee=temp
 
        maxi = maximum(sizee) |> Int64
        cviallanglereduced = Array{Float16}(undef,maxi)
        cviallanglereduced .= BLANK
        tr = sizee |> Int64
        cviallanglereduced[1:tr] .= Dataprep.delete_allnotvalue(cviallangle[:,:],BLANK)
        println("Start saving")
        Dataprep.write_fits("$(FITSPATH)/$FITSNAME","CVI$(DIFFTYPE)_$(SAVENAME)_$(METH)","$(PATHTOSAVE)/Data/",cvimap_averaged,(DATADIMENSION[1],DATADIMENSION[2]),BLANK,finished=true,overwrite=OVERWRITE,more=["LAG",LAG])
        println("CVI map saved in $(PATHTOSAVE)/Data/CVI$(DIFFTYPE)_$(SAVENAME)_$(METH)_NumberOfFilesWithTheSameNameAsPrefixe.fits as a fits.")

        #Dataprep.write_fits("$(FITSPATH)/$FITSNAME","CVI$(DIFFTYPE)_$(SAVENAME)_$(METH)","$(PATHTOSAVE)/Data/",cvimap_averaged,(DATADIMENSION[1],DATADIMENSION[2]),BLANK,finished=true,overwrite=OVERWRITE,more=["LAG",LAG])
        #println("CVI map saved in $(PATHTOSAVE)/Data/CVI$(DIFFTYPE)_$(SAVENAME)_$(METH)_NumberOfFilesWithTheSameNameAsPrefixe.fits as a fits.")
        Dataprep.write_fits("$(FITSPATH)/$FITSNAME","CVI$(DIFFTYPE)_$(SAVENAME)_allangle_$(METH)","$(PATHTOSAVE)/Data",cviallangle,(DATADIMENSION[1]*DATADIMENSION[2],size(cviallangle)[2]),BLANK,finished=true,overwrite=OVERWRITE,more=["LAG",LAG],cvi=true)
        println("CVI map with all angles values saved in the $(PATHTOSAVE)/Data/CVI$(DIFFTYPE)_$(SAVENAME)_allangle_$(METH)_NumberOfFilesWithTheSameNameAsPrefixe.fits as a fits.")
    end

end #function cvi





"""
    fitspl(VARFILEPATH)

Do the Sp(l) adjustments on multiple lag ranges. Produce a plot of the compensated Sp(l) with the computed exponents. Plot the exponents obtained at chosen lag ranges.


Use this script in a julia terminal with :
    julia>Unveil.fitspl("VARFILEPATH")
"""
function fitspl(VARFILEPATH)
    DATPATH,DATNAME,PATHTOSAVE,SAVENAME,NORD,NCOL,NROW,LAGTOFITT,OVERWRITE = Dataprep.read_var_files(VARFILEPATH)
    file = readdlm("$DATPATH/$DATNAME",comment_char='#',skipstart=1)
    LAG = file[1,:]
    LAGTOFIT = [parse(Int, ss) for ss in split(LAGTOFITT,",")]
    #LAGTOFIT = LAGTOFIT[1]:LAGTOFIT[2]
    f=Graphic.adjustspless(file,NORD,NCOL,NROW,LAGTOFIT)
    newname = "$(SAVENAME)_adjustspless_fitrange$(LAGTOFIT[1])_$(LAGTOFIT[2])px"
    if (OVERWRITE==false && isfile("$(PATHTOSAVE)/Figures/$(SAVENAME)_adjustspless_fitrange$(LAGTOFIT[1])_$(LAGTOFIT[2])px.pdf")==true)
        println("THE GIVEN FILE NAME ALREADY EXIST. AN INDICE WILL BE ADDED AT THE END OF THE GIVEN NAME, EQUAL TO THE NUMBER OF FILES WITH THE SAME NAME +1 ")
        count = 0
        for ix=1:size((findall.("$(SAVENAME)_adjustspless_fitrange$(LAGTOFIT[1])_$(LAGTOFIT[2])px",readdir("$(PATHTOSAVE)/Figures/"))))[1]
            if size(findall("$(SAVENAME)_adjustspless_fitrange$(LAGTOFIT[1])_$(LAGTOFIT[2])px",readdir("$(PATHTOSAVE)/Figures/")[ix]))[1]!=0
                count += 1
            end
        end
        newname = "$(newname)_$(count)"

    end #if
    Makie.save("$PATHTOSAVE/Figures/$(SAVENAME)_adjustspless_fitrange$(LAGTOFIT[1])_$(LAGTOFIT[2])px.pdf",f)

    println("LAG=$(LAG[2:end])")
    println("Now : $(LAG[LAGTOFIT[1]]):$(LAG[LAGTOFIT[end]])px")
    println("More ranges to fit ? Y/N")
    CONT = readline()
    while CONT=="Y" || CONT=="y"
        println("LAG=$(LAG[2:end])")
        println("First lag for the new fit ? (position, not lag value)")
        LAGTOFIT[1] = parse(Int, readline())
        println("Last lag for the new fit ? (position, not lag value)")
        LAGTOFIT[2] = parse(Int, readline())
        f=Graphic.adjustspless(file,NORD,NCOL,NROW,LAGTOFIT)
        newname = "$(SAVENAME)_adjustspless_fitrange$(LAGTOFIT[1])_$(LAGTOFIT[2])px"
        if (OVERWRITE==false && isfile("$(PATHTOSAVE)/Figures/$(SAVENAME)_adjustspless_fitrange$(LAGTOFIT[1])_$(LAGTOFIT[2])px.pdf")==true)
            println("THE GIVEN FILE NAME ALREADY EXIST. AN INDICE WILL BE ADDED AT THE END OF THE GIVEN NAME, EQUAL TO THE NUMBER OF FILES WITH THE SAME NAME +1 ")
            count = 0
            for ix=1:size((findall.("$(SAVENAME)_adjustspless_fitrange$(LAGTOFIT[1])_$(LAGTOFIT[2])px",readdir("$(PATHTOSAVE)/Figures/"))))[1]
                if size(findall("$(SAVENAME)_adjustspless_fitrange$(LAGTOFIT[1])_$(LAGTOFIT[2])px",readdir("$(PATHTOSAVE)/Figures/")[ix]))[1]!=0
                    count += 1
                end
            end
            newname = "$(newname)_$(count)"
        end #if
        Makie.save("$PATHTOSAVE/Figures/$(SAVENAME)_adjustspless_fitrange$(LAGTOFIT[1])_$(LAGTOFIT[2])px.pdf",f)
        println("More ranges to fit ?")
        CONT = readline()
    end

    println("LAG=$(LAG[2:end])")
    println("Now : $(LAG[LAGTOFIT[1]]):$(LAG[LAGTOFIT[end]])")
    println("First lag for the fit ? (position, not lag value)")
    LAGTOFIT[1] = parse(Int, readline())
    println("Last lag for the fit ? (position, not lag value)")
    LAGTOFIT[2] = parse(Int, readline())
    f = Graphic.expo(file,NORD,LAGTOFIT,LAG)
    if (OVERWRITE==false && isfile("$(PATHTOSAVE)/Figures/$(SAVENAME)_expo.pdf")==true)
        println("THE GIVEN FILE NAME ALREADY EXIST. AN INDICE WILL BE ADDED AT THE END OF THE GIVEN NAME, EQUAL TO THE NUMBER OF FILES WITH THE SAME NAME +1 ")
        count = 0
        for ix=1:size((findall.("$(SAVENAME)_expo",readdir("$(PATHTOSAVE)/Figures/"))))[1]
            if size(findall("$(SAVENAME)_expo",readdir("$(PATHTOSAVE)/Figures/")[ix]))[1]!=0
                count += 1
            end
        end
        newname = "$(newname)_$(count)"

    end #if
    Makie.save("$PATHTOSAVE/Figures/$(SAVENAME)_expo.pdf",f)

    println("LAG=$(LAG[2:end])")
    println("Now : $(LAG[LAGTOFIT[1]]):$(LAG[LAGTOFIT[end]])")
    println("Add other fit ranges ? Y/N")
    CONT = readline()
    while CONT=="y" || CONT=="Y"
        println("LAG=$(LAG[2:end])")
        println("Now : $(LAG[LAGTOFIT[1]]):$(LAG[LAGTOFIT[end]])")

        println("First lag for the new fit ? (position, not lag value)")
        LAGTOFIT[1] = parse(Int, readline())
        println("Last lag for the new fit ? (position, not lag value)")
        LAGTOFIT[2] = parse(Int, readline())
        f=Graphic.expoadd(file,NORD,LAGTOFIT,LAG,f)
        println("Add other fit ranges ? Y/N")
        CONT = readline()
        newname = "$(SAVENAME)_expo"
        if (OVERWRITE==false && isfile("$(PATHTOSAVE)/Figures/$(SAVENAME)_expo.pdf")==true)
            println("THE GIVEN FILE NAME ALREADY EXIST. AN INDICE WILL BE ADDED AT THE END OF THE GIVEN NAME, EQUAL TO THE NUMBER OF FILES WITH THE SAME NAME +1 ")
            count = 0
            for ix=1:size((findall.("$(SAVENAME)_expo",readdir("$(PATHTOSAVE)/Figures/"))))[1]
                if size(findall("$(SAVENAME)_expo",readdir("$(PATHTOSAVE)/Figures/")[ix]))[1]!=0
                    count += 1
                end
            end
            newname = "$(newname)_$(count)"
    
        end #if
    end
    
    Makie.save("$PATHTOSAVE/Figures/$(SAVENAME)_expo.pdf",f)
end







"""
    pca(VARFILEPATH)

Use a PCA (Principal Component Analysis) process on a cube with N PCs given as input. The package MultivariateStats is used to compute PCA. A velocity channel with the values from each pixel is an observation (see the README for a detailled explanation). A '.txt' file should be used accordingly as an input : use the function 'Unveil.prodallvarfile' to produce it.

OUTPUTS : save a cube reconstructed by the number of PC asked in the '.txt' file

Use this function in a julia terminal with :
    julia> Unveil.pca(VARFILEPATH)
"""
function pca(VARFILEPATH)
    FITSPATH,FILENAME,PATHTOSAVE,SAVENAME,UNITVELOCITY,NBPC,BLANK,OVERWRITE = Dataprep.read_var_files(VARFILEPATH)
    (NBPC == 0) && (NBPC="raw")


    # Read the fits from the path. Return the data, the VelocityVector, the dimension, the velocity_increment, and the header.
    cube,VELOCITYVECTOR,DATADIMENSION,VELOCITYINCREMENT,HEAD = Dataprep.read_fits_ppv("$(FITSPATH)/$(FILENAME)",UNITVELOCITY ; check=false)

    # Prepare directories where plots and data will be saved.
    Dataprep.directory_prep(PATHTOSAVE)

    # Replace any NaN value into a missing value and deleted them (can't do PCA on missing values with that package)
    #SIGMAMAP= Analysis.rms_cube(cube,NOISECAN)[1]
    cube = Dataprep.replace_nantomissing(cube)
    cube = Dataprep.replace_blanktomissing(cube,BLANK)
    
    #cube = Dataprep.replace_nosignal(cube,DATADIMENSION,VELOCITYVECTOR,BLANK,SIGMAMAP)


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
    #println(size(Yt))
    Dataprep.write_fits("$(FITSPATH)/$(FILENAME)","Yt_$(NBPC)PC","$PATHTOSAVE/Data/",Yt,(NBPC,DATADIMENSION[3]),BLANK,overwrite=OVERWRITE,more=["NBPC",NBPC,"VARPERC",VARPERCENT[NBPC]*100,"METHOD","PCA"])
    println("Matrix of PCs saved in in $(PATHTOSAVE)/Data/Yt_$(NBPC)PC_NumberOfFilesWithTheSameNameAsPrefixe.fits as a fits.")

    if ismis == 1
        mmean = Dataprep.addblank(M.mean,missingplaces2D[:,1],BLANK,(DATADIMENSION[1],DATADIMENSION[2]))
        #projec            = Dataprep.addblank(PCA.proj(M),missingplaces2D[:,1:NBPC],BLANK,DATADIMENSION) 
        mmean = reshape(mmean,(DATADIMENSION[1],DATADIMENSION[2]))
    else
        mmean = reshape(M.mean,(DATADIMENSION[1],DATADIMENSION[2]))
    end
    Dataprep.write_fits("$(FITSPATH)/$(FILENAME)","mmean_$(NBPC)PC","$PATHTOSAVE/Data/",mmean,(DATADIMENSION[1],DATADIMENSION[2]),BLANK,overwrite=OVERWRITE,more=["NBPC",NBPC,"VARPERC",VARPERCENT[NBPC]*100,"METHOD","PCA"])
    s = open("$(PATHTOSAVE)/Data/Yt_$(NBPC)PC.bin", "w+")
    write(s,Yt)
    close(s)
    Ytpath = "$(PATHTOSAVE)/Data/Yt_$(NBPC)PC.bin"


    proj = PCA.proj(M)
    if ismis == 1
        proj = Dataprep.addblank(proj,missingplaces2D[:,1:NBPC],BLANK,(DATADIMENSION[1],DATADIMENSION[2],NBPC))
    end
    proj = reshape(proj,(DATADIMENSION[1],DATADIMENSION[2],NBPC))
    Dataprep.write_fits("$(FITSPATH)/$FILENAME","$(SAVENAME)_projectionmatrix","$(PATHTOSAVE)/Data/",proj,(DATADIMENSION_NOMISSING[1],DATADIMENSION_NOMISSING[2],NBPC),BLANK,finished=true,overwrite=OVERWRITE)
    xvector = range(1,NBPC)#[1:HIGHESTPC]

    # Cleaning memory
    cube = 0.0 
    Yt   = 0.0
    head = 0.0
    GC.gc()
    #println(size(PCA.proj(M)))
    #println(size(missingplaces2D))
    #println(PCA.proj(M)[:,2])
    if ismis == 1
        cubereconstructed = Dataprep.addblank(cubereconstructed,missingplaces2D,BLANK,DATADIMENSION)
        #projec            = Dataprep.addblank(PCA.proj(M),missingplaces2D[:,1:NBPC],BLANK,DATADIMENSION) 
    end
    cubereconstructed = reshape(cubereconstructed,DATADIMENSION)
    #projec            = reshape(PCA.proj(M),DATADIMENSION)
    #cubereconstructed = Dataprep.blank_equal(cubereconstructed,BLANK,0)
    println("Saving Fits")
    Dataprep.write_fits("$(FITSPATH)/$(FILENAME)","RECONSTRUCTED_$(SAVENAME)_$(NBPC)PC","$PATHTOSAVE/Data/",cubereconstructed,DATADIMENSION,BLANK,overwrite=OVERWRITE,more=["NBPC",NBPC,"VARPERC",VARPERCENT[NBPC]*100,"METHOD","PCA"])
    println("Data reconstructed from PCA saved in $(PATHTOSAVE)/Data/RECONSTRUCTED_$(SAVENAME)_$(NBPC)PC_NumberOfFilesWithTheSameNameAsPrefixe.fits as a fits.")

    #Dataprep.write_fits("$(FITSPATH)/$(FILENAME)","PROJMATRIX_$(SAVENAME)_$(NBPC)PC","$PATHTOSAVE/Data/",projec,DATADIMENSION,BLANK,overwrite=OVERWRITE,more=["NBPC",NBPC,"VARPERC",VARPERCENT[NBPC]*100,"METHOD","PCA"])
    #println("Projection matrix from PCA saved in $(PATHTOSAVE)/Data/PROJMATRIX_$(SAVENAME)_$(NBPC)PC_NumberOfFilesWithTheSameNameAsPrefixe.fits as a fits.")
    cubereconstructed = 0
    GC.gc()

    Graphic.pratio(M,true,NBPC,"Variance reproduction with the number of PC used")
    newname = "pratio_$(SAVENAME)_$(NBPC)PC"
    if (OVERWRITE==false && isfile("$(PATHTOSAVE)/Figures/pratio_$(SAVENAME)_$(NBPC)PC.pdf")==true)
        println("THE GIVEN FILE NAME ALREADY EXIST. AN INDICE WILL BE ADDED AT THE END OF THE GIVEN NAME, EQUAL TO THE NUMBER OF FILES WITH THE SAME NAME +1 ")
        count = 1
        for ix=1:size((findall.("pratio_$(SAVENAME)_$(NBPC)PC.pdf",readdir("$(PATHTOSAVE)/Figures/"))))[1]
            if size(findall("pratio_$(SAVENAME)_$(NBPC)PC.pdf",readdir("$(PATHTOSAVE)/Figures/")[ix]))[1]!=0
                count += 1
            end
        end
        newname = "pratio_$(SAVENAME)_$(NBPC)PC_$(count)" 
    end
    Plots.savefig("$(PATHTOSAVE)/Figures/$(newname).pdf")

    println("Figure of pratio saved in $(PATHTOSAVE)/Figures/pratio_$(SAVENAME)_$(NBPC)PC.pdf")
end   #pca





"""
    prodallvarfile(;PATH=".",com=true)

Produces all .txt varfile needed to run Unveil.jl. Can choose the path (local by default, variable name is PATH), and the explanation for each variable which will be written as commentary (default is true, variable name is com).
Use this script in a julia terminal with :
    julia>Unveil.prodallvarfile(PATH="",com=false/true)
"""
function prodallvarfile(;PATH=".",com=true)
    Dataprep.prodvarfile(PATH="$PATH",com=com)
end





"""
    structure_functions(VARFILEPATH; meth="moninyaglom")

Compute the structure functions ``S_p(l)`` of a Centroid Velocity Increment cube. By default will use the definition in Monin & Yaglom+75, but can compute the one in  NOT WORKING ANYMORE NEED DEBUGING Hily-Blant+2008 by changing the option 'meth=hily' NOT WORKING ANYMORE NEED DEBUGING. A '.txt' file should be used accordingly as an input : use the function 'Unveil.prodallvarfile' to produce it. The cube given as input should be a CVI (use function **Unveil.cvi** if needed). Prefer a cube with rotations of every lags than azimutal average, like that : (Pixel positions,angles,lag). 


Use this function in a julia terminal with :
    julia> Unveil.structure_functions(VARFILEPATH)
"""
function structure_functions(VARFILEPATH ; meth="moninyaglom", plot=true,LAGMAX=0,NWIN=10)
    FITSPATH,FILENAME,PATHTOSAVE,SAVENAME,ORDERSTXT,BLANK,NROW,NCOL,OVERWRITE = Dataprep.read_var_files(VARFILEPATH)
    ORDERS = [parse(Int, ss) for ss in split(ORDERSTXT,",")]


    # Read the fits from the path. Return the data, the VelocityVector, the dimension, the velocity_increment, and the header.
    cvicube,DATADIMENSION,HEAD = Dataprep.read_fits_cvi("$(FITSPATH)/$(FILENAME)" ; check=false)
    haskey(HEAD,"THRESH") && (THRESHOLD = HEAD["THRESH"])
    haskey(HEAD,"THRESH") || (THRESHOLD = 0)
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

    if plot==true
        for lx in ProgressBar(1:size(LAG)[1])
            f=Graphic.integrantspl(cvicube[:,lx],NCOL,NROW)
            newname = "$(SAVENAME)_intsplLAG$(LAG[lx])px"
            if (OVERWRITE==false && isfile("$(PATHTOSAVE)/Figures/$(SAVENAME)_intsplLAG$(LAG[lx])px.pdf")==true)
                println("THE GIVEN FILE NAME ALREADY EXIST. AN INDICE WILL BE ADDED AT THE END OF THE GIVEN NAME, EQUAL TO THE NUMBER OF FILES WITH THE SAME NAME +1 ")
                count = 0
                for ix=1:size((findall.("$(SAVENAME)_intsplLAG$(LAG[lx])px",readdir("$(PATHTOSAVE)/Figures/"))))[1]
                    if size(findall("$(SAVENAME)_intsplLAG$(LAG[lx])px",readdir("$(PATHTOSAVE)/Figures/")[ix]))[1]!=0
                        count += 1
                    end
                end
                newname = "$(newname)_$(count)"

            end #if
            Makie.save("$PATHTOSAVE/Figures/$(SAVENAME)_intspl$(LAG[lx])px.pdf",f)
        end


        for lx in ProgressBar(1:size(LAG)[1])
            f=Graphic.convintegrantspl(cvicube[:,lx],NCOL,NROW,NWIN)
            newname = "$(SAVENAME)_convsplLAG$(LAG[lx])px"
            if (OVERWRITE==false && isfile("$(PATHTOSAVE)/Figures/$(SAVENAME)_convsplLAG$(LAG[lx])px.pdf")==true)
                println("THE GIVEN FILE NAME ALREADY EXIST. AN INDICE WILL BE ADDED AT THE END OF THE GIVEN NAME, EQUAL TO THE NUMBER OF FILES WITH THE SAME NAME +1 ")
                count = 0
                for ix=1:size((findall.("$(SAVENAME)_convsplLAG$(LAG[lx])px",readdir("$(PATHTOSAVE)/Figures/"))))[1]
                    if size(findall("$(SAVENAME)_convsplLAG$(LAG[lx])px",readdir("$(PATHTOSAVE)/Figures/")[ix]))[1]!=0
                        count += 1
                    end
                end
                newname = "$(newname)_$(count)"

            end #if
            Makie.save("$PATHTOSAVE/Figures/$(SAVENAME)_convspl$(LAG[lx])px.pdf",f)
        end
    end
    cvicube = Dataprep.replace_nantomissing(cvicube)
    cvicube = Dataprep.replace_blanktomissing(cvicube,BLANK)
    cvicube = Dataprep.replace_blanktomissing(cvicube,0)

    #println(cvicube)
    if meth=="moninyaglom"
        sct = Structure_functions.fct_sct(cvicube,LAG,ORDERS)  
    elseif meth=="hily"
        sct = Structure_functions.fct_sct_int(cvicube,LAG,ORDERS) 
    else
        error("The method given as an option when calling the structure_functions function is not correct. Please, use 'moninyaglom' or 'hily'")
    end
    nl = cat(0,LAG ; dims=1)
    nsct = cat(reshape(nl,(1,size(nl)[1])),cat(ORDERS,sct ; dims=2);dims=1)

    Dataprep.write_dat(nsct,"$(PATHTOSAVE)/Data/","$(SAVENAME)_Sp(l)_$(METH)", more=["METHOD $(METH) ; FILE : $(SAVENAME) ; Intensity threshold during CV computation : $(THRESHOLD). Each column is a lag, each row an order. First column give the orders, first row the lags. For information : p=$(ORDERS), and l=$(LAG) ;  " ], overwrite=OVERWRITE)

    if plot==true
        f=Graphic.spless(nsct,ORDERS[end],NCOL,NROW,LAG=LAGMAX)
        newname = "$(SAVENAME)_splESS"
        if (OVERWRITE==false && isfile("$(PATHTOSAVE)/Figures/$(SAVENAME)_splESS.pdf")==true)
            println("THE GIVEN FILE NAME ALREADY EXIST. AN INDICE WILL BE ADDED AT THE END OF THE GIVEN NAME, EQUAL TO THE NUMBER OF FILES WITH THE SAME NAME +1 ")
            count = 0
            for ix=1:size((findall.("$(SAVENAME)_splESS",readdir("$(PATHTOSAVE)/Figures/"))))[1]
                if size(findall("$(SAVENAME)_splESS",readdir("$(PATHTOSAVE)/Figures/")[ix]))[1]!=0
                    count += 1
                end
            end
            newname = "$(newname)_$(count)"

        end #if
        Makie.save("$PATHTOSAVE/Figures/$(SAVENAME)_splESS.pdf",f)

       
    end             #if plot==true




end #function structure_function










"""
    swo(VARFILEPATH ; meth="swo")

Use a SWO (Spectral Window Optimisation) process on a cube. See the README for a detailled working process. Two methods can be computed : the one presented in the README (meth="swo", by default) or the one presented in Pety+03 (meth="pety" DEPRECATED, NEED TO BE REWORKED) A '.txt' file should be used accordingly as an input : use the function 'Unveil.prodallvarfile' to produce it.

Use this function in a julia terminal with :
    julia> Unveil.swo(VARFILEPATH)

"""
function swo(VARFILEPATH ; meth="swo")   
    FITSPATH,FILENAME,PATHTOSAVE,SAVENAME,UNITVELOCITY,BLANK,NOISECANTXT,OVERWRITE = Dataprep.read_var_files(VARFILEPATH)
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
    
    if meth=="swo"
       maskinterv,mask,posimap = SWO.newswo(cube,DATADIMENSION_NOMISSING,VELOCITYVECTOR,NOISECAN)
    elseif meth=="pety"
        maskinterv,mask = SWO.petysnr(cube,DATADIMENSION_NOMISSING,VELOCITYVECTOR,NOISECAN)
    else
        error("The method given as an option when calling the swo function is not correct. Please, use 'swo' or 'pety'")
    end
    #maskinterv,maskpety = SWO.petysnr(cube,DATADIMENSION_NOMISSING,VELOCITYVECTOR,NOISECAN)


    if ismis == 1
        maskinterv = Dataprep.addblank(maskinterv,missingplaces2D,BLANK,DATADIMENSION)
        posimapinf = Dataprep.addblank(posimap[:,1],missingplaces2D[:,1],BLANK,(DATADIMENSION[1],DATADIMENSION[2]))
        posimapsup = Dataprep.addblank(posimap[:,2],missingplaces2D[:,1],BLANK,(DATADIMENSION[1],DATADIMENSION[2]))
        #maskintervpety = Dataprep.addblank(maskintervpety,missingplaces2D,BLANK,DATADIMENSION)
    else 
        posimapinf = posimap[:,1]
        posimapsup = posimap[:,2]
    end

    maskinterv = reshape(maskinterv,DATADIMENSION)
    maskinterv = Dataprep.blank_equal(maskinterv,BLANK,0)
    posimap = Array{Float64}(undef, (DATADIMENSION[1],DATADIMENSION[2],2))
    posimap .= BLANK
    posimap[:,:,1] .= reshape(posimapinf,(DATADIMENSION[1],DATADIMENSION[2])) 
    posimap[:,:,2] .= reshape(posimapsup,(DATADIMENSION[1],DATADIMENSION[2]))
    posimap = Dataprep.replace_blanktomissing(posimap,BLANK)
    int = 10
    #posimap = Dataprep.replace_blanktomissing(posimap,0)
    #p = plot(layout=2)
    #p = heatmap!(p[1],posimap[:,:,1],clims=(1,100))
    #p = heatmap!(p[2],posimap[:,:,2],clims=(110,170))
    #display(p)


    # These next ~40 rows allow to treat spectra for which SWO didn't found an optimised window. Will average the positions that SWO found for the 5x5 pixels around the pixel without a window. Thus, can't work if multiples spectra doesn't have a window around them, but it is logical. 
    # Would be better to include it in a dedicated function.
    for px=1:size(maskinterv)[1]
        for py=1:size(maskinterv)[2]
            if maskinterv[px,py,3]!=0 && maskinterv[px,py,3]!=BLANK
                if px>int && px<DATADIMENSION[1]-int && py>int && py<DATADIMENSION[2]-int
                    posi = floor(moment(collect(skipmissing(posimap[px-int:px+int,py-int:py+int,1])),1,0)) |> Int64
                    posf = floor(moment(collect(skipmissing(posimap[px-int:px+int,py-int:py+int,2])),1,0))+1 |> Int64
                elseif px<int && py>int && py<DATADIMENSION[2]-int
                    posi = floor(moment(collect(skipmissing(posimap[px:px+int,py-int:py+int,1])),1,0)) |> Int64
                    posf = floor(moment(collect(skipmissing(posimap[px:px+int,py-int:py+int,2])),1,0))+1 |> Int64
                elseif py<int && px>int &&  px<DATADIMENSION[1]-int 
                    posi = floor(moment(collect(skipmissing(posimap[px-int:px+int,py:py+int,1])),1,0)) |> Int64
                    posf = floor(moment(collect(skipmissing(posimap[px-int:px+int,py:py+int,2])),1,0))+1 |> Int64
                elseif px>DATADIMENSION[1]-int && py>int && py<DATADIMENSION[2]-int
                    posi = floor(moment(collect(skipmissing(posimap[px-int:px,py-int:py+int,1])),1,0)) |> Int64
                    posf = floor(moment(collect(skipmissing(posimap[px-int:px,py-int:py+int,2])),1,0))+1 |> Int64
                elseif py>DATADIMENSION[2]-int && px>int &&  px<DATADIMENSION[1]-int
                    posi = floor(moment(collect(skipmissing(posimap[px-int:px+int,py-int:py,1])),1,0)) |> Int64
                    posf = floor(moment(collect(skipmissing(posimap[px-int:px+int,py-int:py,2])),1,0))+1 |> Int64
                elseif px==int && py>int && py<DATADIMENSION[2]-int
                    posi = floor(moment(collect(skipmissing(posimap[px-int+1:px+int,py-int:py+int,1])),1,0)) |> Int64
                    posf = floor(moment(collect(skipmissing(posimap[px-int+1:px+int,py-int:py+int,2])),1,0))+1 |> Int64
                elseif px==DATADIMENSION[1]-int && py>int && py<DATADIMENSION[2]-int
                    posi = floor(moment(collect(skipmissing(posimap[px-int:px+int-1,py-int:py+int,1])),1,0)) |> Int64
                    posf = floor(moment(collect(skipmissing(posimap[px-int:px+int-1,py-int:py+int,2])),1,0))+1 |> Int64
                elseif py==int && px>int && px<DATADIMENSION[1]-int
                    posi = floor(moment(collect(skipmissing(posimap[px-int:px+int,py-int+1:py+int,1])),1,0)) |> Int64
                    posf = floor(moment(collect(skipmissing(posimap[px-int:px+int,py-int+1:py+int,2])),1,0))+1 |> Int64
                elseif py==DATADIMENSION[2]-int && px>int && px<DATADIMENSION[1]-int
                    posi = floor(moment(collect(skipmissing(posimap[px-int:px+int,py-int:py+int-1,1])),1,0)) |> Int64
                    posf = floor(moment(collect(skipmissing(posimap[px-int:px+int,py-int:py+int-1,2])),1,0))+1 |> Int64
                elseif px==int && py==int
                    posi = floor(moment(collect(skipmissing(posimap[px-int+1:px+int,py-int+1:py+int,1])),1,0)) |> Int64
                    posf = floor(moment(collect(skipmissing(posimap[px-int+1:px+int,py-int+1:py+int,2])),1,0))+1 |> Int64
                elseif px==int && py==DATADIMENSION[2]-int
                    posi = floor(moment(collect(skipmissing(posimap[px-int+1:px+int,py-int:py+int-1,1])),1,0)) |> Int64
                    posf = floor(moment(collect(skipmissing(posimap[px-int+1:px+int,py-int:py+int-1,2])),1,0))+1 |> Int64
                elseif px==DATADIMENSION[2]-int && py==int
                    posi = floor(moment(collect(skipmissing(posimap[px-int:px+int-1,py-int+1:py+int,1])),1,0)) |> Int64
                    posf = floor(moment(collect(skipmissing(posimap[px-int:px+int-1,py-int+1:py+int,2])),1,0))+1 |> Int64
                elseif px==DATADIMENSION[2]-int && py==DATADIMENSION[2]-int
                    posi = floor(moment(collect(skipmissing(posimap[px-int:px+int-1,py-int:py+int-1,1])),1,0)) |> Int64
                    posf = floor(moment(collect(skipmissing(posimap[px-int:px+int-1,py-int:py+int-1,2])),1,0))+1 |> Int64
                elseif py>=DATADIMENSION[2]-int && px<=int 
                    posi = floor(moment(collect(skipmissing(posimap[px:px+int,py-int+1:py,1])),1,0)) |> Int64
                    posf = floor(moment(collect(skipmissing(posimap[px:px+int,py-int+1:py,2])),1,0))+1 |> Int64
                elseif py>=DATADIMENSION[2]-int && px>=DATADIMENSION[2]-int
                    posi = floor(moment(collect(skipmissing(posimap[px-int:px,py-int+1:py,1])),1,0)) |> Int64
                    posf = floor(moment(collect(skipmissing(posimap[px-int:px,py-int+1:py,2])),1,0))+1 |> Int64
                elseif py<=int && px<=int 
                    posi = floor(moment(collect(skipmissing(posimap[px:px+int,py:py+int,1])),1,0)) |> Int64
                    posf = floor(moment(collect(skipmissing(posimap[px:px+int,py:py+int,2])),1,0))+1 |> Int64
                elseif py<=int && px>=DATADIMENSION[2]-int
                    posi = floor(moment(collect(skipmissing(posimap[px-int:px,py:py+int,1])),1,0)) |> Int64
                    posf = floor(moment(collect(skipmissing(posimap[px-int:px,py:py+int,2])),1,0))+1 |> Int64
                else 
                    posi = 2
                    posf = DATADIMENSION[3]-1
                end
                
                if posf<=0
                    posf = DATADIMENSION[3]-1
                end
                
                if posi<=0
                    posi = 2
                end
                maskinterv[px,py,1:posi] .= 0
                maskinterv[px,py,posf:end] .=0
            end
        end
    end
    
    #maskintervpety = reshape(maskintervpety,DATADIMENSION)
    #maskintervpety = Dataprep.blank_equal(maskintervpety,BLANK,0)

    Dataprep.write_fits("$(FITSPATH)/$FILENAME","RECONSTRUCTED_$(SAVENAME)_SWO","$PATHTOSAVE/Data/",maskinterv,DATADIMENSION,BLANK,overwrite=OVERWRITE,more=["METHOD","SWO"])
    println("Data reconstructed from SWO method saved in $(PATHTOSAVE)/Data/RECONSTRUCTED_$(SAVENAME)_SWO_NumberOfFilesWithTheSameNameAsPrefixe.fits as a fits.")
    #Dataprep.write_fits("$(FITSPATH)/$FILENAME","RECONSTRUCTED_$(SAVENAME)_SWOpety","$PATHTOSAVE/Data/",maskintervpety,DATADIMENSION,BLANK,overwrite=OVERWRITE,more=["METHOD","SWO"])

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


function cvcviOT(VARFILEPATH ; thresh=0, add="0")
    FITSPATH,FITSNAME,PATHTOSAVE,FITSOURCE,SAVENAME,THRESHOLD,NOISECANTXT,UNITVELOCITY,REMOVE,BLANK,LAG,DIFFTYPE,OVERWRITE = Dataprep.read_var_files(VARFILEPATH)
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


    cubesource,VEL,DATADIMSOURCE,INCREM,HEAD = Dataprep.read_fits_ppv("$FITSOURCE",UNITVELOCITY ; check=false)


    if thresh==0
        THRESHOLD = BLANK
        SIGMAT = 1
        
    else
        THRESHOLD = thresh
        cubesource = Dataprep.replace_nantomissing(cubesource)
        cubesource = Dataprep.replace_blanktomissing(cubesource,BLANK)
        if any(ismissing,cubesource) 
            cubesource = Dataprep.pca_prep(cubesource,DATADIMSOURCE)[1]
            cubesource = convert(Array{Float64},cubesource)
        else
            cubesource  = reshape(cubesource,DATADIMSOURCE[1]*DATADIMSOURCE[2],DATADIMSOURCE[3])
            cubesource  = convert(Array{Float64},cubesource)
        end
        SIGMAT     = Analysis.rms_cube(cubesource,NOISECAN)[2]

    end


    # Replace any NaN value into a missing value and deleted them (can't do PCA on missing values with that package).


    cube = Dataprep.replace_nantomissing(cube)
    cube = Dataprep.replace_blanktomissing(cube,BLANK)
    
    if REMOVE==true
        SIGMAMAP= Analysis.rms_cube(cubesource,NOISECAN)[1]

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


    #Dataprep.write_fits("$(FITSPATH)/$FITSNAME","CV_$(SAVENAME)_$(METH)","$(PATHTOSAVE)/Data/",cvmap,(DATADIMENSION_NOMISSING[1],DATADIMENSION_NOMISSING[2]),BLANK,finished=true,overwrite=OVERWRITE)
    Dataprep.write_fits("$(FITSPATH)/$FITSNAME","CV_$(SAVENAME)_$(METH)_$add","$(PATHTOSAVE)/Data/",cvmap,(DATADIMENSION_NOMISSING[1],DATADIMENSION_NOMISSING[2]),BLANK,finished=true,overwrite=OVERWRITE)
    #cvmap = 0.0

    GC.gc()
    println("CV map saved in $(PATHTOSAVE)/Data/CV_$(SAVENAME)_$(METH)_NumberOfFilesWithTheSameNameAsPrefixe.fits as a fits.")

    cvmap = Dataprep.replace_nantomissing(cvmap)
    cvmap = Dataprep.replace_blanktomissing(cvmap,BLANK)

    powerspec,karr = Analysis.power_spectra(cvmap,DATADIMENSION[1])
    #Graphic.energyspec(powerspec,karr,DATADIMENSION[1],PATHTOSAVE,SAVENAME="$(SAVENAME)")

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


    #Dataprep.write_fits("$(FITSPATH)/$FITSNAME","CVI$(DIFFTYPE)_$(SAVENAME)_$(METH)_$(THRESHOLD*1e2)","$(PATHTOSAVE)/Data/",cvimap_averaged,(DATADIMENSION[1],DATADIMENSION[2],size(LAG)[1]),BLANK,finished=true,overwrite=OVERWRITE,more=["LAG",LAG])
    #println("CVI map saved in $(PATHTOSAVE)/Data/CVI$(DIFFTYPE)_$(SAVENAME)_$(METH)_NumberOfFilesWithTheSameNameAsPrefixe.fits as a fits.")
#
    #Dataprep.write_fits("$(FITSPATH)/$FITSNAME","CVI$(DIFFTYPE)_$(SAVENAME)_allangle_$(METH)_$(THRESHOLD*1e2)","$(PATHTOSAVE)/Data",cviallangle,(DATADIMENSION[1]*DATADIMENSION[2],maximum(NANGLE),size(LAG)[1]),BLANK,finished=true,overwrite=OVERWRITE,more=["LAG",LAG])
    #println("CVI map with all angles values saved in the $(PATHTOSAVE)/Data/CVI$(DIFFTYPE)_$(SAVENAME)_allangle_$(METH)_NumberOfFilesWithTheSameNameAsPrefixe.fits as a fits.")
    
    Dataprep.write_fits("$(FITSPATH)/$FITSNAME","CVI$(DIFFTYPE)_$(SAVENAME)_$(METH)_$add","$(PATHTOSAVE)/Data/",cvimap_averaged,(DATADIMENSION[1],DATADIMENSION[2],size(LAG)[1]),BLANK,finished=true,overwrite=OVERWRITE,more=["LAG",LAG])
    println("CVI map saved in $(PATHTOSAVE)/Data/CVI$(DIFFTYPE)_$(SAVENAME)_$(METH)_NumberOfFilesWithTheSameNameAsPrefixe.fits as a fits.")

    Dataprep.write_fits("$(FITSPATH)/$FITSNAME","CVI$(DIFFTYPE)_$(SAVENAME)_allangle_$(METH)_$add","$(PATHTOSAVE)/Data",cviallangle,(DATADIMENSION[1]*DATADIMENSION[2],maximum(NANGLE),size(LAG)[1]),BLANK,finished=true,overwrite=OVERWRITE,more=["LAG",LAG])
    println("CVI map with all angles values saved in the $(PATHTOSAVE)/Data/CVI$(DIFFTYPE)_$(SAVENAME)_allangle_$(METH)_NumberOfFilesWithTheSameNameAsPrefixe.fits as a fits.")
    #
    #Dataprep.write_fits("$(FITSPATH)/$FITSNAME","CVI$(DIFFTYPE)_$(SAVENAME)_$(METH)","$(PATHTOSAVE)/Data/",cvimap_averaged,(DATADIMENSION[1],DATADIMENSION[2],size(LAG)[1]),BLANK,finished=true,overwrite=OVERWRITE,more=["LAG",LAG])
    #println("CVI map saved in $(PATHTOSAVE)/Data/CVI$(DIFFTYPE)_$(SAVENAME)_$(METH)_NumberOfFilesWithTheSameNameAsPrefixe.fits as a fits.")
#
    #Dataprep.write_fits("$(FITSPATH)/$FITSNAME","CVI$(DIFFTYPE)_$(SAVENAME)_allangle_$(METH)","$(PATHTOSAVE)/Data",cviallangle,(DATADIMENSION[1]*DATADIMENSION[2],maximum(NANGLE),size(LAG)[1]),BLANK,finished=true,overwrite=OVERWRITE,more=["LAG",LAG])
    #println("CVI map with all angles values saved in the $(PATHTOSAVE)/Data/CVI$(DIFFTYPE)_$(SAVENAME)_allangle_$(METH)_NumberOfFilesWithTheSameNameAsPrefixe.fits as a fits.")
end #function cvcvi






function newcvi(fits,LAG,DLAG,BLANK)
    CVMAP,HEAD,DATADIMENSION = Dataprep.read_fits_pp(fits)
    #cviall = CVI.newcvicalc(cvmap,STEP,-1000)
    cviall = CVI.nncvi(CVMAP,LAG,DLAG,BLANK)
    Dataprep.write_fits(fits,"NEWCVI","/home/delcamps/Prog/test/",cviall,(size(cviall)[1],size(cviall)[2]),BLANK,finished=true,overwrite=false,more=["LAG",LAG])

end

end # module Unveil
