###################################################################
# Plot the heatmap of the convergence matrix produced during the PCA process
# Use this script in a julia terminal with :
#     julia>include("plotcov.jl")
###################################################################
include("../src/Data_preparation.jl") # Read and write fits
include("../src/Functionforpca.jl")   # Calculations of PCA
include("../src/Graphic.jl")   # Calculations of PCA

using FITSIO                                                  # Read Fits
using MultivariateStats, Statistics, StatsBase, Distributions # Statistic
using Profile, BenchmarkTools                                 # Benchmark
using Mmap, DelimitedFiles                                    # Read and write files (.bin and .txt)
using Plots

using .Graphic
using .Functionforpca
using .Data_preparation


println("")
println(" Path to the variable file ? (txt file containing all the informations relevant to read and work on the data)") 
VARFILEPATH = "../varfiles/pca.txt"
GC.gc()
FITSPATH,PATHTOSAVE,UNITVELOCITY,NBPC,BLANK = read_var_files(VARFILEPATH)
(NBPC == 0) && (NBPC="raw")


# Read the fits from the path. Return the data, the VelocityVector, the dimension, the velocity_increment, and the header.
cube,VELOCITYVECTOR,DATADIMENSION,VELOCITYINCREMENT,HEAD = Data_preparation.read_fits_ppv(FITSPATH,UNITVELOCITY ; check=false)


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
plotly()
using Measures
p = plot(plot(mean(cube,dims=1)[:],showlegend=false,ylabel="Mean emission (a.u.)",color=:black),
    heatmap(cov(cube),xlabel="Velocity channel indice",ylabel="Velocity channel indice",grid=true,foreground_color_grid=:white),
      layout=(2,1), link=:x ,size=(1000,500))
display(p)
      # Perform the first PCA analysis (same notation as in the MultivariateStats doc)
#println("Perform PCA")
#M, Yt, VARPERCENT,cubereconstructed = Functionforpca.pca(cube,NBPC)
