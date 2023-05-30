###################################################################
# Plot a map with spectra averaged by region. 
# ADAPTED FOR A SPECIFIC DATASET (Blagrave, HI). Can be adapted easily by modifying the path for save.
# Use this script in a julia terminal with :
#     julia>include("pl_spectra.jl")
###################################################################


include("../src/Data_preparation.jl") #Read and write fits
include("../src/Functionforpca.jl") #Calculations of PCA
include("../src/Functionforcvi.jl") #Calculations of CVI
include("../src/Graphic.jl")
include("../src/Structure_fct.jl")
include("../src/Data_analysis.jl")
include("../src/Functionforfbm.jl")
using FITSIO, MultivariateStats, Plots, Statistics, Format, Profile, BenchmarkTools, Distributions, GaussianMixtures, StatsBase, ShiftedArrays
using Measures, StatsPlots, KernelDensity, InvertedIndices

using .Functionforcvi
using .Functionforpca
using .Data_preparation
using .Graphic
using .Structure_fct
using .Data_analysis
using .Functionforfbm

gr()
blank = -1000

#println(" Path to the fits file ? ") #/home/delcamps/Prog/Centroid-Velocities-Increment/data/alldata/DHIGLS_DF_Tb_extract.fits
pc = 20
fits_path = "/home/delcamps/Prog/CVI_alone/Data/Blagrave_HI_DF/$(pc)_pc.fits"
fits_pathRaw = "/home/delcamps/Prog/Centroid-Velocities-Increment/data/alldata/DHIGLS_DF_Tb_extract.fits"

vel_units = "km/s"

#println("Directory name for plots and fits ?")
#direc = readline()
direc = "Blagrave_HI_DF"

# Prepare directories where plots and data will be saved.
(isdir("/home/delcamps/Prog/CVI_alone/Plots/$(direc)"))==0  && mkdir("/home/delcamps/Prog/CVI_alone/Plots/$(direc)")
(isdir("/home/delcamps/Prog/CVI_alone/Data/$(direc)"))==0   && mkdir("/home/delcamps/Prog/CVI_alone/Data/$(direc)")

Rawdat,VelocityVector,DataDimension,VelocityIncrement,head = Data_preparation.read_fits_ppv(fits_pathRaw,vel_units ; check=false)
pcdat,VelocityVector,DataDimension,VelocityIncrement,head = Data_preparation.read_fits_ppv(fits_path,vel_units ; check=false)

Graphic.pixels_averaged_spectrum(Rawdat,pcdat,5,5,"$(pc)",(50,50),blank,VelocityVector)
savefig("/home/delcamps/Prog/CVI_alone/Plots/$(direc)/$(direc)_meanspectra_$(pc)2.pdf")
