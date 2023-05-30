include("../src/Data_preparation.jl") #Read and write fits
include("../src/Functionforpca.jl") #Calculations of PCA
#include("../src/Functionforcvi.jl") #Calculations of CVI
#include("../src/Graphic.jl")
#include("../src/Data_analysis.jl")
include("fct_rnd_svd.jl")
using FITSIO, MultivariateStats, Plots, Statistics, Format, Profile, BenchmarkTools, Distributions, GaussianMixtures, StatsBase, ShiftedArrays
using Measures, StatsPlots, KernelDensity, InvertedIndices
using Mmap
using RandomizedLinAlg
#using .Functionforcvi
using .Functionforpca
using .Data_preparation
#using .Graphic
#using .Data_analysis
using .fct_rnd_svd


gr()



pathplot = "/home/delcamps/Prog/CVI_alone/Plots/"
pathdat  = "/home/delcamps/Prog/CVI_alone/Data/"

fits_path =   "/home/delcamps/Data/Blagrave/DHIGLS_DF_Tb.fits" #"/home/delcamps/Data/Blagrave/DHIGLS_DF_Tb.fits" #"/home/delcamps/Prog/Centroid-Velocities-Increment/data/exp_raw/polaris2008.fits" #readline()
#"/home/delcamps/Prog/Centroid-Velocities-Increment/data/exp_raw/polaris2008.fits" #readline()
#/home/delcamps/Data/Blagrave/DHIGLS_DF_Tb.fits
DataNameTitle = "Blagrave_TEST"#"Polaris2008_little_TEST"#"Polaris2008_little_TEST"
vel_units = "km/s"

direc = "TEST_PCA_CONV1"

# Read the fits from the path. Return the data, the VelocityVector, the dimension, the velocity_increment, and the header.
Rawdat,VelocityVector,DataDimension,VelocityIncrement,head = Data_preparation.read_fits_ppv(fits_path,vel_units ; check=false)
Rawdat = reshape(Rawdat,DataDimension[1]*DataDimension[2],DataDimension[3])
#rawdatresh,missing1D,missing2D = Data_preparation.pca_prep(rawdatresh,DataDimension)
#rawdatresh = convert(Array{Float64},rawdatresh)
#DataDimension_nomissing = Data_preparation.read_dim(rawdatresh)
#=
PC = 50
L = PC+2
I = trunc(Int,(DataDimension[3]-PC)/L-1)
ROW = DataDimension_nomissing[1]

# FIRST step : Compute the G matrix
G = fct_rnd_svd.formG(DataDimension[3],L)

# SECOND step : Compute H using A and G
H = fct_rnd_svd.formH(G,rawdatresh,I,L,ROW)
H = convert(AbstractMatrix{Float64},H)

# 3rd step : Form Q using H 
Q = fct_rnd_svd.formQ(H)

# 4th step : Form T using A and Q
T = fct_rnd_svd.formT(rawdatresh,Q)

# 5th step : Form an SVD of T producing Vt, Epst, W
svdobject = RandomizedLinAlg.rsvd_fnkz(T,(I+1)*L)

# 6th step : Compute Ut using Q and W
Ut = fct_rnd_svd.formUt(svdobject,Q)

# 7th step : Retrieve U, V and Eps using Ut, Vt and Epst

# 8th step :  Reconstruct the data with U, V and Eps

=#

# USING RandomizedLinAlg

#Rawdat,missing1D,missing2D = Data_preparation.pca_prep(Rawdat,DataDimension)
Rawdat = convert(Array{Float64},Rawdat)
println("GO randomized svd")
svdobject = RandomizedLinAlg.rsvd_fnkz(Rawdat,50)
Rawdat = 0.0
noisecanals = (1:20,180:250)

l = @layout [a{0.00001h} ; grid(2,2)]
p = plot(layout=l, size=(1000,700),dpi=1000) #,size=(size(rawdat)[1]*25,size(rawdat)[2]*25)
p = plot!(p[1],title="Moments of mean intensities in range [$(noisecanals)] e.g. noise velocity canals \n $(DataNameTitle)",showaxis=:hide,foreground_color_text=:white,seriestype=:scatter, grid=:false,link=:x,titlefontsize=15)
println("GO convergence 1")
Functionforpca.SVDRAND_intensity_moments_specific_canals_withPCA(svdobject,50,(1:30,180:250),DataDimension,p)

#=
test = transpose(svdobject.S.*transpose(svdobject.U))*svdobject.Vt
println("adding blank values")
#test  = Data_preparation.addblank(test,missing2D,-10,DataDimension)
test = reshape(test,DataDimension)
plot(test[50,50,:])
#plot!(reshape(Rawdat,DataDimension)[50,50,:])

h = heatmap(test[:,:,50])
display(h)
sleep(0.5)

#h = heatmap(reshape(Rawdat,DataDimension)[:,:,50])
display(h)
=#