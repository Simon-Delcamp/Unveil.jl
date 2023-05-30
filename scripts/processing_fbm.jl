###################################################################
# Script for first looking at CV, CVI and structure function of fBms cubes. 
# WORK IN PROGRESS
# Use this script in a julia terminal with :
#     julia>include("processing_fbm.jl")
###################################################################

include("../src/Data_preparation.jl") #Read and write fits
include("../src/Functionforpca.jl") #Calculations of PCA
include("../src/Functionforcvi.jl") #Calculations of CVI
include("../src/Graphic.jl")
include("../src/Structure_fct.jl")
include("../src/Data_analysis.jl")
include("../src/Functionforfbm.jl")
using FITSIO, MultivariateStats, Plots, Statistics, Format, Profile, BenchmarkTools, Distributions, GaussianMixtures, StatsBase, ShiftedArrays
using Measures, StatsPlots, LsqFit, LaTeXStrings
using FFTW, BenchmarkTools

using .Functionforcvi
using .Functionforpca
using .Data_preparation
using .Graphic
using .Structure_fct
using .Data_analysis
using .Functionforfbm

gr()

imsize   = 500
Np       = trunc(Int,imsize/2)
powerlaw = -3
blank    = -10000
Lag      = [3,8,16,24,40,60,100,200]
StcFctOrder = [1,2,3]
DirectoryName = "fbm_sim_imsize$(imsize)_powerlaw$powerlaw"
DataNameTitle  = "$(DirectoryName) | $(imsize)imagesize "

# Data reading
Data_preparation.directory_prep("$(DirectoryName)") 

output,map,xvecfreq  = Functionforfbm.fbm2D(imsize,powerlaw)
datadim = (imsize,imsize)
Data_preparation.write_cv_fits(datadir("sims/polaris_hera/CVI/","10pc_lag3.fits"),datadir("sims/$(DirectoryName)/","cv_simulated_$(imsize)imsize_$(powerlaw)powerlaw"),map,datadim,blank)
heatmap(map,aspect_ratio=:equal,title = "Map obtained by our algo, represent cv results. \n $(DataNameTitle)",titlefontsize=9)
savefig(plotsdir("$(DirectoryName)/","cvmap_imsize$(imsize)_powerlaw$powerlaw.png"))
savefig(plotsdir("$(DirectoryName)/","cvmap_imsize$(imsize)_powerlaw$powerlaw.pdf"))

Functionforfbm.power_spectra(map,xvecfreq,imsize)
savefig(plotsdir("$(DirectoryName)/","powerspectra_imsize$(imsize)_powerlaw$powerlaw.pdf"))

map_nomissing = Data_preparation.delete_allnotvalue(map,blank)
pdf_mine = Structure_fct.pdf_normed(map_nomissing,30,(-15,15),blank)
Graphic.cvi_pdf_norm(pdf_mine,[-7,7],[1e-5,1],0.8,3," $(DataNameTitle)","")
savefig(plotsdir("$(DirectoryName)/","pdf_normed_imsize$(imsize)_powerlaw$powerlaw.pdf"))



cvimap_allangle_alllag = Array{Union{Missing,Float64},3}(undef,imsize*imsize,192,size(Lag)[1])
cvimap_mean = Array{Union{Missing,Float64},2}(undef,size(cvimap_allangle_alllag)[1],size(Lag)[1])
cvimap_allangle = Array{Union{Missing,Float64},3}(undef,imsize,imsize,192)


cvimap_allangle_alllag = Functionforcvi.cv_increment!(map,Lag,192,cvimap_allangle_alllag,cvimap_allangle,periodic=true,diff="absolute")


struct_lag = Array{Float64}(undef,size(StcFctOrder)[1],size(Lag)[1])
struct_lag = Structure_fct.construct_fctstruct!(cvimap_allangle_alllag,StcFctOrder,Lag,100,(-1,1),blank,struct_lag)

plot(struct_lag[3,:],struct_lag[2,:],seriestype=:scatter,xaxis=:log,yaxis=:log,label="",xlabel=L"S_3(l)", ylabel=L"S_2(l)",title="Second order structure function \n $(DataNameTitle)",titlefontsize=10)
model(x,xhi) = xhi[2].*x.^xhi[1]    
FitFbm = curve_fit(model, struct_lag[3,2:end],struct_lag[2,2:end], [1.,1e-6])
FitFbm.param
plot!(struct_lag[3,2:end],FitFbm.param[2].*struct_lag[3,2:end].^FitFbm.param[1],xaxis=:log,yaxis=:log,legend=:topleft,label="Fit $(trunc(FitFbm.param[2],sigdigits=3))*Lag^$(trunc(FitFbm.param[1],sigdigits=3))")
savefig(plotsdir("$(DirectoryName)/","stcfct2_imsize$(imsize)_powerlaw$powerlaw.pdf"))



