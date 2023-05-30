###################################################################
# First script used to test the efficacity of the code.
###################################################################
include("../src/Data_preparation.jl") # Read and write fits
include("../src/Functionforcvi.jl")   # Calculations of CVI
include("../src/Graphic.jl")
include("../src/Data_analysis.jl")

using FITSIO                                                    # Read Fits
using MultivariateStats, Statistics, StatsBase, Distributions   # Statistic
using Profile, BenchmarkTools                                   # Benchmark
using Mmap, DelimitedFiles                                      # Read and write files (.bin and .txt)
using ShiftedArrays                                             # Shifted Arrays for CVI
using Plots
using PairPlots
using Printf
using KernelDensity
using Measures
using LaTeXStrings

using .Functionforcvi
using .Data_preparation
using .Graphic
using .Data_analysis

gr()

FITSPATH      = "/home/delcamps/Data/Simulated/OnlyGaussian/"
PATHTOSAVE    = "/home/delcamps/Prog/CVI_alone/"
DIRECTORYNAME = "OnlyGaussian"
UNITVELOCITY  = "m/s"

# Prepare directories where plots and data will be saved.
(isdir("$(PATHTOSAVE)/Plots"))==0  && mkdir("$(PATHTOSAVE)/Plots")
(isdir("$(PATHTOSAVE)/Data"))==0   && mkdir("$(PATHTOSAVE)/Data")



cubebase,HEAD,DATADIMENSION = Data_preparation.read_fits_pp("$(FITSPATH)"*"onegauslarge_val.fits")
cubebase = reshape(cubebase,DATADIMENSION[1]*DATADIMENSION[2])

snr = Vector{Float64}(undef,6)
cvmaps = Array{Float64}(undef,DATADIMENSION[1]*DATADIMENSION[2],6)

names = ["0.00","0.02","0.04","0.06","0.08","0.10"]
namesontwo = ["0.02","0.06","0.10"]
for nois=1:6
    global namesontwo
    cube,VELOCITYVECTOR = Data_preparation.read_fits_ppv("$(FITSPATH)"*"onegauslarge_$(names[nois]).fits",UNITVELOCITY ; check=false)[1:2]
    cvmap = Functionforcvi.moment_one_field(cube,VELOCITYVECTOR.*1000) # Calculate the first velocity moment order on cube
    snr[nois] = Data_analysis.snr_allfield(reshape(cube,DATADIMENSION[1]*DATADIMENSION[2],size(VELOCITYVECTOR)[1]),(98:150,98:150))
    cvmaps[:,nois] = reshape(cvmap,DATADIMENSION[1]*DATADIMENSION[2],1)
end


a = convert(Vector{Float64},cvmaps[:,1])
b = convert(Vector{Float64},cvmaps[:,2])
c = convert(Vector{Float64},cvmaps[:,3])
d = convert(Vector{Float64},cvmaps[:,4])
e = convert(Vector{Float64},cvmaps[:,5])
f = convert(Vector{Float64},cvmaps[:,6])
#Graphic.corner_cvimap(cubebase,cvmaps,DATADIMENSION,"label1",namesontwo)



asnr = @sprintf("%5.2e",snr[1])
bsnr = @sprintf("%5.2e",snr[2])
csnr = @sprintf("%5.2e",snr[3])
dsnr = @sprintf("%5.2e",snr[4])
esnr = @sprintf("%5.2e",snr[5])
fsnr = @sprintf("%5.2e",snr[6])


#=
dens_basea = kde((cubebase,a))
p = plot(layout=grid(7,7),legend = false,link=:x,size=(1300,1300))#, size=(500,300),grid=:false,dpi=500) 
p = histogram!(p[1],cubebase[:],fillcolor=:white,xlims=[-6,10],nbins=14,aspect_ratio=1/59,ylims=[-10,1000 ])
for ix=2:7
    global p
    p = plot!(p[ix],showaxis=:hide,foreground_color_text=:white,grid=:false)
end
p = plot!(p[8],cubebase,a,aspect_ratio=:equal,xlims=[-6,10],ylims=[-6,10],ylabel=L"\sigma =0.81E+09",clims=(0.034,0.035),color=:blue)
p = histogram!(p[9],a[:],orientation=:horizontal, fillcolor=:white,ylims=[-6,10],nbins=14,aspect_ratio=59,xlims=[-10,1000])
for ix=10:14
    p = plot!(p[ix],showaxis=:hide,foreground_color_text=:white,grid=:false)
end
dens_baseb = kde((cubebase,b))
dens_ab = kde((a,b))
p = plot!(p[15],dens_baseb,ylabel=L"\sigma =0.30E+02")
p = plot!(p[16],dens_ab,xlims=[-6,10],ylims=[-6,10])
p = histogram(p[17],b[:],orientation=:horizontal, fillcolor=:white,ylims=[-6,10],nbins=14,aspect_ratio=59,xlims=[-10,1000])
for ix=18:21
    p = plot!(p[ix],showaxis=:hide,foreground_color_text=:white,grid=:false)
end

dens_basec = kde((cubebase,c))
dens_ac = kde((a,c))
dens_bc = kde((b,c))
p = plot!(p[22],dens_baseb,ylabel=L"\sigma =0.16E+02")
p = plot!(p[23],dens_ac,xlims=[-6,10],ylims=[-6,10])
p = plot!(p[24],dens_bc,xlims=[-6,10],ylims=[-6,10])
p = histogram(p[25],c[:],orientation=:horizontal, fillcolor=:white,ylims=[-6,10],nbins=14,aspect_ratio=59,xlims=[-10,1000])

for ix=26:28
    p = plot!(p[ix],showaxis=:hide,foreground_color_text=:white,grid=:false)
end
dens_based = kde((cubebase,d))
dens_ad = kde((a,d))
dens_bd = kde((b,d))
dens_cd = kde((c,d))
p = plot!(p[29],dens_based,ylabel=L"\sigma =0.12E+02")
p = plot!(p[30],dens_ad,xlims=[-6,10],ylims=[-6,10])
p = plot!(p[31],dens_bd,xlims=[-6,10],ylims=[-6,10])
p = plot!(p[32],dens_cd,xlims=[-6,10],ylims=[-6,10])
p = histogram(p[33],d[:],orientation=:horizontal, fillcolor=:white,ylims=[-6,10],nbins=14,aspect_ratio=59,xlims=[-10,1000])
for ix=34:35
    p = plot!(p[ix],showaxis=:hide,foreground_color_text=:white,grid=:false)
end

dens_basee = kde((cubebase,e))
dens_ae = kde((a,e))
dens_be = kde((b,e))
dens_ce = kde((c,e))
dens_de = kde((d,e))
p = plot!(p[36],dens_basee,ylabel=L"\sigma =0.92E+01")
p = plot!(p[37],dens_ae,xlims=[-6,10],ylims=[-6,10])
p = plot!(p[38],dens_be,xlims=[-6,10],ylims=[-6,10])
p = plot!(p[39],dens_ce,xlims=[-6,10],ylims=[-6,10])
p = plot!(p[40],dens_de,xlims=[-6,10],ylims=[-6,10])
p = histogram(p[41],e[:],orientation=:horizontal, fillcolor=:white,ylims=[-6,10],nbins=14,aspect_ratio=59,xlims=[-10,1000])
p = plot!(p[42],showaxis=:hide,foreground_color_text=:white,grid=:false)


dens_basef = kde((cubebase,f))
dens_af = kde((a,f))
dens_bf = kde((b,f))
dens_df = kde((d,f))
dens_cf = kde((c,f))
dens_df = kde((d,f))
dens_ef = kde((e,f))
p = plot!(p[43],dens_basef,ylabel=L"\sigma =0.78E+01")
p = plot!(p[44],dens_af,xlims=[-6,10],ylims=[-6,10])
p = plot!(p[45],dens_bf,xlims=[-6,10],ylims=[-6,10])
p = plot!(p[46],dens_cf,xlims=[-6,10],ylims=[-6,10])
p = plot!(p[47],dens_df,xlims=[-6,10],ylims=[-6,10])
p = plot!(p[48],dens_ef,xlims=[-6,10],ylims=[-6,10])
p = histogram(p[49],f[:],orientation=:horizontal, fillcolor=:white,ylims=[-6,10],nbins=14,aspect_ratio=59,xlims=[-10,1900])

display(p)
=#


p = plot(layout=grid(2,4),legend = false,link=:both,size=(870,450),dpi=1000)
p = plot!(p[1],showaxis=:hide,foreground_color_text=:white,grid=:false)
p = histogram!(p[2],a[:],fillcolor=:white,xlims=[-6,10],nbins=14,aspect_ratio=1/100,ylims=[-10,1600],xlabel=L"snr=0.81E+09")
p = histogram!(p[3],b[:],fillcolor=:white,xlims=[-6,10],nbins=14,aspect_ratio=1/100,ylims=[-10,1600],xlabel=L"snr=0.30E+02")
p = histogram!(p[4],c[:],fillcolor=:white,xlims=[-6,10],nbins=14,aspect_ratio=1/100,ylims=[-10,1600],xlabel=L"snr=0.16E+02")
p = histogram!(p[5],cubebase[:],orientation=:horizontal, fillcolor=:white,ylims=[-6,10],nbins=14,aspect_ratio=120,xlims=[-10,1600],xflip=:true,ylabel="Expected Values")
p = plot!(p[6],cubebase,a,seriestype=:scatter,markersize=0.5,markershape=:+,alpha=0.5,xlims=[-6,10],ylims=[-6,10],color=:black,xlabel="CV values (a.u.)")
p = plot!(p[7],cubebase,b,seriestype=:scatter,markersize=0.5,markershape=:+,alpha=0.5,xlims=[-6,10],ylims=[-6,10],color=:black,xlabel="CV values (a.u.)")
p = plot!(p[8],cubebase,c,seriestype=:scatter,markersize=0.5,markershape=:+,alpha=0.5,xlims=[-6,10],ylims=[-6,10],color=:black,xlabel="CV values (a.u.)")
for ix=6:8
    global p
    p = plot!(p[ix],[-5,10],[-5,10],seriestype=:line,color=:blue)
end
display(p)


savefig("$(PATHTOSAVE)/Plots/3first_histo.png")

p = plot(layout=grid(2,4),legend = false,link=:both,size=(870,450),dpi=1000) #link=:y,
p = plot!(p[1],showaxis=:hide,foreground_color_text=:white,grid=:false)
p = histogram!(p[2],d[:],fillcolor=:white,xlims=[-6,10],nbins=14,aspect_ratio=1/100,ylims=[-10,1600],xlabel=L"snr = 0.12E+02")
p = histogram!(p[3],e[:],fillcolor=:white,xlims=[-6,10],nbins=14,aspect_ratio=1/100,ylims=[-10,1600],xlabel=L"snr = 0.92E+01")
p = histogram!(p[4],f[:],fillcolor=:white,xlims=[-6,10],nbins=14,aspect_ratio=1/100,ylims=[-10,1600],xlabel=L"snr = 0.78E+01")
p = histogram!(p[5],cubebase[:],orientation=:horizontal, fillcolor=:white,ylims=[-6,10],nbins=14,aspect_ratio=120,xlims=[-10,1600],xflip=:true,ylabel="Expected Values")
p = plot!(p[6],cubebase,d,seriestype=:scatter,markersize=0.5,markershape=:+,alpha=0.5,xlims=[-6,10],ylims=[-6,10],color=:black,xlabel="CV values (a.u.)")
p = plot!(p[7],cubebase,e,seriestype=:scatter,markersize=0.5,markershape=:+,alpha=0.5,xlims=[-6,10],ylims=[-6,10],color=:black,xlabel="CV values (a.u.)")
p = plot!(p[8],cubebase,f,seriestype=:scatter,markersize=0.5,markershape=:+,alpha=0.5,xlims=[-6,10],ylims=[-6,10],color=:black,xlabel="CV values (a.u.)")
for ix=6:8
    global p
    p = plot!(p[ix],[-5,10],[-5,10],seriestype=:line,color=:blue)
end
display(p)
savefig("$(PATHTOSAVE)/Plots/3last_histo.png")
#=
data = (;cubebase,a,b,c,d,e,f)
corner(data,["Real values",asnr,bsnr,csnr,dsnr,esnr,fsnr],plotcontours=false,scatter_kwargs=(;xlims=[-6.,3],ylims=[-6,3],color=:red),contour_kwargs=(;xlims=[-6,3],ylims=[-6,3]),hist_kwargs=(;xlims=[-10,2],ylims=[0,300]),hist2d_kwargs=(;color=:inferno),title="Corner plot of CV calculated on simple ppv cube. Names are mean snr values of the cube.",titlefontsize=10)#,contour_kwargs=(xlims=[-10.,10.],ylims=[-10.,10.]))#,hist2d_kwargs=(;xlims=[-10,10]),hist_kwargs=(;xlims=[-10,10]))#,contour_kwargs=(xlims=[-8,8],ylims=[-8,8]))


savefig("$(PATHTOSAVE)/Plots/cornerplotlarge.pdf")
=#