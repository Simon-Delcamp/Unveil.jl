###################################################################
# Second script used to test the efficacity of the code.
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
using PrettyTables
using .Functionforcvi
using .Data_preparation
using .Graphic
using .Data_analysis

gr()

FITSPATH_CV      = "/home/delcamps/Prog/CVI_light/Data/WingsMove/"
FITSPATH_RAW     = "/home/delcamps/Data/Simulated/WingsMove/"
FITSPATH_PC      = "$FITSPATH_CV"
PATHTOSAVE       = "/home/delcamps/Prog/CVI_light/"
UNITVELOCITY     = "m/s"

# Prepare directories where plots and data will be saved.
(isdir("$(PATHTOSAVE)/"))==0  && mkdir("$(PATHTOSAVE)/")
(isdir("$(PATHTOSAVE)/Data"))==0   && mkdir("$(PATHTOSAVE)/Data/")



cubebase,HEAD,DATADIMENSION = Data_preparation.read_fits_pp("$(FITSPATH_CV)/0_00/CV/cv_multlag_rawPC.fits")
cubebase = reshape(cubebase,DATADIMENSION[1]*DATADIMENSION[2])

snr = Vector{Float64}(undef,12)
cvmaps = Array{Float64}(undef,DATADIMENSION[1]*DATADIMENSION[2],12)
uncertcvmaps = Array{Float64}(undef,DATADIMENSION[1]*DATADIMENSION[2],12)

names = ["0_00","0_02","0_04","0_06","0_08","0_10","0_00","0_02_nois","0_04_nois","0_06_nois","0_08_nois","0_10_nois"]
namesraw = ["0_00","0_02","0_04","0_06","0_08","0_10","0_00","0_02","0_04","0_06","0_08","0_10"]

namesrawcube = ["0.00","0.02","0.04","0.06","0.08","0.10","0.00","0.02","0.04","0.06","0.08","0.10"]

noisecanals = [200:256,200:256]

function uncertaintyofmoment(k,S,CV,uncertCV)
    meanCV = moment(CV,1,0)
    #uncert = k*sqrt(S+1)*sum(uncertCV.*(CV.-meanCV).^(k-1))*S^(-3/2)
    uncertA = k^2. *uncertCV.^2 .*(CV.-meanCV).^(2*k-2)
    uncertB = 1+(sum(uncertCV))^2/S^2
    uncert = sqrt(sum(uncertA*uncertB)*S^-2)
    return(uncert)
end

for nois=1:12
    global names,namesraw,snr,cvmap,cvmaps
    if nois==1
        plac = "rawPC"
    elseif nois==7
        plac = "rawPC"
    elseif nois>=8
        plac = "rawPC"
    else
        plac = "recon6PC"
    end
    cvmap,HEAD,DATADIMENSION = Data_preparation.read_fits_pp("$(FITSPATH_CV)/$(namesraw[nois])/CV/cv_multlag_$(plac).fits")
    
    #cvmap = Functionforcvi.moment_one_field(cube,VELOCITYVECTOR.*1000) # Calculate the first velocity moment order on cube
    cube,VELOCITYVECTOR = Data_preparation.read_fits_ppv("$(FITSPATH_RAW)"*"fake_$(namesrawcube[nois]).fits",UNITVELOCITY ; check=false)[1:2]
    cube = reshape(cube,DATADIMENSION[1]*DATADIMENSION[2],size(VELOCITYVECTOR)[1])
    uncertcvmap = Data_analysis.rms_analytic_field(cube,VELOCITYVECTOR,noisecanals)
    snr[nois] = Data_analysis.snr_allfield(cube,(180:250,180:250))
    cvmaps[:,nois] = reshape(cvmap,DATADIMENSION[1]*DATADIMENSION[2],1)
    uncertcvmaps[:,nois] = reshape(uncertcvmap,DATADIMENSION[1]*DATADIMENSION[2],1)
end


a = convert(Vector{Float64},cvmaps[:,1])
b = convert(Vector{Float64},cvmaps[:,2])
c = convert(Vector{Float64},cvmaps[:,3])
d = convert(Vector{Float64},cvmaps[:,4])
e = convert(Vector{Float64},cvmaps[:,5])
f = convert(Vector{Float64},cvmaps[:,6])
a_noisy = convert(Vector{Float64},cvmaps[:,7])
b_noisy = convert(Vector{Float64},cvmaps[:,8])
c_noisy = convert(Vector{Float64},cvmaps[:,9])
d_noisy = convert(Vector{Float64},cvmaps[:,10])
e_noisy = convert(Vector{Float64},cvmaps[:,11])
f_noisy = convert(Vector{Float64},cvmaps[:,12])

#Graphic.corner_cvimap(cubebase,cvmaps,DATADIMENSION,"label1",namesontwo)


asnr = @sprintf("%5.2e",snr[1])
bsnr = @sprintf("%5.2e",snr[2])
csnr = @sprintf("%5.2e",snr[3])
dsnr = @sprintf("%5.2e",snr[4])
esnr = @sprintf("%5.2e",snr[5])
fsnr = @sprintf("%5.2e",snr[6])
a_noisysnr=@sprintf("%.1f",snr[7])
b_noisysnr=@sprintf("%.1f",snr[8])
c_noisysnr=@sprintf("%.1f",snr[9])
d_noisysnr=@sprintf("%.1f",snr[10])
e_noisysnr=@sprintf("%.1f",snr[11])
f_noisysnr=@sprintf("%.1f",snr[12])

mu_uncertmoment  = Array{Float64}(undef,5)
sig_uncertmoment = Array{Float64}(undef,5)
gam_uncertmoment = Array{Float64}(undef,5)
kap_uncertmoment = Array{Float64}(undef,5)
for ix=2:6
    mu_uncertmoment[ix-1] =  1/(DATADIMENSION[1]*DATADIMENSION[2])*sqrt(sum(uncertcvmaps[:,ix].^2))
    sig_uncertmoment[ix-1] = uncertaintyofmoment(2,DATADIMENSION[1]*DATADIMENSION[2],cvmaps[:,ix],uncertcvmaps[:,ix])
    gam_uncertmoment[ix-1] = uncertaintyofmoment(3,DATADIMENSION[1]*DATADIMENSION[2],cvmaps[:,ix],uncertcvmaps[:,ix])
    kap_uncertmoment[ix-1] = uncertaintyofmoment(4,DATADIMENSION[1]*DATADIMENSION[2],cvmaps[:,ix],uncertcvmaps[:,ix])
end

mu_uncertmomentNOISE  = Array{Float64}(undef,5)
sig_uncertmomentNOISE = Array{Float64}(undef,5)
gam_uncertmomentNOISE = Array{Float64}(undef,5)
kap_uncertmomentNOISE = Array{Float64}(undef,5)

for ix=2:6
    mu_uncertmomentNOISE[ix-1] = 1/(DATADIMENSION[1]*DATADIMENSION[2])*sqrt(sum(uncertcvmaps[:,ix+6].^2))
    sig_uncertmomentNOISE[ix-1] = uncertaintyofmoment(2,DATADIMENSION[1]*DATADIMENSION[2],cvmaps[:,ix+6],uncertcvmaps[:,ix+6])
    gam_uncertmomentNOISE[ix-1] = uncertaintyofmoment(3,DATADIMENSION[1]*DATADIMENSION[2],cvmaps[:,ix+6],uncertcvmaps[:,ix+6])
    kap_uncertmomentNOISE[ix-1] = uncertaintyofmoment(4,DATADIMENSION[1]*DATADIMENSION[2],cvmaps[:,ix+6],uncertcvmaps[:,ix+6])
end

mu_uncertmomentBASE = 1/(DATADIMENSION[1]*DATADIMENSION[2])*sqrt(sum(uncertcvmaps[:,1].^2)) 
sig_uncertmomentBASE = uncertaintyofmoment(2,DATADIMENSION[1]*DATADIMENSION[2],cvmaps[:,1],uncertcvmaps[:,1])
gam_uncertmomentBASE = uncertaintyofmoment(3,DATADIMENSION[1]*DATADIMENSION[2],cvmaps[:,1],uncertcvmaps[:,1])
kap_uncertmomentBASE = uncertaintyofmoment(4,DATADIMENSION[1]*DATADIMENSION[2],cvmaps[:,1],uncertcvmaps[:,1])


mu = [moment(cubebase[:],1,0),moment(b[:],1,0),moment(c[:],1,0),moment(d[:],1,0),moment(e[:],1,0),moment(f[:],1,0),moment(b_noisy[:],1,0),moment(c_noisy[:],1,0),moment(d_noisy[:],1,0),moment(e_noisy[:],1,0),moment(f_noisy[:],1,0)] 

sig = [moment(cubebase[:],2),moment(b[:],2),moment(c[:],2),moment(d[:],2),moment(e[:],2),moment(f[:],2),moment(b_noisy[:],2),moment(c_noisy[:],2),moment(d_noisy[:],2),moment(e_noisy[:],2),moment(f_noisy[:],2)] 

gam = [moment(cubebase[:],3),moment(b[:],3),moment(c[:],3),moment(d[:],3),moment(e[:],3),moment(f[:],3),moment(b_noisy[:],3),moment(c_noisy[:],3),moment(d_noisy[:],3),moment(e_noisy[:],3),moment(f_noisy[:],3)] 

kap = [moment(cubebase[:],4),moment(b[:],4),moment(c[:],4),moment(d[:],4),moment(e[:],4),moment(f[:],4),moment(b_noisy[:],4),moment(c_noisy[:],4),moment(d_noisy[:],4),moment(e_noisy[:],4),moment(f_noisy[:],4)] 
 


xaxis = [b_noisysnr,c_noisysnr,d_noisysnr,e_noisysnr,f_noisysnr]
p = plot(layout=grid(2,2),legend = false,link=:x,size=(700,500),dpi=800,xlabelfontsize=10,ylabelfontsize=12)

p = plot!(p[1],xaxis,mu[2:6],seriestype=:scatter,markershape=:x,markersize=3,legend=legend=:bottomleft,label="PCA",yerror=mu_uncertmoment,msc=:blue)
p = plot!(p[1],xaxis,mu[7:11],seriestype=:scatter,markershape=:x,markersize=3,label="No PCA",yerror=mu_uncertmomentNOISE,msc=:peru)
p = plot!(p[1],[xaxis[1],xaxis[5]],[mu[1],mu[1]],label="No noise",ribbon=mu_uncertmomentBASE)

p = plot!(p[2],xaxis,sig[2:6],seriestype=:scatter,markershape=:x,markersize=3,yerror=sig_uncertmoment,msc=:blue)
p = plot!(p[2],xaxis,sig[7:11],seriestype=:scatter,markershape=:x,markersize=3,yerror=sig_uncertmomentNOISE,msc=:peru)
p = plot!(p[2],[xaxis[1],xaxis[5]],[sig[1],sig[1]],ribbon=sig_uncertmomentBASE)

p = plot!(p[3],xaxis,gam[2:6],seriestype=:scatter,markershape=:x,markersize=3,yerror=gam_uncertmoment,msc=:blue)
p = plot!(p[3],xaxis,gam[7:11],seriestype=:scatter,markershape=:x,markersize=3,yerror=gam_uncertmomentNOISE,msc=:peru)
p = plot!(p[3],[xaxis[1],xaxis[5]],[gam[1],gam[1]],ribbon=gam_uncertmomentBASE)

p = plot!(p[4],xaxis,kap[2:6],seriestype=:scatter,markershape=:x,markersize=3,yerror=kap_uncertmoment,msc=:blue)
p = plot!(p[4],xaxis,kap[7:11],seriestype=:scatter,markershape=:x,markersize=3,yerror=kap_uncertmomentNOISE,msc=:peru)
p = plot!(p[4],[xaxis[1],xaxis[5]],[kap[1],kap[1]],ribbon=kap_uncertmomentBASE)

p = plot!(p[4],xlabel="S/N value of noisy cube")
p = plot!(p[3],xlabel="S/N value of noisy cube")

p = plot!(p[1],yaxis=L"\mu")
p = plot!(p[2],yaxis=L"\sigma")
p = plot!(p[3],yaxis=L"\gamma")
p = plot!(p[4],yaxis=L"\kappa")
display(p)
savefig("$(PATHTOSAVE)/Plots/comparmoments.png")


header = [L"\mu" L"\sigma" L"\gamma" L"\kappa";]
data = hcat(mu,sig,gam,kap)
pretty_table(data,header,header_crayon = crayon"yellow bold", tf = tf_mysql)

open("$(PATHTOSAVE)/Data/sig_mu.txt", "w") do f
    pretty_table(f,data,header,header_crayon = crayon"yellow bold", tf = tf_mysql)
end



p = plot(layout=grid(1,3),legend = false,size=(800,300),dpi=1000,leftmargins=0.25cm,bottommargins=0.1cm,xticks=-3.:0.5:-0,xminorticks=0.1,labelfontsize=10,titlefontsize=10)
limx = [-3,0]
limy = [-10,900]
rap  = 300
labx = "CV values"
laby = "Number of iterations"

p = histogram!(p[1],cubebase[:],aspect_ratio=1/rap,xlims=limx,ylims=limy,color=:green,alpha=0.3)
p = histogram!(p[1],cvmaps[:,12],aspect_ratio=1/rap,xlims=limx,ylims=limy,color=:peru,alpha=0.3)
p = plot!(p[1],ylabel="Number of CV",xlabel="CV values")
p = annotate!(p[1],[limx[1]+(limx[2]-limx[1])/5],[limy[end]-limy[end]/8],text("No noise",:green,10))
p = annotate!(p[1],[limx[2]-(limx[2]-limx[1])/5],[limy[end]-limy[end]/8],text("No PCA",:peru,10))

p = histogram!(p[2],cvmaps[:,12],aspect_ratio=1/rap,xlims=limx,ylims=limy,color=:peru,alpha=0.3)
p = histogram!(p[2],cvmaps[:,6],aspect_ratio=1/rap,xlims=limx,ylims=limy,color=:blue,alpha=0.2)
p = plot!(p[2],xlabel="CV values")
p = annotate!(p[2],[limx[1]+(limx[2]-limx[1])/5],[limy[end]-limy[end]/8],text("PCA",:blue,10))
p = annotate!(p[2],[limx[2]-(limx[2]-limx[1])/5],[limy[end]-limy[end]/8],text("No PCA",:peru,10))


p = histogram!(p[3],cubebase[:],aspect_ratio=1/rap,xlims=limx,ylims=limy,color=:green,alpha=0.3)
p = histogram!(p[3],cvmaps[:,6],aspect_ratio=1/rap,xlims=limx,ylims=limy,color=:blue,alpha=0.2)
p = plot!(p[3],xlabel="CV values")
p = annotate!(p[3],[limx[1]+(limx[2]-limx[1])/5],[limy[end]-limy[end]/8],text("PCA",:blue,10))
p = annotate!(p[3],[limx[2]-(limx[2]-limx[1])/5],[limy[end]-limy[end]/8],text("No noise",:green,10))
savefig("$(PATHTOSAVE)/Plots/comparhisto.png")

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

#=
# PRODUCE COMPARISON BETWEEN CV WITH PC AND WITHOUT
histogram(cubebase,alpha=0.3,nbins=10)
histogram!(f,alpha=0.3,nbins=15)
histogram!(f_noisy,alpha=0.3,nbins=50)
limx = [-2.6,0]
limy = [-10,1000]
rap  = 500
labx = "CV values"
laby = "Number of iterations"
p = plot(layout=grid(2,3),legend = false,size=(900,500),dpi=1000,leftmargins=0.25cm,bottommargins=0.1cm,xticks=-2.:0.5:-0,xminorticks=0.1)

p = histogram!(p[1],a,xlims=limx,nbins=20,aspect_ratio=1/rap,ylims=limy,title=L"snr=0.81E+09",alpha=0.3)
p = histogram!(p[1],a_noisy,xlims=limx,nbins=20,aspect_ratio=1/rap,ylims=limy,alpha=0.2)
p = annotate!(p[1],[limx[1]+(limx[2]-limx[1])/5],[limy[end]-limy[end]/5],text(L"\sigma = "*"$(trunc(sqrt(std(a[:])); sigdigits=2))"*"\n"*L"\mu = "*"$(trunc(mean(a[:]); sigdigits=3))",:blue,8))
p = annotate!(p[1],[limx[2]-(limx[2]-limx[1])/6],[limy[end]-limy[end]/5],text(L"\sigma = "*"$(trunc(sqrt(std(a_noisy[:])); sigdigits=2))"*"\n"*L"\mu = "*"$(trunc(mean(a_noisy[:]); sigdigits=3))",:peru,8))

p = histogram!(p[2],b,xlims=limx,nbins=20,aspect_ratio=1/rap,ylims=limy,title=L"snr=0.30E+02",alpha=0.3)
p = histogram!(p[2],b_noisy,xlims=limx,nbins=20,aspect_ratio=1/rap,ylims=limy,alpha=0.2)
p = annotate!(p[2],[limx[1]+(limx[2]-limx[1])/5],[limy[end]-limy[end]/5],text(L"\sigma = "*"$(trunc(sqrt(std(b[:])); sigdigits=2))"*"\n"*L"\mu = "*"$(trunc(mean(b[:]); sigdigits=3))",:blue,8))
p = annotate!(p[2],[limx[2]-(limx[2]-limx[1])/6],[limy[end]-limy[end]/5],text(L"\sigma = "*"$(trunc(sqrt(std(b_noisy[:])); sigdigits=2))"*"\n"*L"\mu = "*"$(trunc(mean(b_noisy[:]); sigdigits=3))",:peru,8))

p = histogram!(p[3],c,xlims=limx,nbins=20,aspect_ratio=1/rap,ylims=limy,title=L"snr=0.16E+02",alpha=0.3)
p = histogram!(p[3],c_noisy,xlims=limx,nbins=20,aspect_ratio=1/rap,ylims=limy,alpha=0.2)
p = annotate!(p[3],[limx[1]+(limx[2]-limx[1])/5],[limy[end]-limy[end]/5],text(L"\sigma = "*"$(trunc(sqrt(std(c[:])); sigdigits=2))"*"\n"*L"\mu = "*"$(trunc(mean(c[:]); sigdigits=3))",:blue,8))
p = annotate!(p[3],[limx[2]-(limx[2]-limx[1])/6],[limy[end]-limy[end]/5],text(L"\sigma = "*"$(trunc(sqrt(std(c_noisy[:])); sigdigits=2))"*"\n"*L"\mu = "*"$(trunc(mean(c_noisy[:]); sigdigits=3))",:peru,8))

p = histogram!(p[4],d,xlims=limx,nbins=20,aspect_ratio=1/rap,ylims=limy,title=L"snr=0.12E+02",alpha=0.3)
p = histogram!(p[4],d_noisy,xlims=limx,nbins=30,aspect_ratio=1/rap,ylims=limy,alpha=0.2)
p = annotate!(p[4],[limx[1]+(limx[2]-limx[1])/5],[limy[end]-limy[end]/5],text(L"\sigma = "*"$(trunc(sqrt(std(d[:])); sigdigits=2))"*"\n"*L"\mu = "*"$(trunc(mean(d[:]); sigdigits=3))",:blue,8))
p = annotate!(p[4],[limx[2]-(limx[2]-limx[1])/6],[limy[end]-limy[end]/5],text(L"\sigma = "*"$(trunc(sqrt(std(d_noisy[:])); sigdigits=2))"*"\n"*L"\mu = "*"$(trunc(mean(d_noisy[:]); sigdigits=3))",:peru,8))

p = histogram!(p[5],e,xlims=limx,nbins=20,aspect_ratio=1/rap,ylims=limy,title=L"snr=0.92E+01",alpha=0.3)
p = histogram!(p[5],e_noisy,xlims=limx,nbins=30,aspect_ratio=1/rap,ylims=limy,alpha=0.2)
p = annotate!(p[5],[limx[1]+(limx[2]-limx[1])/5],[limy[end]-limy[end]/5],text(L"\sigma = "*"$(trunc(sqrt(std(e[:])); sigdigits=2))"*"\n"*L"\mu = "*"$(trunc(mean(e[:]); sigdigits=3))",:blue,8))
p = annotate!(p[5],[limx[2]-(limx[2]-limx[1])/6],[limy[end]-limy[end]/5],text(L"\sigma = "*"$(trunc(sqrt(std(e_noisy[:])); sigdigits=2))"*"\n"*L"\mu = "*"$(trunc(mean(e_noisy[:]); sigdigits=3))",:peru,8))

p = histogram!(p[6],f,xlims=limx,nbins=20,aspect_ratio=1/rap,ylims=limy,title=L"snr=0.78E+01",alpha=0.3)
p = histogram!(p[6],f_noisy,xlims=limx,nbins=30,aspect_ratio=1/rap,ylims=limy,alpha=0.2)
p = annotate!(p[6],[limx[1]+(limx[2]-limx[1])/5],[limy[end]-limy[end]/5],text(L"\sigma = "*"$(trunc(sqrt(std(f[:])); sigdigits=2))"*"\n"*L"\mu = "*"$(trunc(mean(f[:]); sigdigits=3))",:blue,8))
p = annotate!(p[6],[limx[2]-(limx[2]-limx[1])/6],[limy[end]-limy[end]/5],text(L"\sigma = "*"$(trunc(sqrt(std(f_noisy[:])); sigdigits=2))"*"\n"*L"\mu = "*"$(trunc(mean(f_noisy[:]); sigdigits=3))",:peru,8))

for ix=1:6
    global p
# p = histogram!(p[ix],cubebase,fillcolor=:green,xlims=limx,ylims=limy,nbins=10,alpha=0.25)
    
end
p = histogram!(p[1],ylabel=laby)
p = histogram!(p[4],ylabel=laby)
for ix=4:6
    global p
    p = histogram!(p[ix],xlabel=labx)
end
display(p)

savefig("$(PATHTOSAVE)/Plots/histo_compar_CVPCAandNOT.png")


# PRODUCE COMPARISON BETWEEN CV PCA AND EXPECTED
histogram(cubebase,alpha=0.3,nbins=10)
histogram!(f,alpha=0.3,nbins=15)
histogram!(f_noisy,alpha=0.3,nbins=50)
limx = [-2.6,0]
limy = [-10,1000]
rap  = 500
labx = "CV values"
laby = "Number of iterations"
p = plot(layout=grid(2,3),legend = false,size=(900,500),dpi=1000,leftmargins=0.25cm,bottommargins=0.1cm,xticks=-2.:0.5:-1,xminorticks=0.25)

p = histogram!(p[1],a,xlims=limx,nbins=20,aspect_ratio=1/rap,ylims=limy,title=L"snr=0.81E+09",alpha=0.3)
p = annotate!(p[1],[limx[1]+(limx[2]-limx[1])/5],[limy[end]-limy[end]/5],text(L"\sigma = "*"$(trunc(sqrt(std(a[:])); sigdigits=2))"*"\n"*L"\mu = "*"$(trunc(mean(a[:]); sigdigits=3))",:blue,8))

p = histogram!(p[2],b,xlims=limx,nbins=20,aspect_ratio=1/rap,ylims=limy,title=L"snr=0.30E+02",alpha=0.3)
p = annotate!(p[2],[limx[1]+(limx[2]-limx[1])/5],[limy[end]-limy[end]/5],text(L"\sigma = "*"$(trunc(sqrt(std(b[:])); sigdigits=2))"*"\n"*L"\mu = "*"$(trunc(mean(b[:]); sigdigits=3))",:blue,8))

p = histogram!(p[3],c,xlims=limx,nbins=20,aspect_ratio=1/rap,ylims=limy,title=L"snr=0.16E+02",alpha=0.3)
p = annotate!(p[3],[limx[1]+(limx[2]-limx[1])/5],[limy[end]-limy[end]/5],text(L"\sigma = "*"$(trunc(sqrt(std(c[:])); sigdigits=2))"*"\n"*L"\mu = "*"$(trunc(mean(c[:]); sigdigits=3))",:blue,8))

p = histogram!(p[4],d,xlims=limx,nbins=20,aspect_ratio=1/rap,ylims=limy,title=L"snr=0.12E+02",alpha=0.3)
p = annotate!(p[4],[limx[1]+(limx[2]-limx[1])/5],[limy[end]-limy[end]/5],text(L"\sigma = "*"$(trunc(sqrt(std(d[:])); sigdigits=2))"*"\n"*L"\mu = "*"$(trunc(mean(d[:]); sigdigits=3))",:blue,8))

p = histogram!(p[5],e,xlims=limx,nbins=20,aspect_ratio=1/rap,ylims=limy,title=L"snr=0.92E+01",alpha=0.3)
p = annotate!(p[5],[limx[1]+(limx[2]-limx[1])/5],[limy[end]-limy[end]/5],text(L"\sigma = "*"$(trunc(sqrt(std(e[:])); sigdigits=2))"*"\n"*L"\mu = "*"$(trunc(mean(e[:]); sigdigits=3))",:blue,8))

p = histogram!(p[6],f,xlims=limx,nbins=20,aspect_ratio=1/rap,ylims=limy,title=L"snr=0.78E+01",alpha=0.3)
p = annotate!(p[6],[limx[1]+(limx[2]-limx[1])/5],[limy[end]-limy[end]/5],text(L"\sigma = "*"$(trunc(sqrt(std(f[:])); sigdigits=2))"*"\n"*L"\mu = "*"$(trunc(mean(f[:]); sigdigits=3))",:blue,8))

for ix=1:6
    global p
   p = histogram!(p[ix],cubebase,fillcolor=:green,xlims=limx,ylims=limy,nbins=20,alpha=0.25)
end
p = histogram!(p[1],ylabel=laby)
p = histogram!(p[4],ylabel=laby)
for ix=4:6
    global p
    p = histogram!(p[ix],xlabel=labx)
end
display(p)

savefig("$(PATHTOSAVE)/Plots/histo_compar_CVPCAandEXPECTED.png")

dif = reshape(a.-f,100,100)
heatmap(dif,xlims=[0,100],ylims=[0,100],c=:viridis,clims=(-0.2,0.2),xlabel="Position (a.u.)",ylabel="Position (a.u)",aspect_ratio=:equal,dpi=1000)
savefig("$(PATHTOSAVE)/Plots/diff_expect-0-10.pdf")


# All compar but with 3 cube 
#histogram(cubebase,alpha=0.3,nbins=10)
#histogram!(f,alpha=0.3,nbins=15)
#histogram!(f_noisy,alpha=0.3,nbins=50)
limx = [-2.5,0]
limy = [-10,1000]
rap  = 500
labx = "CV values"
laby = "Number of iterations"
p = plot(layout=grid(3,3),legend = false,size=(600,500),dpi=1000,leftmargins=0.25cm,bottommargins=0.1cm,xticks=-2.5:0.5:0,xminorticks=0.25)



p = histogram!(p[1],b_noisy,xlims=limx,nbins=20,aspect_ratio=1/rap,ylims=limy,alpha=0.2,fillcolor=:peru,title=L"snr=0.30E+02",xaxis=:false,xticks=:false)
p = annotate!(p[1],[limx[2]-(limx[2]-limx[1])/6],[limy[end]-limy[end]/20],text(L"\sigma = "*"$(trunc(sqrt(std(b_noisy[:])); sigdigits=2))"*"\n"*L"\mu = "*"$(trunc(mean(b_noisy[:]); sigdigits=3))"*"\n"*L"\gamma = "*"$(trunc(moment(b_noisy[:],3); sigdigits=2))",:peru,8))
p = annotate!(p[1],[limx[1]+(limx[2]-limx[1])/5],[limy[end]-limy[end]/20],text(L"\sigma = "*"$(trunc(sqrt(std(cubebase[:])); sigdigits=2))"*"\n"*L"\mu = "*"$(trunc(mean(cubebase[:]); sigdigits=3))"*"\n"*L"\gamma = "*"$(trunc(moment(cubebase[:],3); sigdigits=2))",:green,8))


p = histogram!(p[2],d_noisy,xlims=limx,nbins=30,aspect_ratio=1/rap,ylims=limy,alpha=0.2,fillcolor=:peru,title=L"snr=0.12E+02",xaxis=:false,xticks=:false,yaxis=:false,yticks=:false)
p = annotate!(p[2],[limx[2]-(limx[2]-limx[1])/6],[limy[end]-limy[end]/20],text(L"\sigma = "*"$(trunc(sqrt(std(d_noisy[:])); sigdigits=2))"*"\n"*L"\mu = "*"$(trunc(mean(d_noisy[:]); sigdigits=3))"*"\n"*L"\gamma = "*"$(trunc(moment(d_noisy[:],3); sigdigits=2))",:peru,8))
p = annotate!(p[2],[limx[1]+(limx[2]-limx[1])/5],[limy[end]-limy[end]/20],text(L"\sigma = "*"$(trunc(sqrt(std(cubebase[:])); sigdigits=2))"*"\n"*L"\mu = "*"$(trunc(mean(cubebase[:]); sigdigits=3))"*"\n"*L"\gamma = "*"$(trunc(moment(cubebase[:],3); sigdigits=2))",:green,8))

p = histogram!(p[3],f_noisy,xlims=limx,nbins=30,aspect_ratio=1/rap,title=L"snr=0.78E+01",ylims=limy,alpha=0.2,fillcolor=:peru,xaxis=:false,xticks=:false,yaxis=:false,yticks=:false)
p = annotate!(p[3],[limx[2]-(limx[2]-limx[1])/6],[limy[end]-limy[end]/20],text(L"\sigma = "*"$(trunc(sqrt(std(f_noisy[:])); sigdigits=2))"*"\n"*L"\mu = "*"$(trunc(mean(f_noisy[:]); sigdigits=3))"*"\n"*L"\gamma = "*"$(trunc(moment(f_noisy[:],3); sigdigits=2))",:peru,8))
p = annotate!(p[3],[limx[1]+(limx[2]-limx[1])/5],[limy[end]-limy[end]/20],text(L"\sigma = "*"$(trunc(sqrt(std(cubebase[:])); sigdigits=2))"*"\n"*L"\mu = "*"$(trunc(mean(cubebase[:]); sigdigits=3))"*"\n"*L"\gamma = "*"$(trunc(moment(cubebase[:],3); sigdigits=2))",:green,8))

for ix=1:3
    global p
   p = histogram!(p[ix],cubebase,fillcolor=:green,xlims=limx,ylims=limy,nbins=20,alpha=0.25)
end

p = histogram!(p[4],b,xlims=limx,nbins=20,aspect_ratio=1/rap,ylims=limy,alpha=0.3,fillcolor=:blue)
p = histogram!(p[4],b_noisy,xlims=limx,nbins=20,aspect_ratio=1/rap,ylims=limy,alpha=0.2,fillcolor=:peru,xaxis=:false,xticks=:false)
p = annotate!(p[4],[limx[1]+(limx[2]-limx[1])/5],[limy[end]-limy[end]/20],text(L"\sigma = "*"$(trunc(sqrt(std(b[:])); sigdigits=2))"*"\n"*L"\mu = "*"$(trunc(mean(b[:]); sigdigits=3))"*"\n"*L"\gamma = "*"$(trunc(moment(b[:],3); sigdigits=2))",:blue,8))
p = annotate!(p[4],[limx[2]-(limx[2]-limx[1])/6],[limy[end]-limy[end]/20],text(L"\sigma = "*"$(trunc(sqrt(std(b_noisy[:])); sigdigits=2))"*"\n"*L"\mu = "*"$(trunc(mean(b_noisy[:]); sigdigits=3))"*"\n"*L"\gamma = "*"$(trunc(moment(b_noisy[:],3); sigdigits=2))",:peru,8))

p = histogram!(p[5],d,xlims=limx,nbins=20,aspect_ratio=1/rap,ylims=limy,alpha=0.3,fillcolor=:blue)
p = histogram!(p[5],d_noisy,xlims=limx,nbins=30,aspect_ratio=1/rap,ylims=limy,alpha=0.2,fillecolor=:peru,xaxis=:false,xticks=:false,yaxis=:false,yticks=:false)
p = annotate!(p[5],[limx[1]+(limx[2]-limx[1])/5],[limy[end]-limy[end]/20],text(L"\sigma = "*"$(trunc(sqrt(std(d[:])); sigdigits=2))"*"\n"*L"\mu = "*"$(trunc(mean(d[:]); sigdigits=3))"*"\n"*L"\gamma = "*"$(trunc(moment(d[:],3); sigdigits=2))",:blue,8))
p = annotate!(p[5],[limx[2]-(limx[2]-limx[1])/6],[limy[end]-limy[end]/20],text(L"\sigma = "*"$(trunc(sqrt(std(d_noisy[:])); sigdigits=2))"*"\n"*L"\mu = "*"$(trunc(mean(d_noisy[:]); sigdigits=3))"*"\n"*L"\gamma = "*"$(trunc(moment(d_noisy[:],3); sigdigits=2))",:peru,8))

p = histogram!(p[6],f,xlims=limx,nbins=20,aspect_ratio=1/rap,ylims=limy,alpha=0.3,fillcolor=:blue)
p = histogram!(p[6],f_noisy,xlims=limx,nbins=30,aspect_ratio=1/rap,ylims=limy,alpha=0.2,fillecolor=:peru,xaxis=:false,xticks=:false,yaxis=:false,yticks=:false)
p = annotate!(p[6],[limx[1]+(limx[2]-limx[1])/5],[limy[end]-limy[end]/8],text(L"\sigma = "*"$(trunc(sqrt(std(f[:]));sigdigits=2))"*"\n"*L"\mu = "*"$(trunc(mean(f[:]); sigdigits=3))"*"\n"*L"\gamma = "*"$(trunc(moment(f[:],3); sigdigits=2))",:blue,8))
p = annotate!(p[6],[limx[2]-(limx[2]-limx[1])/6],[limy[end]-limy[end]/20],text(L"\sigma = "*"$(trunc(sqrt(std(f_noisy[:])); sigdigits=2))"*"\n"*L"\mu = "*"$(trunc(mean(f_noisy[:]); sigdigits=3))"*"\n"*L"\gamma = "*"$(trunc(moment(f_noisy[:],3); sigdigits=2))",:peru,8))

p = histogram!(p[7],b,xlims=limx,nbins=20,aspect_ratio=1/rap,ylims=limy,alpha=0.3,fillcolor=:blue)
p = annotate!(p[7],[limx[1]+(limx[2]-limx[1])/5],[limy[end]-limy[end]/20],text(L"\sigma = "*"$(trunc(sqrt(std(b[:])); sigdigits=2))"*"\n"*L"\mu = "*"$(trunc(mean(b[:]); sigdigits=3))"*"\n"*L"\gamma = "*"$(trunc(moment(b[:],3); sigdigits=2))",:blue,8))
p = annotate!(p[7],[limx[2]-(limx[2]-limx[1])/6],[limy[end]-limy[end]/20],text(L"\sigma = "*"$(trunc(sqrt(std(cubebase[:])); sigdigits=2))"*"\n"*L"\mu = "*"$(trunc(mean(cubebase[:]); sigdigits=3))"*"\n"*L"\gamma = "*"$(trunc(moment(cubebase[:],3); sigdigits=2))",:green,8))

p = histogram!(p[8],d,xlims=limx,nbins=20,aspect_ratio=1/rap,ylims=limy,alpha=0.3,fillcolor=:blue,yaxis=:false,yticks=:false)
p = annotate!(p[8],[limx[1]+(limx[2]-limx[1])/5],[limy[end]-limy[end]/20],text(L"\sigma = "*"$(trunc(sqrt(std(d[:])); sigdigits=2))"*"\n"*L"\mu = "*"$(trunc(mean(d[:]); sigdigits=3))"*"\n"*L"\gamma = "*"$(trunc(moment(d[:],3); sigdigits=2))",:blue,8))
p = annotate!(p[8],[limx[2]-(limx[2]-limx[1])/6],[limy[end]-limy[end]/20],text(L"\sigma = "*"$(trunc(sqrt(std(cubebase[:])); sigdigits=2))"*"\n"*L"\mu = "*"$(trunc(mean(cubebase[:]); sigdigits=3))"*"\n"*L"\gamma = "*"$(trunc(moment(cubebase[:],3); sigdigits=2))",:green,8))

p = histogram!(p[9],f,xlims=limx,nbins=20,aspect_ratio=1/rap,ylims=limy,alpha=0.3,fillcolor=:blue,yaxis=:false,yticks=:false)
p = annotate!(p[9],[limx[1]+(limx[2]-limx[1])/5],[limy[end]-limy[end]/20],text(L"\sigma = "*"$(trunc(sqrt(std(f[:])); sigdigits=2))"*"\n"*L"\mu = "*"$(trunc(mean(f[:]); sigdigits=3))"*"\n"*L"\gamma = "*"$(trunc(moment(f[:],3); sigdigits=2))",:blue,8))
p = annotate!(p[9],[limx[2]-(limx[2]-limx[1])/6],[limy[end]-limy[end]/20],text(L"\sigma = "*"$(trunc(sqrt(std(cubebase[:])); sigdigits=2))"*"\n"*L"\mu = "*"$(trunc(mean(cubebase[:]); sigdigits=3))"*"\n"*L"\gamma = "*"$(trunc(moment(cubebase[:],3); sigdigits=2))",:green,8))
for ix=7:9
    global p
   p = histogram!(p[ix],cubebase,fillcolor=:green,xlims=limx,ylims=limy,nbins=20,alpha=0.25)
end
p = histogram!(p[4],ylabel=laby)

for ix=7:9
    global p
    p = histogram!(p[ix],xlabel=labx)
end
display(p)
#= WORKING WELL
limx = [-2.2,-1.2]
limy = [-10,4000]
rap  = 4010
p = plot(layout=grid(2,4),legend = false,link=:both,size=(900,450),dpi=1000,leftmargins=0.25cm,bottommargins=0.25cm,topmargin=0.01cm)
p = plot!(p[1],showaxis=:hide,foreground_color_text=:white,grid=:false)
p = histogram!(p[2],a[:],fillcolor=:white,xlims=limx,nbins=10,aspect_ratio=1/rap,ylims=limy,xlabel=L"\sigma=0.81E+09")
p = histogram!(p[3],b[:],fillcolor=:white,xlims=limx,nbins=14,aspect_ratio=1/rap,ylims=limy,xlabel=L"\sigma=0.30E+02")
p = histogram!(p[4],c[:],fillcolor=:white,xlims=limx,nbins=14,aspect_ratio=1/rap,ylims=limy,xlabel=L"\sigma=0.16E+02")
p = histogram!(p[5],cubebase[:],orientation=:horizontal, fillcolor=:white,ylims=limx,nbins=10,aspect_ratio=rap+600,xlims=limy,xflip=:true,ylabel="Expected Values")
p = plot!(p[6],cubebase,a,seriestype=:scatter,markersize=0.5,markershape=:+,alpha=0.2,xlims=limx,ylims=limx,color=:black,xlabel="CV values (a.u.)")
p = plot!(p[7],cubebase,b,seriestype=:scatter,markersize=0.5,markershape=:+,alpha=0.2,xlims=limx,ylims=limx,color=:black,xlabel="CV values (a.u.)")
p = plot!(p[8],cubebase,c,seriestype=:scatter,markersize=0.5,markershape=:+,alpha=0.2,xlims=limx,ylims=limx,color=:black,xlabel="CV values (a.u.)")
for ix=6:8
    global p
    p = plot!(p[ix],limx,limx,seriestype=:line,color=:blue,alpha=0.2)
end
display(p)


savefig("$(PATHTOSAVE)/Plots/3first_histo_reconstructed.png")

p = plot(layout=grid(2,4),legend = false,link=:both,size=(900,450),dpi=1000,leftmargins=0.25cm,bottommargins=0.25cm,topmargin=0.01cm) #link=:y,
p = plot!(p[1],showaxis=:hide,foreground_color_text=:white,grid=:false)
p = histogram!(p[2],d[:],fillcolor=:white,xlims=limx,nbins=14,aspect_ratio=1/rap,ylims=limy,xlabel=L"\sigma = 0.12E+02")
p = histogram!(p[3],e[:],fillcolor=:white,xlims=limx,nbins=14,aspect_ratio=1/rap,ylims=limy,xlabel=L"\sigma = 0.92E+01")
p = histogram!(p[4],f[:],fillcolor=:white,xlims=limx,nbins=20,aspect_ratio=1/rap,ylims=limy,xlabel=L"\sigma = 0.78E+01")
p = histogram!(p[5],cubebase[:],orientation=:horizontal, fillcolor=:white,ylims=limx,nbins=10,aspect_ratio=rap+600,xlims=limy,xflip=:true,ylabel="Expected Values")
p = plot!(p[7],cubebase,e,seriestype=:scatter,markersize=0.5,markershape=:+,alpha=0.2,xlims=limx,ylims=limx,color=:black,xlabel="CV values (a.u.)")
p = plot!(p[6],cubebase,d,seriestype=:scatter,markersize=0.5,markershape=:+,alpha=0.2,xlims=limx,ylims=limx,color=:black,xlabel="CV values (a.u.)")
p = plot!(p[8],cubebase,f,seriestype=:scatter,markersize=0.5,markershape=:+,alpha=0.2,xlims=limx,ylims=limx,color=:black,xlabel="CV values (a.u.)")
for ix=6:8
    global p
    p = plot!(p[ix],limx,limx,seriestype=:line,color=:blue,alpha=0.2)
end
display(p)
savefig("$(PATHTOSAVE)/Plots/3last_histo_reconstructed.png")


=#
#=
data = (;cubebase,a,b,c,d,e,f)
corner(data,["Real values",asnr,bsnr,csnr,dsnr,esnr,fsnr],plotcontours=false,scatter_kwargs=(;xlims=[-6.,3],ylims=[-6,3],color=:red),contour_kwargs=(;xlims=[-6,3],ylims=[-6,3]),hist_kwargs=(;xlims=[-10,2],ylims=[0,300]),hist2d_kwargs=(;color=:inferno),title="Corner plot of CV calculated on simple ppv cube. Names are mean snr values of the cube.",titlefontsize=10)#,contour_kwargs=(xlims=[-10.,10.],ylims=[-10.,10.]))#,hist2d_kwargs=(;xlims=[-10,10]),hist_kwargs=(;xlims=[-10,10]))#,contour_kwargs=(xlims=[-8,8],ylims=[-8,8]))


savefig("$(PATHTOSAVE)/Plots/cornerplotlarge.pdf")
=#
=#