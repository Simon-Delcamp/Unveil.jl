###################################################################
# Compare multiple CVI cubes.
###################################################################

include("../src/Data_preparation.jl") # Read and write fits
include("../src/Functionforcvi.jl")   # Calculations of CVI
include("../src/Graphic.jl")
#include("../../Julia_prog/OutPackages/PairPlots.jl/src/PairPlots.jl")

using StatsPlots

using .Data_preparation
using .Functionforcvi
using .Graphic
#using .PairPlots

cvi0,DATADIMENSION,HEADER   = Data_preparation.read_fits_ppv("/home/delcamps/Data/Pipe_nebula/TestPCA/EvolvCVWind/cvirelative_multlag_0.fits","km/s",check=false) 
cvi0  = reshape(cvi0,DATADIMENSION[1]*DATADIMENSION[2],DATADIMENSION[3])  
cvi1   = reshape(Data_preparation.read_fits_ppv("/home/delcamps/Data/Pipe_nebula/TestPCA/EvolvCVWind/cvirelative_multlag_1.fits","km/s",check=false)[1],DATADIMENSION[1]*DATADIMENSION[2],DATADIMENSION[3])  
cvi2   = reshape(Data_preparation.read_fits_ppv("/home/delcamps/Data/Pipe_nebula/TestPCA/EvolvCVWind/cvirelative_multlag_2.fits","km/s",check=false)[1],DATADIMENSION[1]*DATADIMENSION[2],DATADIMENSION[3])   
cvi3   = reshape(Data_preparation.read_fits_ppv("/home/delcamps/Data/Pipe_nebula/TestPCA/EvolvCVWind/cvirelative_multlag_3.fits","km/s",check=false)[1],DATADIMENSION[1]*DATADIMENSION[2],DATADIMENSION[3])   
cvi4   = reshape(Data_preparation.read_fits_ppv("/home/delcamps/Data/Pipe_nebula/TestPCA/EvolvCVWind/cvirelative_multlag_4.fits","km/s",check=false)[1],DATADIMENSION[1]*DATADIMENSION[2],DATADIMENSION[3])   
cvi5   = reshape(Data_preparation.read_fits_ppv("/home/delcamps/Data/Pipe_nebula/TestPCA/EvolvCVWind/cvirelative_multlag_5.fits","km/s",check=false)[1],DATADIMENSION[1]*DATADIMENSION[2],DATADIMENSION[3])   
cvi6   = reshape(Data_preparation.read_fits_ppv("/home/delcamps/Data/Pipe_nebula/TestPCA/EvolvCVWind/cvirelative_multlag_6.fits","km/s",check=false)[1],DATADIMENSION[1]*DATADIMENSION[2],DATADIMENSION[3])   
cvi7   = reshape(Data_preparation.read_fits_ppv("/home/delcamps/Data/Pipe_nebula/TestPCA/EvolvCVWind/cvirelative_multlag_7.fits","km/s",check=false)[1],DATADIMENSION[1]*DATADIMENSION[2],DATADIMENSION[3])   
cvi8   = reshape(Data_preparation.read_fits_ppv("/home/delcamps/Data/Pipe_nebula/TestPCA/EvolvCVWind/cvirelative_multlag_8.fits","km/s",check=false)[1],DATADIMENSION[1]*DATADIMENSION[2],DATADIMENSION[3])   
cvi9   = reshape(Data_preparation.read_fits_ppv("/home/delcamps/Data/Pipe_nebula/TestPCA/EvolvCVWind/cvirelative_multlag_9.fits","km/s",check=false)[1],DATADIMENSION[1]*DATADIMENSION[2],DATADIMENSION[3])   
cvi10  = reshape(Data_preparation.read_fits_ppv("/home/delcamps/Data/Pipe_nebula/TestPCA/EvolvCVWind/cvirelative_multlag_10.fits","km/s",check=false)[1],DATADIMENSION[1]*DATADIMENSION[2],DATADIMENSION[3])  
cvi11  = reshape(Data_preparation.read_fits_ppv("/home/delcamps/Data/Pipe_nebula/TestPCA/EvolvCVWind/cvirelative_multlag_11.fits","km/s",check=false)[1],DATADIMENSION[1]*DATADIMENSION[2],DATADIMENSION[3])  
cvi12  = reshape(Data_preparation.read_fits_ppv("/home/delcamps/Data/Pipe_nebula/TestPCA/EvolvCVWind/cvirelative_multlag_12.fits","km/s",check=false)[1],DATADIMENSION[1]*DATADIMENSION[2],DATADIMENSION[3])  
cvi13  = reshape(Data_preparation.read_fits_ppv("/home/delcamps/Data/Pipe_nebula/TestPCA/EvolvCVWind/cvirelative_multlag_13.fits","km/s",check=false)[1],DATADIMENSION[1]*DATADIMENSION[2],DATADIMENSION[3])  
cvi14  = reshape(Data_preparation.read_fits_ppv("/home/delcamps/Data/Pipe_nebula/TestPCA/EvolvCVWind/cvirelative_multlag_14.fits","km/s",check=false)[1],DATADIMENSION[1]*DATADIMENSION[2],DATADIMENSION[3])  

#cvi0 = Data_preparation.replace_nantoblank(cvi0,-1000)
#cvi1 = Data_preparation.replace_nantoblank(cvi1,-1000)
#cvi2 = Data_preparation.replace_nantoblank(cvi2,-1000)


M = Array{Float64}(undef,size(cvi0)[1],15)
M[:,1] = cvi0[:,1]
M[:,2] = cvi1[:,1]
M[:,3] = cvi2[:,1]
M[:,4] = cvi3[:,1]
M[:,5] = cvi4[:,1]
M[:,6] = cvi5[:,1]
M[:,7] = cvi6[:,1]
M[:,8] = cvi7[:,1]
M[:,9] = cvi8[:,1]
M[:,10] = cvi9[:,1]
M[:,11] = cvi10[:,1]
M[:,12] = cvi11[:,1]
M[:,13] = cvi12[:,1]
M[:,14] = cvi13[:,1]
M[:,15] = cvi14[:,1]
cornerplot(M[:,1:6],compact=true)

#table = (; cvi1 = cvi1[:,1], cvi2 = cvi2[:,1])

#PairPlots.pairplot(table)


