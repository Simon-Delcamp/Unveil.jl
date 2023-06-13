###################################################################
# Function used for PCA computations processes. 
# Called these functions by calling (or writting on a script ):
#       >include("../src/Functionforpca.jl")
#       >using .Functionforpca
#       >output = Functionforpca.NameOfTheFunction(input)
###################################################################

module Functionforpca

include("Functionforcvi.jl") #Calculations of CVI
include("Dataprep.jl") #Read and write fits

using .Dataprep
using .Functionforcvi

export cv_convergence
export convergence_cvdistribution!
export convergence_moments_byintegration
export convergence_pcproj!
export intensity_moments_specific_channels
export intensity_moments_specific_channels_withPCA
export intensity_moments_specific_channels_withPCASVDRAND
export multiplepca
export multiple_moment
export multiple_moment!
export pca
export pca_nomorecalc
export variance_evolution

using MultivariateStats, Plots, Statistics,StatsBase
using Mmap



"""
        cv_convergence(pcmax,M,Yt,velocity_vector,datadim)

Return a cube reconstructed with different number of PC, its first and its second velocity moment order. The cube considered is the one processed by PCA, that gave M and Yt as a result.
The arrays in output are in 3D, where the third dimensions gives results for different number of PC used. Datadim should be the dimension of the 3D initial PPV cube. M is a PCA type, produced by the function pca. Yt is a vector or a matrix where each column gives the principal components of an observation.
"""
function cv_convergence(pcmax,M,Yt,velocity_vector,datadim)
        data_multiplepc     = multiple_pca(M,Yt,pcmax,datadim)
        data_multiplemoment_one,data_multiplemoment_two = Functionforcvi.multiple_moment(data_multiplepc,velocity_vector)
        return(data_multiplepc,data_multiplemoment_one,data_multiplemoment_two)
end



function convergence_cvdistribution!(M,Yt,HIGHESTPC,DATADIM,VELOCITYVECTOR,BLANK,cvcube)
        for pcx=1:HIGHESTPC
                RECON = pca_nomorecalc(M,Yt,1,pcx,DATADIM)
                cvcube[:,pcx] = Functionforcvi.moment_one_field(RECON,VELOCITYVECTOR,BLANK)
        end#for pcx
        return(cvcube)
end #function convergence_cvdistribution



"""
        convergence_moments_byintegration(data2D,integration_range,pcnumber,velocity_increment,M,Yt,datadim,velocity_vector,IM)

Return the I-th moment order of the integrated by sections spectras differences, from arrays reconstructed with N PC and arrays without PCA. Also return the new velocity vector where each value is centered on the integrated sections. The function integrate by section (the integration_range) spectra from reconstructed and non-reconstructed data, then compute the given moment order of their differences. Data should be in 2 dimensions (PV cube), and datadim should refers to the dimension of the 3D data cube (PPV). Pcnumber should be a vector with all the number of PC you want (example Vector(1:10)).
"""
function convergence_moments_byintegration(data2D_path,integration_range,pcnumber,M,Yt_path,datadim,velocity_vector,IM)
        velocity_increment       = abs.(velocity_vector[2]-velocity_vector[1])
        numberofsection          = trunc(abs(datadim[2]*velocity_increment)/integration_range) |> Int
        sectionsize_indices      = trunc(datadim[2]/numberofsection)|> Int
        (sectionsize_indices==0) && (sectionsize_indices=1)
        (numberofsection==0)     && (numberofsection=1)
        sectionsize_indices%2==0 && (velocity_newvector = velocity_vector[1:trunc(Int,sectionsize_indices):end])
        velocity_newvector = velocity_vector[1:trunc(Int,(sectionsize_indices)):end]
        size(velocity_newvector)[1]<=numberofsection && push!(velocity_newvector,velocity_vector[end])
        size(velocity_newvector)[1]>=numberofsection && deleteat!(velocity_newvector,size(velocity_newvector)[1])
        mom               =  zeros(Float64,numberofsection,size(pcnumber)[1])
        Threads.@threads  for pc=pcnumber[1]:pcnumber[end]
                s = open("$(Yt_path)")
                Yt = Mmap.mmap(s,Matrix{Float64},(pcnumber[end],datadim[2]))
                close(s)
                data_reconstructed = pca_nomorecalc(M,Yt,1,pc,datadim)
                Yt = 0.0
                data_reconstructed = Data_preparation.blank_inf(data_reconstructed,-80,0)
                data_reconstructed = Data_preparation.blank_sup(data_reconstructed,52,0)
                GC.gc()
                #Threads.@threads for ix=1:numberofsection
                for ix=1:numberofsection
                        s = open("$(data2D_path)")
                        data2D = Mmap.mmap(s,Matrix{Float64},(datadim[1],datadim[2]))
                        close(s)
                        if ((ix+1)*sectionsize_indices)>datadim[2]
                                summedrecon = sum(data_reconstructed[:,ix*sectionsize_indices:end],dims=2).*abs(velocity_increment)
                                summedraw = sum(data2D[:,ix*sectionsize_indices:end],dims=2).*abs(velocity_increment)
                                data2D = 0.0
                                if IM==1 
                                        mom[ix,pc-pcnumber[1]+1] = moment(abs.(summedraw.-summedrecon),1,0)
                                else
                                        mom[ix,pc-pcnumber[1]+1] = moment(abs.(summedraw.-summedrecon),IM)
                                end
                        else
                                summedrecon = sum(data_reconstructed[:,ix*sectionsize_indices:(ix+1)*sectionsize_indices],dims=2).*abs(velocity_increment)
                                summedraw = sum(data2D[:,ix*sectionsize_indices:(ix+1)*sectionsize_indices],dims=2).*abs(velocity_increment)
                                data2D = 0.0
                                if IM==1
                                        mom[ix,pc-pcnumber[1]+1] = moment(abs.(summedraw.-summedrecon),1,0)
                                else
                                        mom[ix,pc-pcnumber[1]+1] = moment(abs.(summedraw.-summedrecon),IM)
                                end
                        end
                end
        end
        return(mom,velocity_newvector)
end
#
#function convergence_moments_byintegration(data2D_path,integration_range,pcnumber,M,Yt_path,datadim,velocity_vector,IM)
#        velocity_increment       = abs.(velocity_vector[2]-velocity_vector[1])
#        numberofsection          = trunc(abs(datadim[2]*velocity_increment)/integration_range) |> Int
#        sectionsize_indices      = trunc(datadim[2]/numberofsection)|> Int
#        (sectionsize_indices==0) && (sectionsize_indices=1)
#        (numberofsection==0)     && (numberofsection=1)
#        sectionsize_indices%2==0 && (velocity_newvector = velocity_vector[1:trunc(Int,sectionsize_indices):end])
#        velocity_newvector = velocity_vector[1:trunc(Int,(sectionsize_indices)):end]
#        size(velocity_newvector)[1]<=numberofsection && push!(velocity_newvector,velocity_vector[end])
#        size(velocity_newvector)[1]>=numberofsection && deleteat!(velocity_newvector,size(velocity_newvector)[1])
#        mom               =  zeros(Float64,numberofsection,size(pcnumber)[1])
#        @inbounds @views  for pc=pcnumber[1]:pcnumber[end]
#                s = open("$(Yt_path)")
#                Yt = Mmap.mmap(s,Matrix{Float64},(pcnumber[end],datadim[2]))
#                close(s)
#                data_reconstructed = pca_nomorecalc(M,Yt,1,pc,datadim)
#                Yt = 0.0
#                data_reconstructed = Data_preparation.blank_inf(data_reconstructed,-80,0)
#                data_reconstructed = Data_preparation.blank_sup(data_reconstructed,52,0)
#                GC.gc()
#                @inbounds @views for ix=1:numberofsection
#                        s = open("$(data2D_path)")
#                        data2D = Mmap.mmap(s,Matrix{Float64},(datadim[1],datadim[2]))
#                        close(s)
#                        if ((ix+1)*sectionsize_indices)>datadim[2]
#                                summedrecon = sum(data_reconstructed[:,ix*sectionsize_indices:end],dims=2).*abs(velocity_increment)
#                                summedraw = sum(data2D[:,ix*sectionsize_indices:end],dims=2).*abs(velocity_increment)
#                                data2D = 0.0
#                                if IM==1 
#                                        mom[ix,pc-pcnumber[1]+1] = moment(abs.(summedraw.-summedrecon),1,0)
#                                else
#                                        mom[ix,pc-pcnumber[1]+1] = moment(abs.(summedraw.-summedrecon),IM)
#                                end
#                        else
#                                summedrecon = sum(data_reconstructed[:,ix*sectionsize_indices:(ix+1)*sectionsize_indices],dims=2).*abs(velocity_increment)
#                                summedraw = sum(data2D[:,ix*sectionsize_indices:(ix+1)*sectionsize_indices],dims=2).*abs(velocity_increment)
#                                data2D = 0.0
#                                if IM==1
#                                        mom[ix,pc-pcnumber[1]+1] = moment(abs.(summedraw.-summedrecon),1,0)
#                                else
#                                        mom[ix,pc-pcnumber[1]+1] = moment(abs.(summedraw.-summedrecon),IM)
#                                end
#                        end
#                end
#        end
#        return(mom,velocity_newvector)
#end
#





function convergence_pcproj!(M,cube,proj)
        matpro = projection(M)
        for px=1:size(proj)[2]
                proj[:,px] = cube.*matpro[:,px]
        end #for px
        return(proj)
end


"""
        intensity_moments_specific_channels(arr_path,channels,datadim)

Will compute the first, second, third and fourth, moments order of the mean of given channels in a data cube. The cube is not given as input but only the path where it is saved as a 'bin' temporary file (saving RAM). Datadim refer to the dimension of the 3D data cube (the original data cube).
"""
function intensity_moments_specific_channels(arr_path,channels,datadim)
        s = open("$(arr_path)")
        arr = Mmap.mmap(s,Matrix{Float64},(datadim))
        close(s)
        mean_one = mean.(arr[:,channels[1]])
        mean_second = mean.(arr[:,channels[2]])
        arr = 0.0
        mom1 = (moment(mean_one,1,0)+moment(mean_second,1,0))/2
        mom2 = (moment(mean_one,2)+moment(mean_second,2))/2
        mom3 = (moment(mean_one,3)+moment(mean_second,3))/2
        mom4 = (moment(mean_one,4)+moment(mean_second,4))/2
        allmoment = []
        push!(allmoment,mom1,mom2,mom3,mom4)
        return(allmoment)
end



"""
        intensity_moments_specific_channels_withPCA(M,Yt_path::String,pcnumber,canals,datadim)

Will compute the first, second, third and fourth moments order of the mean of specific
canals in the data reconstructed with different number of PC to check if they are converging.
Datadim refer to the dimension of the 3D data cube (the original data cube). M is a PCA type, produced by the function pca. Yt_path is the path where Yt is saved as a 'bin' temporary fil. Yt is a vector or a matrix where each column gives the principal components of an observation. M and Yt are obtained with the function pca.
Moment order can be any Int from 1 to 4.
Return one 1D array.
"""
function intensity_moments_specific_channels_withPCA(M,Yt_path::String,pcnumber,canals,datadim)
        mom1_allpc = zeros(Float64,pcnumber[end],1)
        mom2_allpc = similar(mom1_allpc)
        mom3_allpc = similar(mom1_allpc)
        mom4_allpc = similar(mom1_allpc)
        moment_pca = []
        for nbpc = 1:pcnumber[end]
                s = open("$(Yt_path)")
                Yt = Mmap.mmap(s,Matrix{Float64},(pcnumber,datadim[2]))
                close(s)
                data_reconstructed = pca_nomorecalc(M,Yt,1,nbpc,datadim)
                Yt = 0.0
                # Calculate the mean of the moment order N
                mean_data_one = mean.(data_reconstructed[:,canals[1]])
                mean_data_second = mean.(data_reconstructed[:,canals[2]])
                data_reconstructed = 0.0
                mom1_allpc[nbpc] = (moment(mean_data_one,1,0)+moment(mean_data_second,1,0))/2
                mom2_allpc[nbpc] = (moment(mean_data_one,2)+moment(mean_data_second,2))/2
                mom3_allpc[nbpc] = (moment(mean_data_one,3)+moment(mean_data_second,3))/2
                mom4_allpc[nbpc] = (moment(mean_data_one,4)+moment(mean_data_second,4))/2

        end
        push!(moment_pca,mom1_allpc)
        push!(moment_pca,mom2_allpc)
        push!(moment_pca,mom3_allpc)
        push!(moment_pca,mom4_allpc)   

end


"""
        intensity_moments_specific_channels_withPCASVDRAND(svdobject,pcnumber,channels,p)  

Same as intensity_moments_specific_channels_withPCA but using a Principal Component Analysis based on random SVD.
"""
function intensity_moments_specific_channels_withPCASVDRAND(svdobject,pcnumber,channels,p)
        mom1_allpc = zeros(Float64,pcnumber[end],1)
        mom2_allpc = similar(mom1_allpc)
        mom3_allpc = similar(mom1_allpc)
        mom4_allpc = similar(mom1_allpc)
        moment_pca = []
        @views @inbounds for nbpc = 1:pcnumber[end]
                data_reconstructed = transpose(svdobject.S[1:nbpc].*transpose(svdobject.U[:,1:nbpc]))*svdobject.Vt[1:nbpc,:]
                # Calculate the mean of the moment order N
                mean_data_one = mean.(data_reconstructed[:,channels[1]])
                mean_data_second = mean.(data_reconstructed[:,channels[2]])
                data_reconstructed = 0.0
                mom1_allpc[nbpc] = (moment(mean_data_one,1,0)+moment(mean_data_second,1,0))/2
                mom2_allpc[nbpc] = (moment(mean_data_one,2)+moment(mean_data_second,2))/2
                mom3_allpc[nbpc] = (moment(mean_data_one,3)+moment(mean_data_second,3))/2
                mom4_allpc[nbpc] = (moment(mean_data_one,4)+moment(mean_data_second,4))/2

        end
        push!(moment_pca,mom1_allpc)
        push!(moment_pca,mom2_allpc)
        push!(moment_pca,mom3_allpc)
        push!(moment_pca,mom4_allpc)   
        @views @inbounds for ix=2:size(moment_pca)[1]+1
                p = plot!(p[ix],moment_pca[ix-1],xlabel="Number of PC",title="$(ix-1) moment order",seriestype=:scatter,markershape=:+,markersize=2.5,label="",minorgrid=true)
        end
end



"""
        multiplepca(M,Yt,pcnumber,datadim)

Do multiple PCA with different number of PC on the same data. Will return a 3D array, with pixel position in the first dimension, velocity in the second dimension, and the data reconstructed with a specific number of pc in the third dimension. Pcnumber should be a vector with all the number of PC you want (example Vector(1:10)).  M is a PCA type, produced by the function pca. Yt is a vector or a matrix where each column gives the principal components of an observation.
"""
function multiplepca(M,Yt,pcnumber,datadim)
        typeof(datadim)==Tuple{Int64,Int64,Int64} && (data_multiplepc  = zeros(Float64,datadim[1]*datadim[2],datadim[3],trunc(Int,size(pcnumber)[1])))
        typeof(datadim)==Tuple{Int64,Int64}       && (data_multiplepc  = zeros(Float64,datadim[1],datadim[2],trunc(Int,size(pcnumber)[1])))
        for (index, pcvalue) in enumerate(pcnumber)
                typeof(datadim)==Tuple{Int64,Int64,Int64} && (data_multiplepc[:,:,index] = reshape(pca_nomorecalc(M,Yt,1,pcvalue,datadim),datadim[1]*datadim[2],datadim[3]))
                typeof(datadim)==Tuple{Int64,Int64}       && (data_multiplepc[:,:,index] = reshape(pca_nomorecalc(M,Yt,1,pcvalue,datadim),datadim[1],datadim[2]))
        end
        return(data_multiplepc)
end



"""
        multiple_moment(velvector,data_dimension,M,Yt,pcmax)

Compute first and second velocity moment order of multiple data cube. The data cubes are reconstructed from PCA with different number of PC. 
M is a PCA type, produced by the function pca. Yt is a vector or a matrix where each column gives the principal components of an observation.
"""
function multiple_moment(velvector,data_dimension,M,Yt,pcmax)
    momentone = zeros(Float64,data_dimension[1],size(pcmax)[1])
    momenttwo = similar(momentone)
    @views @inbounds for ix = 1:pcmax[end]-pcmax[1]
        arr_reconstructed = pca_nomorecalc(M,Yt,1,pcmax[ix],data_dimension)
        @views @inbounds for pix = 1:data_dimension[1]
            momentone[pix,ix] = moment_one(arr_reconstructed[pix,:],velvector)
            momenttwo[pix,ix] = moment_two(arr_reconstructed[pix,:],velvector)
        end
    end
    return(momentone,momenttwo)
end



"""
        multiple_moment!

Same as multiple_moment but on pre-allocated data.
        """
function multiple_moment!(velvector,data_dimension,M,Yt,pcmax,momentone,momenttwo,arr_reconstructed)
        indexing = Vector{Int}(1:size(arr_reconstructed)[3])
        @inbounds @views for ix in eachindex(indexing)
                arr_reconstructed = pca_nomorecalc(M,Yt,1,pcmax[ix],data_dimension)
                @inbounds @views for pix = 1:data_dimension[1]
                        momentone[pix,ix] = moment(velvector,1,aweights(arr_reconstructed[pix,:,ix]),0)
                        momenttwo[pix,ix] = moment(velvector,2,aweights(arr_reconstructed[pix,:,ix]))
                end
        end
        return(momentone,momenttwo)
end


"""
        pca(data2D,pc::Integer,path::String ; percent=1.0)

Produce a Principal Component Analysis on a data, using 'pc' number of principal component. 'percent' gives the wanted percentage of variance reconstructed (default=1). 'path' is the path where results will be saved (in 'bin' temporary files, named as mmap_Yt.bin, mmap_Mproj.bin and mmap_Mmean.bin)
Data are converted in an Array{Float64,2} to do the PCA ( type as Union{Missing} will produce an error). Data should be given in 2D (pixels per velocities).
Will return the path to the projections of the data2D on the PC axis, the mean of these projections, the path to the matrix composed of the principal components in columns and the percentage of explained variance. These results will allow to reconstruct the cube with a given number of PC (using pca_nomorecalc).
"""
function pca(data2D,pc::Integer,path::String ; percent=1.0)
        println("M fit")
        M          = fit(PCA, data2D ; maxoutdim=pc , pratio=percent)
        println("Yt ")
        Yt         = transform(M,data2D)
        data2D     = 0
        s = open("$(path)/mmap_Yt.bin", "w+")
        write(s,Yt)
        close(s)
        Yt = 0
        Yt_path = "$(path)/mmap_Yt.bin"
        varpercent = cumsum(principalvars(M)/tprincipalvar(M))
        # Stock useful matrix to reconstruct the data without using M 
        println("Stock useful matrix to reconstruct the data without using M ")
        s = open("$(path)/mmap_Mproj.bin", "w+")
        write(s,M.proj)
        close(s)
        Mproj_path = "$(path)/mmap_Mproj.bin"
        s = open("$(path)/mmap_Mmean.bin", "w+")
        write(s,M.mean)
        close(s)
        M = 0
        Mmean_path = "$(path)/mmap_Mmean.bin"
        return  Mproj_path,Mmean_path,Yt_path,varpercent
end

"""
        pca(data2D,pc::Integer ; percent=1.0)

Produce a Principal Component Analysis on a data, using 'pc' number of principal component. 'percent' gives the wanted percentage of variance reconstructed (default=1). 
Data are converted in an Array{Float64,2} to do the PCA ( type as Union{Missing} will produce an error). Data should be given in 2D (pixels per velocities).
Will return the pca fit M (a type created by the MultivariateStats package), the matrix composed of the principal components in columns, the percentage of explained variance and the data reconstructed with 'pc' number of principal component. 
"""
function pca(data2D,pc::Integer ; percent=1.0)
        M          = fit(PCA, data2D ; maxoutdim=pc , pratio=percent)
        Yt         = transform(M,data2D)
        data2D     = 0.0
        GC.gc()
        varpercent = cumsum(principalvars(M)/tprincipalvar(M))
        data_reconstructed = reconstruct(M,Yt) # reconstructed data without missing values
        return  M,Yt,varpercent,data_reconstructed
end



"""
        pca_nomorecalc(M,Yt,pcfirst,pclast,Datadim)

Reconstruct data with a certain number of PC, without redo all of the PCA calculation. Necessary to do a first pca calculation to run this one, because needed M (a new type related to PCA constructed by MultivariateStats package) and Yt (matrix with each column a PC). Datadim should be a tuple with two values (PV cube, 2D) or three (PPV cube, 3D).
"""
function pca_nomorecalc(M,Yt,pcfirst,pclast,datadim)
        typeof(datadim)==Tuple{Int64,Int64,Int64} && return(decentralize(M.proj[1:size(M.proj)[1],pcfirst:pclast] * Yt[pcfirst:pclast,1:datadim[3]], M.mean))
        typeof(datadim)==Tuple{Int64,Int64}       && return(decentralize(M.proj[1:size(M.proj)[1],pcfirst:pclast] * Yt[pcfirst:pclast,1:datadim[2]], M.mean))
        error("Not a valid dimension ; need 2 (PV cube) or 3 (PPV cube) values for datadim")
end 



"""
        pca_nomorecalc(Mproj_path::String,Yt_path::String,Mmean_path::String,pcfirst,pclast,datadim)

Same as pca_nomorecalc(M,Yt,pcfirst,pclast,Datadim) but using here the path to the temporary files containing the informations (M and Yt). To be used after pca(data2D,pc::Integer,path::String ; percent=1.0) in order to produce a reconstructed cube.
"""
function pca_nomorecalc(Mproj_path::String,Yt_path::String,Mmean_path::String,pcfirst,pclast,datadim)
        s = open("$(Yt_path)")
        Yt = Mmap.mmap(s,Matrix{Float64},(datadim[2]-1,datadim[2]))
        close(s)

        s = open("$(Mproj_path)")
        Mproj = Mmap.mmap(s,Matrix{Float64},(datadim))
        close(s)

        s = open("$(Mmean_path)")
        Mmean = Mmap.mmap(s,Vector{Float64},(datadim[1]))
        close(s)
        
        if typeof(datadim)==Tuple{Int64,Int64,Int64}
                reconstructed = decentralize(Mproj[1:size(Mproj)[1],pcfirst:pclast] * Yt[pcfirst:pclast,1:datadim[3]],Mmean)
                Yt = 0.0
                Mproj = 0.0
                Mmean = 0.0
                s = open("/tmp/mmap_reconstructed.bin", "w+")
                write(s,reconstructed)
                close(s)
                reconstructed = 0.0
                return("/tmp/mmap_reconstructed.bin")
        elseif typeof(datadim)==Tuple{Int64,Int64}  
                reconstructed = decentralize(Mproj[1:size(Mproj)[1],pcfirst:pclast] * Yt[pcfirst:pclast,1:datadim[2]], Mmean)
                Yt = 0.0
                Mproj = 0.0
                Mmean = 0.0
                s = open("/tmp/mmap_reconstructed.bin", "w+")
                write(s,reconstructed)
                close(s)
                reconstructed = 0.0
                return("/tmp/mmap_reconstructed.bin")
        end
        error("Not a valid dimension ; need 2 (PV cube) or 3 (PPV cube) values for datadim")
end #pca_nomorecalc



"""
        variance_evolution(M,threshold)

Return an Array with the numbers of PC explaining a percentage of the data variance tinier (or equal) than the threshold given in input.
"""
function variance_evolution(M,threshold)
        explained_percentage = principalvars(M)./tprincipalvar(M)
        tot                  = cumsum(explained_percentage)
        pcnumber             = Int64[]
        for ix=1:size(explained_percentage)[1]
                if     (100 .-tot[ix].*100)>=threshold
                           push!(pcnumber,ix)
                elseif (100 .-tot[ix].*100)<threshold
                           (100 .-tot[ix+1].*100)<threshold && return(pcnumber)
                           (100 .-tot[ix+1].*100)>=threshold && push!(pcnumber,ix)
                else
                           (100 .-tot[ix].*100)>=threshold/100 && push!(pcnumber,ix)
                end
        end
        return(pcnumber)
end


# <<<<<<<<<<<<<<<<< TESTING >>>>>>>>>>>>>>>>>>>>>>>>
function fit_pca(m,n,k)
	# matrix to encode
	A = randn(m,k)*randn(k,n)
	loss = QuadLoss()
	r = ZeroReg()
	glrm = GLRM(A,loss,r,r,k)
	X,Y,ch = fit!(glrm)
	println("Convergence history:",ch.objective)
	return A,X,Y,ch
end


end
