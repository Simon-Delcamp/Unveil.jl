module PCA

export pca
export pca_nomorecalc
export proj

using MultivariateStats

"""
        pca(data2D,pc::Integer,path::String ; percent=1.0)

Produce a Principal Component Analysis on a data, using 'pc' number of principal component. 'percent' gives the wanted percentage of variance reconstructed (default=1). 'path' is the path where results will be saved (in 'bin' temporary files, named as mmap_Yt.bin, mmap_Mproj.bin and mmap_Mmean.bin)
Data are converted in an Array{Float64,2} to do the PCA ( type as Union{Missing} will produce an error). Data should be given in 2D (pixels per velocities).
Will return the path to the projections of the data2D on the PC axis, the mean of these projections, the path to the matrix composed of the principal components in columns and the percentage of explained variance. These results will allow to reconstruct the cube with a given number of PC (using pca_nomorecalc).
"""
function pca(data2D,pc::Integer,path::String ; percent=1.0)
        println("M fit")
        M          = fit(MultivariateStats.PCA, data2D ; maxoutdim=pc , pratio=percent)
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
        M          = fit(MultivariateStats.PCA, data2D ; maxoutdim=pc , pratio=percent)
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


function proj(M)
    return(projection(M))
end #proj


end #module