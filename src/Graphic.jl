###################################################################
# Functions for plotting, graphical purposes
# Called these functions by calling (or writting on a script ):
#       >include("../src/Graphic.jl")
#       >using .Graphic
#       >output = Graphic.NameOfTheFunction(input)
###################################################################
module Graphic


include("Data_preparation.jl") #Read and write fits
include("Structure_functions.jl")
include("Functionforpca.jl") #Calculations of PCA

using .Functionforpca

using .Structure_functions
using .Data_preparation


using Plots,MultivariateStats, Statistics, StatsBase, Distributions, LinearAlgebra
using Formatting, LaTeXStrings, KernelDensity, StatsPlots,LsqFit, Colors
using Makie
using GLMakie
using Measures

import StatsPlots.plot
import StatsPlots.plot!
import StatsPlots.text
import StatsPlots.@recipe
import Makie.heatmap!
import StatsPlots.heatmap

export animate_ppvcube
export animate_spectra
export checkwindowopti
export contourmap
export corner_cvimap
export cvi_pdf
export cvi_pdf_norm
export distribcv_multipc
export heatmap_compare_multiplepc
export heatmap_subplot_two
export heatmap_subplot_same
export ifdiff
export meanpc
export mom_conv
export moment_byintegration
export moment_multiplepc
export moment_multiplepc_withraw
export moment_specific_canals
export pixels_averaged_spectrum
export plpc
export ploptiwind
export pratio
export region_spectrum
export region_spectrum_twodata
export sci_not
export struct_function
export plot_pdf
export StcFct
export StcFctWithFit
export StcFctExponent
export StcFctExponentWithError


import Measures:mm
import PyPlot



"""
    animate_ppvcube(cube,data_name)

Plot a heatmap of a 3D cube along third dimension, and save it as a gif file.
"""
function animate_ppvcube(cube,data_name,velocity_vector)
    anim = @animate for ix=1:5:size(velocity_vector)[1]
        heatmap(cube[:,:,ix],clims=(-1,maximum(cube)),aspect_ratio=:equal)
    end
    gif(anim, plotsdir("$(data_name)","ppvcube_from_fbm_truc.gif"),fps=6)
end



"""
    animate_spectra(arr,positionnumber,direction)

Make an animation (.gif) plotting all the spectrum of a line of pixels (dims=1)
or a row of pixels (dims=2). Array has to be 3D.
Used the array, the position of the line or the row, and the direction (line=1 or row=2).
"""
function animate_spectra(arr,positionnumber,direction,frames,field)
    anim = Animation()
    p = plot(arr[:,1],layout=1, ylims=[-1,30],xformatter=:scientific, xlabel="Velocity (km/s)" , ylabel="Emission (Jy/beam)")
    direction==1 && (until=size(arr)[2])
    direction==2 && (until=size(arr)[1])
    if direction==1
        for ix=1:until
            titl = "Column position $(ix)"
            p = plot(arr[positionnumber,ix,:],title=titl)
        end
    elseif direction==2
        for ix=1:until
            titl = "Row position $(ix)"
            p = plot(arr[ix,positionnumber,:],title=titl)
        end
    end
    frame(anim)
    direction==1 && gif(anim, plotsdir("field","spectra_allcolumn_row$(positionnumber).gif"),fps=frames)
    direction==2 && gif(anim, plotsdir("field","spectra_$(positionnumber)columns_allrow.gif"),fps=frames)
end


"""

"""
function checkwindowopti(sourcecube,cubewo,VELOCITYVECTOR,NBROW,NBCOL;limx=(minimum(VELOCITYVECTOR),maximum(VELOCITYVECTOR)))
    l = @layout [grid(NBROW,NBCOL)]
    p = plot(layout=l,legend = false, xformatter=_->"", yformatter=_->"", grid=:false,primary=false,link=:both )#,size=(1000,700),dpi=1000)
    xsize = size(sourcecube)[1]
    posx = []
    for ix = 1:NBROW*NBCOL
        append!(posx,rand(1:xsize))
        p = plot!(p[ix],VELOCITYVECTOR,sourcecube[posx[ix],:],color=:navy,xlim=limx)
        p = plot!(p[ix],VELOCITYVECTOR,cubewo[posx[ix],:],color=:red,title="$(posx[ix])",titlefontsize=5,xlim=limx)
    end
    #p = plot!(p[1],showaxis=:hide,foreground_color_text=:white,titlefontsize=1)

    display(p)
end




"""



"""
function pcacompwo(dif,VELRES,listpc,facsig)
    l = @layout [grid(2,size(listpc)[1])]
    p = plot(layout=l,legend = false, xformatter=_->"", yformatter=_->"", grid=:false,primary=false,link=:y )#,size=(1000,700),dpi=1000)
    for ix = 1:size(listpc)[1]
        p = plot!(p[ix],dif[:,:,ix],aspect_ratio=:equal,clims=(-facsig*VELRES,facsig*VELRES),seriestype=:heatmap)
    end
    dif = reshape(dif,size(dif)[1]*size(dif)[2],size(dif)[3])
    for ix=size(listpc)[1]+1:2*size(listpc)[1]
        p = histogram!(p[ix],dif[:,ix-4],axisfontsize=5,xlims=(-facsig*VELRES,facsig*VELRES),fill=:white,norm=:true,title="$(listpc[ix-4])PC",titlefontsize=9,tickfontsize=5)#)
    end
    display(p)
end




"""
    contourmap(arr,arraytocontour,lev,xvec,yvec)

Plot the data arr in heatmap type and contours in front at levels=lev.
Contours can be provided from another data.
Arr has to be 2D.
"""
function contourmap(arr,arraytocontour,lev,xvec,yvec)
    plotly()
    p = plot(xvec,yvec,arr,seriestype=:heatmap,aspect_ratio=1)
    p = contour!(xvec,yvec,arraytocontour,levels=lev)
    display(p)
    gr()
end



"""
    corner_cvimap(array1,array2,datadim,labelarray1,labelarray2)

Make a cornerplot between each element of the third dimension of array2, and also
compare with array1.
Array1 has to be 1D (cvi pixel values), and array2 has to be 2D (cvi pixel values
and different cvi maps obtained (with different PC for example)).
Labelarray1 has to be 1 String, Labelarray2 has to be an array os String, with as elements
as the third dimension of array2.
["dollars(ix) PC" for ix=2:4] to label the named of the third dimension of array2 has number of PC used
"""
function corner_cvimap(array1,array2,datadim,labelarray1,labelarray2)
    M = Matrix{Float64}(undef,datadim[1]*datadim[2],size(array2)[2]+1)
    for ix=2:size(array2)[2]
        M[:,ix-1] = array2[:,ix]
    end
    M[:,size(M)[2]] = array1[:]
    push!(labelarray2,labelarray1)
    p = cornerplot(M,label = labelarray2,compact=false,markersize=3,color=:black)
    display(p)
end



"""
    cvi_pdf(data,para0,bin,ylim1,ylim2,alph,Lag,data_name_title,legend ; add=false)

Plot the pdf of the cvi, WIP : and add a gaussian fit on the data. By default, will do a new plot. If add=true, will add this plot on the last figure created.
"""
function cvi_pdf(data,para0,bin,ylim1,ylim2,alph,Lag,data_name_title,title_sup,leg ; add=false)
    gr()

    if add==false
        p = plot(collect(skipmissing(data)),seriestype=:scatterhist_noxerror,alpha=alph,xlims=[-10,10],ylims=[ylim1,ylim2],yaxis=:log,normalized=:pdf,nbins=bin,ylabel=L"p(\Delta V_r)",xlabel=L"\Delta V_r",title="Log-liner plot of the pdf of the longitudinal increments of velocity, lag = r = $(Lag) pix. \n $(data_name_title) \n $(title_sup)",label=leg,titlefontsize=10,framestyle=:box)
    else
        p = plot!(collect(skipmissing(data)),seriestype=:scatterhist_noxerror,alpha=alph,xlims=[-10,10],ylims=[ylim1,ylim2],yaxis=:log,normalized=:pdf,nbins=bin,ylabel=L"p(\Delta V_r)",xlabel=L"\Delta V_r",title="Log-liner plot of the pdf of the longitudinal increments of velocity, lag = r = $(Lag) pix. \n $(data_name_title)\n $(title_sup)",label=leg,titlefontsize=10,framestyle=:box)
    end

    model(t, q) = exp.(-(t.-q[1]).^2/2.0/q[2]^2)./(q[2]sqrt(2pi))
    xar = Array{Float64}(-10:20/size(collect(skipmissing(data)))[1]:10)[1:end-1]
    gaus_fitted = curve_fit(model,xar,collect(skipmissing(data)),para0)
    par = gaus_fitted.param
    p = plot!(xar,model(xar,par),yaxis=:log,xlims=[-10,10])
    display(p) #,ylims=[0.0001,1000000000]
    return par,xar
end



"""
    cvi_pdf_norm(data,sigm,bin,ylim1,ylim2,alph,Lag,data_name_title,legend ; add=false)

Plot the pdf of the cvi, WIP : and add a gaussian fit on the data. By default, will do a new plot. If add=true, will add this plot on the last figure created.
"""
function cvi_pdf_norm(dat::UnivariateKDE,xlim,ylim,alph,Lag,data_name_title,leg ; add=false)
    gr()
    mu       = moment(dat.x,1,aweights(dat.density),0)
    if add==false
        p = plot(dat.x.-mu,dat.density,seriestype=:scatter,alpha=alph,xlims=[xlim[1],xlim[2]],ylims=[ylim[1],ylim[2]],yaxis=:log,ylabel=L"p(\Delta V_r)",xlabel=L"\Delta V_r / \sigma(\Delta V_r)",title="Log-liner pdf of the longitudinal \n increments of velocity at lag l. \n $(data_name_title) ",label=leg,titlefontsize=10,framestyle=:box)
    else
        p = plot!(dat.x.-mu,dat.density,seriestype=:scatter,alpha=alph,xlims=[xlim[1],xlim[2]],ylims=[ylim[1],ylim[2]],yaxis=:log,ylabel=L"p(\Delta V_r)",xlabel=L"\Delta V_r / \sigma(\Delta V_r)",title="Log-liner pdf of the longitudinal \n increments of velocity at lag l. \n $(data_name_title)",label=leg,titlefontsize=10,framestyle=:box)
    end
    display(p) #,ylims=[0.0001,1000000000]
end




function distribcv_multipc(mom1,mom2,mom3,mom4,metric,xvector)
    gr()
    l = @layout [grid(2,2);a{0.4h} ]
    p = plot(layout=l,legend = false,  grid=:true,link=:x,leftmargins=0.3cm,labelfontsize=7)#,primary=false)#,link=:y )#,size=(1000,700),dpi=1000)

    #ytick = collect((minimum(mom1):size(mom1)[1]/10*(maximum(mom1)-minimum(mom1))/size(mom1)[1]:maximum(mom1)))
    #yticklabel = [sci_not(x,2) for x in ytick]
    p = plot!(p[1],xvector,mom1,st=:scatter,shape=:cross,ms=1.5,ylabel=L"\mu",c=:black,alpha=0.5,tickfontsize=5,yformatter=:scientific,xaxis=(0:20:xvector[end]),yguidefontrotation=-90)
    p = plot!(p[2],xvector,mom2,st=:scatter,shape=:cross,ms=1.5,ylabel=L"\sigma",c=:black,alpha=0.5,tickfontsize=5,yformatter=:scientific,xaxis=(0:20:xvector[end]),yguidefontrotation=-90)#,yticks=(ytick,yticklabel))
    p = plot!(p[3],xvector,mom3,st=:scatter,shape=:cross,ms=1.5,ylabel=L"\gamma",c=:black,alpha=0.5,tickfontsize=5,xaxis=(0:20:xvector[end]),yguidefontrotation=-90)#,xlabel="Number of PC")#,yticks=(ytick,yticklabel))
    p = plot!(p[4],xvector,mom4,st=:scatter,shape=:cross,ms=1.5,ylabel=L"\kappa",c=:black,alpha=0.5,tickfontsize=5,xaxis=(0:20:xvector[end]),yguidefontrotation=-90)#,xlabel="Number of PC")#,yticks=(ytick,yticklabel))

    metric = log10.(metric)
    p = plot!(p[5],xvector,metric,st=:scatter,shape=:cross,ms=1.5,c=:black,alpha=0.5,ylabel="Log(metric)",xlabel=L"Number\ of\ PC",tickfontsize=5,xaxis=(0:10:xvector[end]),minorgrid=true)

    display(p)

end #distribcv_multipc


function distribcv_multiow(mom1,mom2,mom3,mom4,metric,xvector,titl;SIGMAT=0)
    gr()
    l = @layout [A{0.01h};grid(2,2);a{0.4h} ]
    p = plot(layout=l,legend = false,  grid=:true,link=:x,leftmargins=0.3cm,labelfontsize=7)#,primary=false)#,link=:y )#,size=(1000,700),dpi=1000)
    p=plot!(p[1],title=titl,legend = false,  grid=:false,showaxis=:hide,foreground_color_text=:white,titlefontsize=8, bottom_margin = -40Plots.px)
    #ytick = collect((minimum(mom1):size(mom1)[1]/10*(maximum(mom1)-minimum(mom1))/size(mom1)[1]:maximum(mom1)))
    #yticklabel = [sci_not(x,2) for x in ytick]
    p = plot!(p[2],xvector,mom1,st=:scatter,shape=:cross,ms=1.5,ylabel=L"\mu",c=:black,tickfontsize=5,yformatter=:scientific,yguidefontrotation=-90)
    if SIGMAT!=0
        p = plot!(p[3],[xvector[1],xvector[end]],[SIGMAT,SIGMAT],ls=:dash,lc=:orange,lw=0.1,ylabel=L"\mu",c=:black,tickfontsize=5,yformatter=:scientific,yguidefontrotation=-90)
    end
    p = plot!(p[3],xvector,mom2,st=:scatter,shape=:cross,ms=1.5,ylabel=L"\sigma",c=:black,tickfontsize=5,yformatter=:scientific,yguidefontrotation=-90)#,yticks=(ytick,yticklabel))
    p = plot!(p[4],xvector,mom3,st=:scatter,shape=:cross,ms=1.5,ylabel=L"\gamma",c=:black,tickfontsize=5,yguidefontrotation=-90)#,xlabel="Number of PC")#,yticks=(ytick,yticklabel))
    p = plot!(p[5],[xvector[1],xvector[end]],[3,3],ls=:dash,lc=:orange,lw=0.1,ylabel=L"\mu",c=:black,tickfontsize=5,yformatter=:scientific,yguidefontrotation=-90)

    p = plot!(p[5],xvector,mom4,st=:scatter,shape=:cross,ms=1.5,ylabel=L"\kappa",c=:black,tickfontsize=5,yguidefontrotation=-90)#,xlabel="Number of PC")#,yticks=(ytick,yticklabel))

    metric = log10.(metric)
    p = plot!(p[6],xvector,metric,st=:scatter,shape=:cross,ms=1.5,c=:black,ylabel="Log(metric)",xlabel="Size of integration \n (number of velocity canal)",tickfontsize=5,minorgrid=true)

    display(p)

end #distribcv_multiow


"""
    compare_rawpc(xvec1,yvec1,array1,dimens,xvec2,yvec2,arraypc,title1::String,saveall::Bool,field)

Plot multiple files, each with two heatmap : the first is a heatmap of the array1 (same for all file), the second a heatmap of the arraypc which change for each file in order to plot all of its third dimension (Data reconstructed with different number of PC for example). Array1 has to be 2D (Pixel*Pixel), arraypc has to be 3D (Pixel*Pixel*PC). The subplots of each files have the same colorscale.
"""
function heatmap_compare_multiplepc(xvec1,yvec1,array1,dimens,xvec2,yvec2,arraypc,title1::String,saveall::Bool,field)
    for ix=1:size(arraypc)[3]
        graphe.heatmap_subplot_two(xvec1,yvec1,array1,xvec2,yvec2,reshape(arraypc[:,:,ix],dimens[1],dimens[2]),(minimum(x->isnan(x) ? -Inf : x,arraypc[:,:,ix]),maximum(x->isnan(x) ? -Inf : x,arraypc[:,:,ix])),(minimum(x->isnan(x) ? -Inf : x,arraypc[:,:,ix]),maximum(x->isnan(x) ? -Inf : x,arraypc[:,:,ix])),title1,"$(field) | $(steppca*ix-(steppca-1))pc",saveall,plotsdir("$(field)/moments","$(field)_$(ix)PC_raw_$(steppca*ix-(steppca-1))pc.pdf"),field)
    end
end



"""
    heatmap_interacted(data,delta_xvec,delta_yvec)

Will produce a heatmap plot of a 2d dataset, with colorscale changing by an interaction of the user on a slider.
"""
function heatmap_interacted(data,delta_xvec,delta_yvec)
    data = Data_preparation.permcolrow(data)
    scene = Scene(resolution=(800,1000))
    fig = Figure(resolution=(800,1000))
    ax = Axis(fig[1, 1])
    lsgrid = labelslidergrid!(fig,"Color scale sup",[0:0.1:10],formats = [x -> "$x" for s in ["K"]],tellwidth = true, width = Auto(),startvalue=58)
    lsgrid2 = labelslidergrid!(fig,"Color scale inf",[-10:0.1:0],formats = [x -> "$x" for s in ["K"]],tellwidth = true, width = Auto(),startvalue=58)

    fig[2,1] = lsgrid.layout
    fig[3,1] = lsgrid2.layout

    sliderobservables = [s.value for s in lsgrid.sliders]
    sliderobservables2 = [s.value for s in lsgrid2.sliders]
    ax.aspect = DataAspect()
    global colorscale_sup
    global colorscale_inf
    colorscale_inf = 0
    colorscale_sup = 10
    #limits!(ax,delta_xvec[end],delta_xvec[1],delta_yvec[1],delta_yvec[end])
    point = lift(sliderobservables...) do slvalues...
        colorscale_sup = slvalues[1]
        heatmap!(ax,reverse(delta_xvec),delta_yvec,data,colormap=cgrad(:curl),colorrange=(colorscale_inf,colorscale_sup))
    end
    point = lift(sliderobservables2...) do slvalues...
        colorscale_inf = slvalues[1]
        heatmap!(ax,reverse(delta_xvec),delta_yvec,data,colormap=cgrad(:curl),colorrange=(colorscale_inf,colorscale_sup))
    end
    #ax.xreversed = true
    display(fig)
    Data_preparation.wait_for_key(prompt="The colorscale choosen will be used in the figure saved. Press any key when ready to continue.")

    return (colorscale_inf,colorscale_sup)
end



"""
    heatmap_subplot_two(xvec1,yvec1,array1,xvec2,yvec2,array2,lim1,lim2,title1,title2,save,field)

Made a subplot in a heatmap form, using two different data (array1 and array2) with a defined colorbar scale (lim1 and lim2), their x and y coordinates (xvec1,yvec1 and xvec2,yvec2) a title for each plot (title1 and title 2), if the plot need to be saved (save) and where.
"""
function heatmap_subplot_two(xvec1,yvec1,array1,xvec2,yvec2,array2,lim1,lim2,title1::String,title2::String,saveall::Bool,savename,field)
    im = plot(layout=2,leg=true,size=(1000,500))
    plot!(im[1],array1,seriestype=:heatmap,clims=lim1,title = title1,aspect_ratio=:equal,xlims=[20,130],ylims=[25,100])
    plot!(im[2],array2,seriestype=:heatmap,clims=lim2,title = title2,aspect_ratio=:equal,xlims=[20,130],ylims=[25,100])

    # plot!(im[1],xvec1,yvec1,array1,seriestype=:heatmap,clims=lim1,title = title1)
    # plot!(im[2],xvec1,yvec1,array2,seriestype=:heatmap,clims=lim2,title = title2)
    saveall==1 && (savefig(savename))
    return(im)
end



"""
    heatmap_subplot_same(arr,j,velocity_range)

Made a subplot in a heatmap form, using the one data but plotting different velocity canals (velocity_range). J is the position of the first plot.
"""
function heatmap_subplot_same(arr,j,velocity_range)
    p = plot(layout=(4,4),leg=false,axis=nothing)
    for i = velocity_range
        j+=1
        plot!(p[j],arr[:,:,i],seriestype=:heatmap,clims=(-1,25))
    end
    display(p)
end



"""
    ifdiff(spectratokeep,array1,reconstructed_array,velocity_vector,datadim,titl,nbrow,nbcol)

Plot the spectra associate with absolute differences in moment calculation (N moment order or CV) higher than a threshold between two data cube. Plot also the differences of the spectra. Array_moment and reconstructed_moment are two 2D arrays, while array1 and reconstructed_array are 3D data cube (PPV).
"""
function ifdiff(spectratokeep,array1,reconstructed_array,velocity_vector,datadim,titl,nbrow,nbcol)
    gr()
    l     = @layout [a{0.00001h} ; grid(nbrow,nbcol)]
    p     = plot(layout=l,legend = false, size=(1000,700),grid=:false,dpi=1000) #,size=(size(rawdat)[1]*25,size(rawdat)[2]*25)
    titl = titl*" | $(size(spectratokeep)[1]) spectra giving that difference | Red : with PCA"
    p     = plot!(p[1],title=titl,showaxis=:hide,foreground_color_text=:white)
    step = 50
    for ix=2:size(spectratokeep)[1]
        if ix>nbrow*nbcol+1
            println("  ")
            println("Total number of spectra giving this difference : $(size(spectratokeep)[1])")
            break
        end
        #step = 2
        p = plot!(p[ix],velocity_vector[:],array1[spectratokeep[ix-1],:],color=:blue,alpha=0.6)
        p = plot!(p[ix],annotate=(maximum(velocity_vector),maximum(array1[spectratokeep[ix-1],:]),text("$(trunc(moment(velocity_vector,1,aweights(array1[spectratokeep[ix-1],:]),0); sigdigits=3))",10,color=:blue)))
        p = plot!(p[ix],[moment(velocity_vector,1,aweights(array1[spectratokeep[ix-1],:]),0),moment(velocity_vector,1,aweights(array1[spectratokeep[ix-1],:]),0)],[minimum(reconstructed_array[spectratokeep[ix-1],:]),maximum(reconstructed_array[spectratokeep[ix-1],:])],color=:blue)
        p = plot!(p[ix],velocity_vector[:],reconstructed_array[spectratokeep[ix-1],:],color=:red,alpha=0.6)
        p = plot!(p[ix],annotate=(minimum(velocity_vector)+1,maximum(array1[spectratokeep[ix-1],:]),text("$(trunc(moment(velocity_vector,1,aweights(reconstructed_array[spectratokeep[ix-1],:]),0); sigdigits=3))",10,color=:red)))
        p = plot!(p[ix],[moment(velocity_vector,1,aweights(reconstructed_array[spectratokeep[ix-1],:]),0),moment(velocity_vector,1,aweights(reconstructed_array[spectratokeep[ix-1],:]),0)],[minimum(reconstructed_array[spectratokeep[ix-1],:]),maximum(reconstructed_array[spectratokeep[ix-1],:])],color=:red)
        p = plot!(p[ix],velocity_vector[:],array1[spectratokeep[ix-1],:].-reconstructed_array[spectratokeep[ix-1],:],color=:green,seriestype=:sticks,alpha=0.9)

    end
    println("  ")
    println("Total number of spectra giving this difference : $(size(spectratokeep)[1])")
    println("  ")
    println("Total number of spectra : $(datadim[1])")

    display(p)
end



"""
    ifdiff_cvmap(cvmap,specdiff_map,xvec,yvec,lim1,lim2,titl)

Add the specdiff_map on the plot of the cvmap . Specdiff_map can be a map where each pixel missing is a spectrq giving a difference in the CV calculation > to a threshold, qnd the others pixels are equal to 0.
"""
function ifdiff_cvmap(cvmap,specdiff_map,xvec,yvec,lim1,lim2,titl)
    gr()
    p = heatmap(xvec,yvec,cvmap.+specdiff_map,clims=(lim1,lim2),title=titl,titlefontsize=10,aspect_ratio=:equal)
    display(p)

end

"""
    meanpc

Plot two images in one plot : the first is the mean of all the spectrums of a data (using directly the rawdat, will be meaned in the function), the second is the data reconstructed with N PC (npc) obtained by a pca method on these data, using M which is a PCA object type introduced in MultivariateStats package. Need the scales of the color bars (scale1 and scale2), the titles of each plot (title1 and title2), a boolean to know if the plot need to be saved (save), the name of the save if needed (filename), and if the plots needs to be display. Arr can be 2D array (PV) or a 3D array (PPV).
"""
function meanpc(arr,M,npc::Int,scale1,scale2,title1::String,title2::String,save::Bool,filename::String,disp::Bool)
        typeof(size(arr))==Tuple{Int64,Int64,Int64} && (arr = reshape(arr,size(arr)[1]*size(arr)[2],size(arr)[3]))
        pc_map   = zeros(Float64,size(arr)[1]*size(arr)[2])
        pc_map  .= projection(M)[:,npc]
        mean_map = mean(th,dims=2)
        pc_map   = reshape(pc_map,size(rawdat)[1],size(rawdat)[2])
        mean_map = reshape(mean_map,size(rawdat)[1],size(rawdat)[2])
        heatmap_subplot_two(mean_map,pc_map,scale1,scale2,title1,title2,save,filename,disp)
end



"""
    mom_conv(momarray,array_specavr,velocity_vector,pcnumb,momentorder,nbrow,nbcol,velocity_range)

Same as the function moment_byintegration, but all plotted in different subplot in the same file.
Pcnumb is a vector with all PC number you work with.
Momentorder should be equal to the moment number you want to plot (used for the title).
"""
function mom_conv(momarray,array_specavr,velocity_newvector,velocity_vector,pcnumb,momentorder::Int,nbrow,nbcol,DATANAMETITLE;xcientific=false)
    momentorder == 1 ? _mom_conv1(momarray,array_specavr,velocity_newvector,velocity_vector,pcnumb,nbrow,nbcol,DATANAMETITLE;xcientific=false) :
    momentorder == 2 ? _mom_conv2(momarray,array_specavr,velocity_newvector,velocity_vector,pcnumb,nbrow,nbcol,DATANAMETITLE;xcientific=false) :
    momentorder == 3 ? _mom_conv3(momarray,array_specavr,velocity_newvector,velocity_vector,pcnumb,nbrow,nbcol,DATANAMETITLE;xcientific=false) :
    momentorder == 4 ? _mom_conv4(momarray,array_specavr,velocity_newvector,velocity_vector,pcnumb,nbrow,nbcol,DATANAMETITLE;xcientific=false) :
    _mom_conv(momarray,velocity_vector,pcnumb)
end

function _mom_conv1(momarray,array_specavr,velocity_newvector,velocity_vector,pcnumb,nbrow,nbcol,DATANAMETITLE;xcientific=false,disp=false)
    gr()
    velocity_range = abs(velocity_newvector[1]-velocity_newvector[2])
    l = @layout [a{0.00001h} ; grid(nbrow,nbcol)]
    p = plot(layout=l,legend = false, size=(2000,1700),grid=:false,fontsize=25)
    titl = "First moment order | Spectra integrated on $(trunc(velocity_range; sigdigits=2)) km/s \n $(DATANAMETITLE)" 
    if (maximum(momarray[:,1])<=1e-15)
        momarray .*= 1e15
        titl = "First moment order multiply by 1e15 | Spectra integrated on $(trunce(velocity_range; sigdigits=2)) km/s \n $(DATANAMETITLE)"
    end
    p = plot!(p[1],title=titl,showaxis=:hide,foreground_color_text=:white,titlefontsize=25)
        for jx = pcnumb[1]:pcnumb[end]
                factor = abs.(minimum(momarray[:,jx])/minimum(array_specavr))
                ytick = collect(minimum(momarray[:,jx]):5*(maximum(momarray[:,jx])-minimum(momarray[:,jx]))/size(momarray[:,jx])[1]:maximum(momarray[:,jx]))
                #ytick = collect(minimum(momarray[:,jx]):0.1*minimum(momarray[:,jx])/maximum(momarray[:,jx]):maximum(momarray[:,jx]))
                yticklabel = [sci_not(y,1) for y in ytick]
                if xcientific==true
                    xtick = collect(minimum(velocity_newvector):5*(maximum(velocity_newvector)-minimum(velocity_newvector))/size(velocity_newvector)[1]:maximum(velocity_newvector))
                    xticklabel = [sci_not(x,1) for x in xtick]
                    p = plot!(p[jx-pcnumb[1]+2],velocity_newvector,momarray[:,jx],xticks=(xtick,xticklabel),yticks=(ytick,yticklabel),ytickfontsize=15,xtickfontsize=15,titlefontsize=15,color=:dodgerblue1,title="Without PCA - $(jx) PC") #,title="$(pctoplot) PC"xrotation=30,
                else
                    p = plot!(p[jx-pcnumb[1]+2],velocity_newvector,momarray[:,jx],yticks=(ytick,yticklabel),ytickfontsize=15,xtickfontsize=15,titlefontsize=15,color=:dodgerblue1,title="Without PCA - $(jx) PC") #,title="$(pctoplot) PC"xrotation=30,
                end
               p = plot!(twinx(p[jx-pcnumb[1]+2]),velocity_vector,array_specavr,color=:green,alpha=0.6,legend=:false,showaxis=:hide,foreground_color_text=:white)
               #p = plot!(p[jx-pcnumb[1]+2],velocity_vector,array_specavr,color=:green,alpha=0.6,legend=:false,showaxis=:hide,foreground_color_text=:white)
               p = plot!(p[jx-pcnumb[1]+2],velocity_newvector,momarray[:,jx],seriestype=:scatter,color=:navy,markershape=:+,markersize=1.5,xrotation=45,xminorticks=-10:1:2)
               if (jx-pcnumb[1]+2)>=(nbrow*nbcol)-nbcol+2
                p = plot!(p[jx-pcnumb[1]+2],xaxis="Velocity (a.u.)",xlabelfontsize=16)
            end
        end
    if disp==true
      display(p)
    end
end

function _mom_conv2(momarray,array_specavr,velocity_newvector,velocity_vector,pcnumb,nbrow,nbcol,DATANAMETITLE;xcientific=false,disp=false)
    gr()
    velocity_range = abs(velocity_newvector[1]-velocity_newvector[2])
    l = @layout [a{0.00001h} ; grid(nbrow,nbcol)]
    p = plot(layout=l,legend = false, size=(2000,1700),grid=:false,fontsize=25)
    titl = "Second moment order | Spectra integrated on $(trunc(velocity_range; sigdigits=2)) km/s \n $(DATANAMETITLE)"
    if (maximum(momarray[:,1])<=1e-15)
        momarray .*= 1e15
        titl = "Second moment order multiplied by 1e15| Spectra integrated on $(trunc(velocity_range; sigdigits=2)) km/s \n $(DATANAMETITLE)"
    end
    p = plot!(p[1],showaxis=:hide,foreground_color_text=:white,titlefontsize=25)#,title=titl)
    for jx = pcnumb[1]:pcnumb[end]
            factor = abs.(minimum(momarray[:,jx])/minimum(array_specavr))
            ytick = collect(minimum(momarray[:,jx]):5*(maximum(momarray[:,jx])-minimum(momarray[:,jx]))/size(momarray[:,jx])[1]:maximum(momarray[:,jx]))
            yticklabel = [sci_not(y,1) for y in ytick]
            if xcientific==true
                xtick = collect(minimum(velocity_newvector):5*(maximum(velocity_newvector)-minimum(velocity_newvector))/size(velocity_newvector)[1]:maximum(velocity_newvector))
                xticklabel = [sci_not(x,1) for x in xtick]
                p = plot!(p[jx-pcnumb[1]+2],velocity_newvector,momarray[:,jx],xticks=(xtick,xticklabel),yticks=(ytick,yticklabel),ytickfontsize=15,xtickfontsize=15,titlefontsize=15,color=:dodgerblue1,title="Without PCA - $(jx) PC") #,title="$(pctoplot) PC"xrotation=30,
            else
                p = plot!(p[jx-pcnumb[1]+2],velocity_newvector,momarray[:,jx],yticks=(ytick,yticklabel),ytickfontsize=15,xtickfontsize=15,titlefontsize=15,color=:dodgerblue1,title="Without PCA - $(jx) PC") #,title="$(pctoplot) PC"xrotation=30,
            end
            p = plot!(twinx(p[jx-pcnumb[1]+2]),velocity_vector,array_specavr,color=:green,alpha=0.6,legend=:false,showaxis=:hide,foreground_color_text=:white)
            p = plot!(p[jx-pcnumb[1]+2],velocity_newvector,momarray[:,jx],seriestype=:scatter,color=:navy,markershape=:+,markersize=1.5,xrotation=45,xminorticks=-10:1:2)
            if (jx-pcnumb[1]+2)>=(nbrow*nbcol)-nbcol+2
                p = plot!(p[jx-pcnumb[1]+2],xaxis="Velocity (a.u.)",xlabelfontsize=16)
            end
    end
    if disp==true
        display(p)
    end
end


function _mom_conv3(momarray,array_specavr,velocity_newvector,velocity_vector,pcnumb,nbrow,nbcol,DATANAMETITLE;xcientific=false,disp=false)
        gr()
        velocity_range = abs(velocity_newvector[1]-velocity_newvector[2])
        l = @layout [a{0.00001h} ; grid(nbrow,nbcol)]
        p = plot(layout=l,legend = false, size=(2000,1800),grid=:false,fontsize=25)
        titl = "Third moment order | Spectra integrated on $(trunc(velocity_range; sigdigits=2)) km/s \n $(DATANAMETITLE)"
        if (maximum(momarray[:,1])<=1e-15)
            momarray .*= 1e15
            titl = "Third moment order multiplied by 1e15| Spectra integrated on $(trunc(velocity_range; sigdigits=2)) km/s \n $(DATANAMETITLE)"
        end
        p = plot!(p[1],title=titl,showaxis=:hide,foreground_color_text=:white,titlefontsize=25)
            for jx = pcnumb[1]:pcnumb[end]
                    factor = abs.(minimum(momarray[:,jx])/minimum(array_specavr))
                    ytick = collect(minimum(momarray[:,jx]):5*(maximum(momarray[:,jx])-minimum(momarray[:,jx]))/size(momarray[:,jx])[1]:maximum(momarray[:,jx]))
                    yticklabel = [sci_not(y,1) for y in ytick]
                    if xcientific==true
                        xtick = collect(minimum(velocity_newvector):5*(maximum(velocity_newvector)-minimum(velocity_newvector))/size(velocity_newvector)[1]:maximum(velocity_newvector))
                        xticklabel = [sci_not(x,1) for x in xtick]
                        p = plot!(p[jx-pcnumb[1]+2],velocity_newvector,momarray[:,jx],xticks=(xtick,xticklabel),yticks=(ytick,yticklabel),ytickfontsize=15,xtickfontsize=15,titlefontsize=15,color=:dodgerblue1,title="Without PCA - $(jx) PC") #,title="$(pctoplot) PC"xrotation=30,
                    else
                        p = plot!(p[jx-pcnumb[1]+2],velocity_newvector,momarray[:,jx],yticks=(ytick,yticklabel),ytickfontsize=15,xtickfontsize=15,titlefontsize=15,color=:dodgerblue1,title="Without PCA - $(jx) PC") #,title="$(pctoplot) PC"xrotation=30,
                    end
                    p = plot!(twinx(p[jx-pcnumb[1]+2]),velocity_vector,array_specavr,color=:green,alpha=0.6,legend=:false,showaxis=:hide,foreground_color_text=:white)
                    p = plot!(p[jx-pcnumb[1]+2],velocity_newvector,momarray[:,jx],seriestype=:scatter,color=:navy,markershape=:+,markersize=1.5,xrotation=45,xminorticks=-10:1:2)
                    if (jx-pcnumb[1]+2)>=(nbrow*nbcol)-nbcol+2
                        p = plot!(p[jx-pcnumb[1]+2],xaxis="Velocity (a.u.)",xlabelfontsize=16)
                    end
            end

            if disp==true
                display(p)
            end
        end

function _mom_conv4(momarray,array_specavr,velocity_newvector,velocity_vector,pcnumb,nbrow,nbcol,DATANAMETITLE;xcientific=false,disp=false)
    gr()
    velocity_range = abs(velocity_newvector[1]-velocity_newvector[2])
    l = @layout [a{0.00001h} ; grid(nbrow,nbcol)]
    p = plot(layout=l,legend = false, size=(2000,1700),grid=:false,fontsize=25)
    titl = "Fourth moment order | Spectra integrated on $(trunc(velocity_range; sigdigits=2)) km/s \n $(DATANAMETITLE)"
    if (maximum(momarray[:,1])<=1e-15)
        momarray .*= 1e15
        titl = "Fourth moment order multiplied by 1e15| Spectra integrated on $(trunc(velocity_range; sigdigits=2)) km/s\n $(DATANAMETITLE)"
    end
    p = plot!(p[1],showaxis=:hide,foreground_color_text=:white,titlefontsize=25)
    for jx = pcnumb[1]:pcnumb[end]
            factor = abs.(minimum(momarray[:,jx])/minimum(array_specavr))
            ytick = collect(minimum(momarray[:,jx]):5*(maximum(momarray[:,jx])-minimum(momarray[:,jx]))/size(momarray[:,jx])[1]:maximum(momarray[:,jx]))
            yticklabel = [sci_not(y,1) for y in ytick]
            if xcientific==true
                xtick = collect(minimum(velocity_newvector):5*(maximum(velocity_newvector)-minimum(velocity_newvector))/size(velocity_newvector)[1]:maximum(velocity_newvector))
                xticklabel = [sci_not(x,1) for x in xtick]
                p = plot!(p[jx-pcnumb[1]+2],velocity_newvector,momarray[:,jx],xticks=(xtick,xticklabel),yticks=(ytick,yticklabel),ytickfontsize=15,xtickfontsize=15,titlefontsize=15,color=:dodgerblue1,title="Without PCA - $(jx) PC") #,title="$(pctoplot) PC"xrotation=30,
            else
                p = plot!(p[jx-pcnumb[1]+2],velocity_newvector,momarray[:,jx],yticks=(ytick,yticklabel),ytickfontsize=15,xtickfontsize=15,titlefontsize=15,color=:dodgerblue1,title="Without PCA - $(jx) PC") #,title="$(pctoplot) PC"xrotation=30,
            end
            p = plot!(twinx(p[jx-pcnumb[1]+2]),velocity_vector,array_specavr,color=:green,alpha=0.6,legend=:false,showaxis=:hide,foreground_color_text=:white)
            p = plot!(p[jx-pcnumb[1]+2],velocity_newvector,momarray[:,jx],seriestype=:scatter,color=:navy,markershape=:+,markersize=1.5,xrotation=0,xminorticks=-10:1:2)
            if (jx-pcnumb[1]+2)>=(nbrow*nbcol)-nbcol+2
                p = plot!(p[jx-pcnumb[1]+2],xaxis="Velocity (a.u.)",xlabelfontsize=16)
            end
    end
    if disp==true
        display(p)
    end
end


"""
    moment_byintegration(momarray,pcnumb,velocity_vector,momentorder)

Plot the 3rd or 4th moment order of the differences of integrated by sections spectras between original data (without PCA) and reconstructed data with PCA. Pcnumb is a vector with all PC number you work with. Momentorder should be equal to the moment number you want to plot (used for the title).
"""
function moment_byintegration(momarray,pcnumb,velocity_vector,momendorder)
    size(pcnumb)[1]>1 && (pcnumb=pcnumb[end])
    if pcnumb==1
         p = plot(velocity_vector,momarray[:,pcnumb:(pcnumb+2)],label=["$(pcnumb)" "$(pcnumb+1)" "$(pcnumb+2)"])
    elseif pcnumb==2
         p = plot(velocity_vector,momarray[:,(pcnumb-1):(pcnumb+2)],label=["$(pcnumb-1)" "$(pcnumb)" "$(pcnumb+1)" "$(pcnumb+2)"])
    elseif pcnumb==size(momarray)[2]
         p = plot(velocity_vector,momarray[:,(pcnumb-15):pcnumb-1])
    elseif pcnumb==(size(momarray)[2]-1)
         p = plot(velocity_vector,momarray[:,(pcnumb-2):pcnumb+1],label=["$(pcnumb-2)" "$(pcnumb-1)" "$(pcnumb)" "$(pcnumb+1)"])
     else
         p = plot(velocity_vector,momarray[:,(pcnumb-15):(pcnumb-1)],label="")
    end
    new_velocity_increment = abs(velocity_vector[1]-velocity_vector[2])
    p = plot!(title="$(momendorder) moment order of the differences between original data and \n reconstructed data on spectra integrated by $(trunc(new_velocity_increment; sigdigits=3)) km/s sections",xlabels="Velocity",ylabel="Emission",titlefontsize=10,legend=false)
    display(p)
end



"""
    moment_multiplepc(moment_multiplepc,pcnumb,nbrow,nbcol,titl,limx,limy;disp=false)

Display scatter subplots comparing moment maps (N moment order or CV) obtain with different number of PC. NPC is compared with (N-1)PC. A linear function y=x is also plotted. Moment_multiplepc is a 3D array where the third dimension gives different maps calculated with a number of PC different.
"""
function moment_multiplepc(moment_multiplepc,pcnumb,nbrow,nbcol,titl,limx,limy;disp=false)
    gr()
    l     = @layout [a{0.01h} ; grid(nbrow,nbcol)]
    p     = plot(layout=l,legend = false, size=(1000,1000),grid=:false,dpi=1000) #,size=(size(rawdat)[1]*25,size(rawdat)[2]*25)
    p     = plot!(p[1],title=titl,showaxis=:hide,foreground_color_text=:white,seriestype=:scatter,legend = false, grid=:false,link=:both,titlefontsize=12,aspect_ratio=:equal,margin=0mm)
    for jx = pcnumb[1]+1:pcnumb[end]
            dens = kde((moment_multiplepc[:,jx],moment_multiplepc[:,jx-1]),boundary=((limx[1]-4*std(moment_multiplepc[:,jx]),limx[2]+4*std(moment_multiplepc[:,jx])),(limy[1]-4*std(moment_multiplepc[:,jx]),limy[2]+4*std(moment_multiplepc[:,jx]))))
            p = plot!(p[jx-pcnumb[1]+1],dens,title="$(jx-1) PC VS $(jx) PC",titlefontsize=10,xlim=limx,ylim=limy)
            p = plot!(p[jx-pcnumb[1]+1],limx,limy,xlim=limx,ylim=limy,aspect_ratio=:equal)
            p = annotate!(p[jx-pcnumb[1]+1],[limx[end]-(limx[2]-limx[1])/4],[limy[1]+(limx[2]-limx[1])/4],text(L"\mu"*"=$(trunc(mean(moment_multiplepc[:,jx]); sigdigits=2))"*"\n"*L"\sigma"*"=$(trunc(std(moment_multiplepc[:,jx]); sigdigits=2))",:red,12))
            p = annotate!(p[jx-pcnumb[1]+1],[limx[1]+(limx[2]-limx[1])/4],[limy[end]-(limx[2]-limx[1])/4],text(L"\mu"*"=$(trunc(mean(moment_multiplepc[:,jx-1]); sigdigits=2))"*"\n"*L"\sigma"*"=$(trunc(std(moment_multiplepc[:,jx-1]); sigdigits=2))",:blue,12))
    end
    if nbrow==2 
        p = plot!(p[2],ylabel="CV values")
        p = plot!(p[4],xlabel="CV values")
        p = plot!(p[5],xlabel="CV values")
    else
        p = plot!(p[2],ylabel="CV values")
        p = plot!(p[5],ylabel="CV values")
        p = plot!(p[8],ylabel="CV values")
        p = plot!(p[8],xlabel="CV values")
        p = plot!(p[9],xlabel="CV values")
        p = plot!(p[10],xlabel="CV values")
    end
    if disp==true 
        display(p)
    end
end


"""
    moment_multiplepc_withraw(moment_multiplepc,pcnumb,array1,nbrow,nbcol,titl)

Same as moment_multiplepc but compare each moment maps obtain with different number of PC with the original data (without PC). A linear function y=x is also plotted.
"""
function moment_multiplepc_withraw(moment_multiplepc,pcnumb,array1,nbrow,nbcol,titl)
    gr()
    l     = @layout [a{0.00001h} ; grid(nbrow,nbcol)]
    p     = plot(layout=l,legend = false, size=(1000,700),grid=:false,dpi=1000) #,size=(size(rawdat)[1]*25,size(rawdat)[2]*25)
    p     = plot!(p[1],title=titl,showaxis=:hide,foreground_color_text=:white,seriestype=:scatter,legend = false, grid=:false,link=:both,titlefontsize=10)
    for jx = 2:nblig*nbcol
            p = plot!(p[jx],array1[:],moment_multiplepc[:,pcnumb[jx-1]],seriestype=:scatter,markersize=0.3,markercolor=:white,alpha=0.3,markershape=:cross,title="Raw VS $(pctoplot) PC",titlefontsize=10)
            p = plot!(p[jx],moment_multiplepc[:,pcnumb[jx-1]],moment_multiplepc[:,pcnumb[jx-1]])
    end
    display(p)
end



"""
    moment_specific_canals(moment_pca,moment_nopca,data_name,pcmax)

Plot a file with 4 subplots, each showing a moment calculation on the mean intensity of specific canals (generally associated with noise). 'moment_pca' are the moments calculated on data reconstructed from PCA, and 'moment_nopca' the moments calculated on the original data set.
"""
function moment_specific_canals(moment_pca,moment_nopca,pcmax)

    p = plot(layout=grid(2,2), size=(1000,700),dpi=1000,leftmargins=0.3cm,link=:x,ylabelfontsize=14) 
    xtick = collect(0:5:pcmax[end])
    xticklabel = [x for x in xtick]
    p = plot!(p[1],moment_pca[1,:],xtick=(xtick,xticklabel),ylabel=L"\mu",seriestype=:scatter,markershape=:+,markersize=4,label="",formatter = :scientific,minorgrid=true)
    p = plot!(p[2],moment_pca[2,:],xtick=(xtick,xticklabel),ylabel=L"\sigma",seriestype=:scatter,markershape=:+,markersize=4,label="",formatter = :auto,minorgrid=true)
    p = plot!(p[3],moment_pca[3,:],xtick=(xtick,xticklabel),ylabel=L"\gamma",seriestype=:scatter,markershape=:+,markersize=4,label="",formatter = :scientific,minorgrid=true)
    p = plot!(p[4],moment_pca[4,:],xtick=(xtick,xticklabel),ylabel=L"\kappa",seriestype=:scatter,markershape=:+,markersize=4,label="",formatter = :auto,minorgrid=true)
    p = plot!(p[1],[0,pcmax[end]],[moment_nopca[1],moment_nopca[1]],label="No PCA",linewidth=2)
    p = plot!(p[2],[0,pcmax[end]],[moment_nopca[2],moment_nopca[2]],label="",linewidth=2)
    p = plot!(p[3],[0,pcmax[end]],[moment_nopca[3],moment_nopca[3]],label="",linewidth=2)
    p = plot!(p[4],[0,pcmax[end]],[moment_nopca[4],moment_nopca[4]],label="",linewidth=2)

    p = plot!(p[3],xlabel="Number of PC")
    p = plot!(p[4],xlabel="Number of PC")
    p = plot!(p[2], yaxis=:log)
    p = plot!(p[4], yaxis=:log)
    display(p)
end


#= TO savedfunction moment_specific_canals(moment_pca,moment_nopca,data_name,pcmax)
    gr()
    l = @layout [a{0.00001h} ; grid(2,2)]
    p = plot(layout=l, size=(1000,900),dpi=1000,leftmargins=0.1cm) #,size=(size(rawdat)[1]*25,size(rawdat)[2]*25)
    p = plot!(p[1],showaxis=:hide,foreground_color_text=:white,seriestype=:scatter, grid=:false,link=:x,titlefontsize=14)#,title="Moments of mean intensity in noise velocity canals \n $(data_name)")
    xtick = collect(0:5:pcmax[end])
    xticklabel = [x for x in xtick]
   @inbounds @views for ix=2:size(moment_pca)[1]+1
        ytick = collect(minimum(moment_pca[ix-1,:]):abs(5*(maximum(moment_pca[ix-1,:])-minimum(moment_pca[ix-1,:]))/size(moment_pca[ix-1,:])[1]):maximum(moment_pca[ix-1,:,]))
        yticklabel = [sci_not(y,2) for y in ytick]        
        p = plot!(p[ix],moment_pca[ix-1,:],xticks=(xtick,xticklabel),yticks=(ytick,yticklabel),ylabel="Moment order : $(ix-1) ",seriestype=:scatter,markershape=:+,markersize=2.5,label="",minorgrid=true,titlefontsize=12)
        (moment_nopca[ix-1]<=maximum(moment_pca[ix-1,:])) && (moment_nopca[ix-1]>=minimum(moment_pca[ix-1,:])) && (p = plot!(p[ix],[0,pcmax[end]],[moment_nopca[ix-1],moment_nopca[ix-1]],label="Without PCA",legend=:bottomright))
        (moment_nopca[ix-1]>maximum(moment_pca[ix-1,:])) && (p = plot!(p[ix],annotate=(trunc(Int,pcmax[end]-pcmax[end]/5),minimum(moment_pca[ix-1,:])+(maximum(moment_pca[ix-1,:])-minimum(moment_pca[ix-1,:]))/6,"Without PCA : $(sci_not(trunc(moment_nopca[ix-1]; sigdigits=3),2))",10)))
        (moment_nopca[ix-1]<minimum(moment_pca[ix-1,:])) && (p = plot!(p[ix],annotate=(trunc(Int,pcmax[end]-pcmax[end]/5),minimum(moment_pca[ix-1,:]),"Without PCA : $(sci_not(trunc(moment_nopca[ix-1]; sigdigits=3),2))",10)))

    end
    p = plot!(p[4],xlabel="Number of PC")
    p = plot!(p[5],xlabel="Number of PC")
    p = plot!(p[3], yaxis=:log)
    p = plot!(p[5], yaxis=:log)
    display(p)
end
=#
"""
    pixels_averaged_spectrum(arr,nbrow,nbcol,titl,pixel_range:Tuple,blank,vel1)

Plot the averaged spectrums of a data. The average is doing on a square of (pixel_range[1]) by (pixel_range[2]) pixels.
"""
function pixels_averaged_spectrum(arr,nbrow,nbcol,titl,pixel_range::Tuple,blank,vel1)
    gr()
    l = @layout [a{0.00001h} ; grid(nbrow,nbcol)]
    p = plot(layout=l,legend = false, size=(1000,700),grid=:false,dpi=1000)
    p = plot!(p[1],title=titl,showaxis=:hide,foreground_color_text=:white)
    plotposition = 1
    for jx = 1:nbrow
        for ix = 1:nbcol
            plotposition += 1
            arr_avr = mean(mean(arr[ix:ix+pixel_range[1],jx:jx+pixel_range[2],:],dims=1),dims=2)
            isequal(arr_avr[1,1,1],missing)==1           && (arr_avr[1,1,:].=blank)
            p = plot!(p[plotposition],vel1,arr_avr[1,1,:],linewidth=0.9)
            ix+pixel_range[1] > size(arr)[2] && continue
        end
        jx+pixel_range[2] > size(arr)[1] && continue
    end
    display(p)
end



"""
    pixels_averaged_spectrum(arr,arr2,nbrow,nbcol,titl,pixel_range::Tuple,blank,vel1)

Plot the averaged spectrums of two cube arr and arr2. The average is doing on a square of (pixel_range)  by (pixel_range) pixels.
"""
function pixels_averaged_spectrum(arr,arr2,nbrow,nbcol,titl,pixel_range::Tuple,blank,vel1)
    gr()
    l = @layout [a{0.00001h} ; grid(nbrow,nbcol)]
    p = plot(layout=l,legend = false, size=(1000,700),grid=:false,dpi=1000)
    p = plot!(p[1],title=titl,showaxis=:hide,foreground_color_text=:white)
    plotposition = 1
    for jx = 1:nbrow
        for ix = 1:nbcol
            plotposition += 1
            arr_avr = mean(mean(arr[ix:ix+pixel_range[1],jx:jx+pixel_range[2],:],dims=1),dims=2)
            isequal(arr_avr[1,1,1],missing)==1           && (arr_avr[1,1,:].=blank)
            p = plot!(p[plotposition],vel1,arr_avr[1,1,:],linewidth=0.7,alpha=0.7)
            arr2_avr = mean(mean(arr2[ix:ix+pixel_range[1],jx:jx+pixel_range[2],:],dims=1),dims=2)
            isequal(arr2_avr[1,1,1],missing)==1           && (arr2_avr[1,1,:].=blank)
            p = plot!(p[plotposition],vel1,arr2_avr[1,1,:],linewidth=0.7,color=:red,alpha=0.7)
            ix+pixel_range[1] > size(arr)[2] && continue
        end
        jx+pixel_range[2] > size(arr)[1] && continue
    end
    display(p)
end

"""
    pixels_averaged_spectrum_twodata(arr,reconstructed_array,nbrow,nbcol,titl,pixel_range,blank,vel1,vel2)

Plot the averaged spectrums of two fits : the rawdata and the data reconstructed from a PCA applied to these rawdata, with two velocity vectors different. The average is doing on a square of (pixel_range)  by (pixel_range) pixels.
"""
function pixels_averaged_spectrum_twodata(arr,reconstructed_array,nbrow,nbcol,titl,pixel_range,blank,vel1,vel2)
    gr()
    l = @layout [a{0.00001h} ; grid(nbrow,nbcol)]
    p = plot(layout=l,legend = false, size=(1000,700),grid=:false,dpi=1000)
    p = plot!(p[1],title=titl,showaxis=:hide,foreground_color_text=:white)
    for j = 1:nblig
        for ix = 1:nbrow*nblig
            arr_avr           = mean(mean(arr[ix:ix+pixel_range,jx:jx+pixel_range,:],dims=1),dims=2)
            reconstructed_avr = mean(mean(reconstructed_array[ix:ix*pixel_range,jx:jx*pixel_range,:],dims=1),dims=2)
            isequal(arr_avr[1,1,1],missing)==1           && (arr_avr[1,1,:].=blank)
            isequal(reconstructed_avr[1,1,1],missing)==1 && (reconstructed_avr[1,1,:].=blank)
            p = plot!(p[j,i],vel1,spec[1,1,:],aspect_ratio=:none,linewidth=0.3,linecolor=:red)
            p = plot!(p[j,i],vel2,spec_rec[1,1,:],aspect_ratio=:none,linewidth=0.3,linecolor=:blue)
            ix+pixel_range > size(rawdat)[2] && continue
        end
        jx+pixel_range > size(rawdat)[1] && continue
        display(p)
    end
end



"""
        plpc(Yt,pc)

Plot some principal components (PC) of data analysed by a PCA. Yt is the transformation of the raw data into principal components.
"""
function plpc(Yt,pc::Integer,titl)
        gr()
        if pc<4
            pc==1 && (l = @layout [a{0.00001h} ; grid(1,1)])
            pc==2 && (l = @layout [a{0.00001h} ; grid(1,2)])
            pc==3 && (l = @layout [a{0.00001h} ; grid(2,2)])
            p = plot(layout=l,legend = false, size=(1000,700),grid=:false,dpi=1000)
            p = plot!(p[1],title=titl,showaxis=:hide,foreground_color_text=:white)
            for j = 2:pc+1
                plot!(p[j],Yt[j-1,:],title="PC number $(j-1)",titlefontsize=10,legend=false)
            end
        end
        if pc>=4
            l = @layout [a{0.00001h} ; grid(2,2)]
            p = plot(layout=l,legend = false, size=(1000,700),grid=:false,dpi=1000)
            p = plot!(p[1],title=titl,showaxis=:hide,foreground_color_text=:white)
            for j = 2:5
                plot!(p[j],Yt[j-1,:],title="PC number $(j-1)",titlefontsize=10,legend=false)
            end
        end
        display(p)
end




"""



"""
function ploptiwind(cubedif,sizenotzero,SIGMAT,BLANK,RANGE)
    l = @layout [a{0.0001h};grid(1,2)]
    p = plot(layout=l,legend = false,  grid=:true,link=:x,leftmargins=0.3cm,labelfontsize=7)#,primary=false)#,link=:y )#,size=(1000,700),dpi=1000)
    p=plot!(p[1],title="$(RANGE) velocity canals used for the integration \n Left colorscale : +- 5* NOISE RMS",legend = false,  grid=:false,showaxis=:hide,foreground_color_text=:white,titlefontsize=8, bottom_margin = -40Plots.px)
    nbins= floor(Int,sqrt(size(cubedif)[1]*size(cubedif)[2]) )
    cubedif = Data_preparation.replace_blanktomissing(cubedif,BLANK)
    sizenotzero = Data_preparation.replace_blanktomissing(sizenotzero,BLANK)

    #histogram!(im[1],reshape(cubedif,(size(cubedif)[1]*size(cubedif)[2])),xlims=(-3*SIGMAT,+3*SIGMAT),bins=(size(cubedif)[1]*size(cubedif)[2]))
    #histogram!(im[2],reshape(sizenotzero,(size(cubedif)[1]*size(cubedif)[2])),xlims=(0,100),bins=nbins)
    p = plot!(p[2],cubedif,seriestype=:heatmap,clims=(-5*SIGMAT,+5*SIGMAT),aspect_ratio=:equal,colorbar=true,tickfontsize=5)
    p = plot!(p[3],sizenotzero,seriestype=:heatmap,aspect_ratio=:equal,clims=(0,100),colorbar=true,tickfontsize=5)
    
    display(p)


end


"""
    pratio(M::PCA,ylog::Bool,pc,titl)

Plot the explained percentage of the data per principal component. If ylog is True, the yaxis will be plotted in log. Pc arguments are used to name the plot.
"""
function pratio(M,ylog::Bool,pc,titl)
    explained_percentage = (principalvars(M)./tvar(M))[1:pc] # each variance from each PC is divided by the sum of all pc variances. tvar(M) is the total variance of the observation. principalvars(M) gives the variance of principal components
    tot = cumsum(explained_percentage)
    if ylog==1
        loga=:log
        leg=:bottomleft
    else
        loga=:none
        leg = :right
    end
    p = plot([1:1:size(explained_percentage)[1]],explained_percentage*100, label="p",title=titl,titlefontsize=10)
    p = plot!([1:1:size(explained_percentage)[1]],explained_percentage*100, markercolor=:blue, seriestype=:scatter, markershape=:+,label=false,markersize=2)
    p = plot!([1:1:(size(explained_percentage)[1])],100 .-cumsum(explained_percentage).*100,linestyle=:solid, linewidth=2,yaxis=loga, label="1-CumulSum(p)",xlabel="Number of principal component",ylabel = "Variance percentage explained", framestyle=:box)
    p = plot!(xaxis=(0:10:pc),xrotation=45,ylims=(1e-2,2e2),yminorticks=10) #ylim=[minimum(explained_percentage*100),110],,yminorticks=10,yticks=5
    display(p)
end


"""
    region_spectrum

Plot spectra from a specific region
"""
function region_spectrum()

end



"""
    region_spectrum_twodata

"""
function region_spectrum_twodata(data2D,data2Dsecond,spectra,velocity_newvector,nbrow,nbcol,titl)
    gr()
    l = @layout [a{0.00001h} ; grid(nbrow,nbcol)]
    p = plot(layout=l,legend = false, size=(1000,700),grid=:false,dpi=1000)
    p = plot!(p[1],title=titl,showaxis=:hide,foreground_color_text=:white)
    for ix = 1:nbrow*nbcol
            p = plot!(p[ix+1],velocity_newvector,data2D[spectra[ix],:],ytickfontsize=15,xtickfontsize=15,titlefontsize=15,color=:dodgerblue1,title="",alpha=0.5) #,title="$(pctoplot) PC"xrotation=30,
            p = plot!(p[ix+1],velocity_newvector,data2D[spectra[ix],:],seriestype=:scatter,color=:navy,markershape=:+,markersize=1.5)

            p = plot!(p[ix+1],velocity_newvector,data2Dsecond[spectra[ix],:],ytickfontsize=15,xtickfontsize=15,titlefontsize=15,color=:red2,title="",alpha=0.5) #,title="$(pctoplot) PC"xrotation=30,
            p = plot!(p[ix+1],velocity_newvector,data2Dsecond[spectra[ix],:],seriestype=:scatter,color=:crimson,markershape=:+,markersize=1.5)

    end
    display(p)
end



"""
    sci_not(x,ndec)

Used to transform axis x labelling in scientific notation with ndec the number of decimals.
"""
function sci_not( x , ndec )
  xchar = strip(Formatting.sprintf1("%17.$(ndec)e",x))
  return xchar
end


function dec_not( x , ndec )
  xchar = strip(Formatting.sprintf1("%.$(ndec)f",x))
  return xchar
end



"""
    StcFct()

Plot structure functions of orders p in function of the structure function of order 3. If add is true, the plot will be added on the actual plot. The figure can be saved with the number of PC used during the process in the name by changing pcfinal.

"""
function StcFct(StcFct,ThirdOrderStcFct,OrderP,DataNameTitle ; add=false, save=true,pcfinal=0)
    my_cgrad = :blues#[:red, :yellow, :blue, :pink, :green, :brown]
    labs = string.(OrderP)
    labs = permutedims(labs)
    cols = distinguishable_colors(length(StcFct), dropseed=true)
    if add==true
        p = plot!(ThirdOrderStcFct,[StcFct[ix,:] for ix=1:size(StcFct)[1]],seriestype=:scatter,yaxis=:log,xaxis=:log,title="Structure function Sp(l) plotted against S3(l) for $(DataNameTitle)",titlefontsize=10,labels=labs,xlabel=L"S_3(l)",ylabel=L"S_p(l)",legend=:bottomright)
    else add=false
        p = plot(ThirdOrderStcFct,StcFct[1,:], c = cols[1],seriestype=:scatter,yaxis=:log,xaxis=:log,title="Structure function Sp(l) plotted against S3(l) for $(DataNameTitle)",titlefontsize=10,labels=labs[1],xlabel=L"S_3(l)",ylabel=L"S_p(l)",legend=:bottomright)
        for ix=2:size(StcFct)[1]
            p = plot!(ThirdOrderStcFct,StcFct[ix,:], c = cols[ix],seriestype=:scatter,yaxis=:log,xaxis=:log,title="Structure function Sp(l) plotted against S3(l) for $(DataNameTitle)",titlefontsize=10,labels=labs[ix])
        end
    end
    if save==false
        return
    end
    display(p)
    #pcfinal!=0 && (savefig(plotsdir("$()/structure_function","struct_p$(OrderP[end])_$(pcfinal)pc.pdf")))
    #pcfinal==0 && (savefig(plotsdir("$()/structure_function","struct_p$(OrderP[end])_nopca.pdf")))

end


"""
    StcFctWithFit()

Plot structure functions of orders p in function of the structure function of order 3 with the fit used to obtain the power-law exponent. If add is true, the plot will be added on the actual plot. The figure can be saved with the number of PC used during the process in the name by changing pcfinal.

"""
function StcFctWithFit(StcFct,ThirdOrderStcFct,OrderP,zeta,fitcan,DataNameTitle ; add=false, save=true,pcfinal=0)
    my_cgrad = :blues#[:red, :yellow, :blue, :pink, :green, :brown]
    labs = string.(OrderP)
    labs = permutedims(labs)
    cols = distinguishable_colors(length(StcFct), dropseed=true)
    if add==true
        p = plot!(ThirdOrderStcFct,[StcFct[ix,:] for ix=1:size(StcFct)[1]],seriestype=:scatter,yaxis=:log,xaxis=:log,title="Structure function Sp(l) plotted against S3(l) for $(DataNameTitle)",titlefontsize=10,labels=labs,xlabel=L"S_3(l)",ylabel=L"S_p(l)",legend=:bottomright)
    else add=false
        p = plot(ThirdOrderStcFct,StcFct[1,:], c = cols[1],seriestype=:scatter,yaxis=:log,xaxis=:log,title="Structure function Sp(l) plotted against S3(l) for $(DataNameTitle)",titlefontsize=10,labels=labs[1],xlabel=L"S_3(l)",ylabel=L"S_p(l)",legend=:bottomright)
        p = plot!(ThirdOrderStcFct[fitcan],zeta[1,2].*ThirdOrderStcFct[fitcan].^zeta[1,1] , c = cols[1],legendfontsize=7,yaxis=:log,xaxis=:log, label="",seriestype=:line)
        for ix=2:size(StcFct)[1]
            p = plot!(ThirdOrderStcFct,StcFct[ix,:], c = cols[ix],seriestype=:scatter,yaxis=:log,xaxis=:log,title="Structure function Sp(l) plotted against S3(l) for $(DataNameTitle)",titlefontsize=10,labels=labs[ix])
            p = plot!(ThirdOrderStcFct[fitcan],zeta[ix,2].*ThirdOrderStcFct[fitcan].^zeta[ix,1], c = cols[ix],legendfontsize=7,yaxis=:log,xaxis=:log, label="",seriestype=:line)
        end
    end
    if save==false
        return
    end
    display(p)
    #pcfinal!=0 && (savefig(plotsdir("structure_function","struct_p$(OrderP[end])_fitted_$(pcfinal)pc.pdf")))
    #pcfinal==0 && (savefig(plotsdir("structure_function","struct_p$(OrderP[end])_fitted_nopca.pdf")))
end



"""
    StcFctExponent(zeta,third_order_zeta,OrderP,xlim,ylim,labs,data_name_title,directories_name; add=false,save=true,pcfinal=0)

Plot the exponent of the structure function power-law fit in function of the order p of the structure function. If add is true, the plot will be added on the actual plot. Also add a y=1/3x line on the plot. The figure can be saved with the number of PC used during the process in the name by changing pcfinal.
"""
function StcFctExponent(zeta,ThirdOrderZeta,OrderP,xlim,ylim,labs,DataNameTitle; add=false,save=true,pcfinal=0,markers=:circle)
    gr()
    #xar = Array{Float64}(xlim[1]:trunc(Int,abs(xlim[1]-xlim[end])/size(OrderP)[1]):xlim[end])
    if add==true
        p = plot!(OrderP,zeta[:,1]./ThirdOrderZeta,seriestype=:scatter,alpha=0.5,label=labs,markershape = markers)
    else add=false
        p = plot(OrderP,zeta[:,1]./ThirdOrderZeta,seriestype=:scatter,alpha=0.5,label=labs,markershape = markers)
        p = plot!(OrderP,OrderP./3,seriestype=:line,label="K41",ylims=ylim,xlims=xlim,xlabel="p",ylabel=L"\zeta(p)/\zeta(3)",title="Power-law exponents of the structure function in function of the order p \n $(DataNameTitle)",titlefontsize=10,xaxis=OrderP,legend=:topleft)
        
        p = plot!(OrderP,OrderP./9 .+2 .*(1 .-(2 ./3) .^(OrderP./3)),linestyle=:dash,label="SL94")
        p = plot!(OrderP,OrderP./9 .+ 1 .-(1 ./3) .^(OrderP./3),linestyle=:dash,label="B02")
    end
    if save==false
        return
    end
    display(p)
    #pcfinal!=0 && (savefig(plotsdir("structure_function","exponent_$(pcfinal)pc.pdf")))
    #pcfinal==0 && (savefig(plotsdir("structure_function","exponent_nopca.pdf")))
end



"""
    StcFctExponentWithError()

    Plot the exponent of the structure function power-law fit in function of the order p of the structure function. Compute and plot also the error of the fit. If add is true, the plot will be added on the actual plot. Also add a y=1/3x line on the plot. The figure can be saved with the number of PC used during the process in the name by changing pcfinal.
"""
function StcFctExponentWithError(zeta,ThirdOrderZeta,StcFct,can,OrderP,sig,xlim,ylim,labs,DataNameTitle; save=true, pcfinal=0)
    gr()
    p0 = ([0.3,500.],[0.7,200.],[1.,3.],[1.3,4.],[1.5,6.],[1.7,6.])
    ZetaUncertainty = Matrix{Float64}(undef,size(OrderP)[1],2)
    for ix=1:size(OrderP)[1]
        ZetaUncertainty[ix,1] = Structure_fct.fit_fctstruct(StcFct[ix,can],StcFct[3,can],p0[ix],sig, confinterv=true)[3][1][1]
        ZetaUncertainty[ix,2] = Structure_fct.fit_fctstruct(StcFct[ix,can],StcFct[3,can],p0[ix],sig, confinterv=true)[3][1][2]
    end
    xar = Array{Float64}(xlim[1]:trunc(Int,abs(xlim[1]-xlim[end])/size(OrderP)[1]):xlim[end])
    p = plot(OrderP,zeta[:,1]./ThirdOrderZeta,seriestype=:scatter,alpha=0.5,label=labs,ribbons=(zeta[:,1]./zeta[3,1].-ZetaUncertainty[:,1],abs.(zeta[:,1]./zeta[3,1].-ZetaUncertainty[:,2])))
    p = plot!(xar,xar./3,seriestype=:line,label="1/3",ylims=ylim,xlims=xlim,xlabel="p",ylabel=L"\zeta(p)/\zeta(3)",title="Power-law exponents of the structure function in function of the order p \n $(DataNameTitle)",titlefontsize=10,xaxis=OrderP,legend=:topleft)
    if save==false
        return
    end
    display(p)
    #pcfinal!=0 && (savefig(plotsdir("structure_function","exponent_werror_$(pcfinal)pc.pdf")))
    #pcfinal==0 && (savefig(plotsdir("structure_function","exponent_werror_nopca.pdf")))
end



# ------------- RECIPES FOR PLOT ------------- #
@recipe function f(::Type{Val{:scatterhist_noxerror}}, x, y, z) #, weights = plotattributes[:weights]
    h = Plots._make_hist((y,), plotattributes[:bins], normed = plotattributes[:normalize], weights = plotattributes[:weights] )
    edge = collect(h.edges[1])
    centers = Plots._bin_centers(edge)
    dx = Plots.diff(h.edges[1])/2
    x := centers
    y := h.weights
    seriestype := :scatter
    @series begin
        seriestype := :yerror
        yerror --> sqrt.(h.weights)/100
        (centers, h.weights)
    end
    ()
end
@recipe function f(::Type{Val{:scatterhist}}, x, y, z) #, weights = plotattributes[:weights]
    h = Plots._make_hist((y,), plotattributes[:bins], normed = plotattributes[:normalize], weights = plotattributes[:weights] )
    edge = collect(h.edges[1])
    centers = Plots._bin_centers(edge)
    dx = Plots.diff(h.edges[1])/2
    x := centers
    y := h.weights
    seriestype := :scatter
    @series begin
        seriestype := :xerror
        xerror --> dx
        (centers, h.weights)
    end
    @series begin
        seriestype := :yerror
        yerror --> sqrt.(h.weights)/100
        (centers, h.weights)
    end
    ()
end





end #module




#= TO BE INTRODUCE IN THE CODE

plot(velocity_vector,data_reconstructed3D[150,300,:],label="8 PC")
plot!(velocity_vector,rawdat[150,300,:],label="Raw data")
plot!(velocity_vector,rawdat[150,300,:].-data_reconstructed3D[150,300,:],seriestype=:sticks,label="Differences")
plot!([moment(velocity_vector,1,aweights(data_reconstructed3D[150,300,:]),0),moment(velocity_vector,1,aweights(data_reconstructed3D[150,300,:]),0)],[0,3],color=:blue,label="CV with 8 PC : $(trunc(moment(velocity_vector,1,aweights(data_reconstructed3D[150,300,:]),0); sigdigits=3)) km/s")
plot!([moment(velocity_vector,1,aweights(collect(skipmissing(rawdat[150,300,:]))),0),moment(velocity_vector,1,aweights(collect(skipmissing(rawdat[150,300,:]))),0)],[0,3],color=:red,label="CV without PCA : $(trunc(moment(velocity_vector,1,aweights(collect(skipmissing(rawdat[150,300,:]))),0); sigdigits=3)) km/s")
plot!(xlabel="Velocity (km/s)",ylabel="Emission (K)")
savefig(plotsdir("$(directories_name)/CV","one_spectra_150_300.pdf"))

=#
