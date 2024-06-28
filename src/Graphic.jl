
module Graphic

include("Analysis.jl")

using .Analysis
using Plots,MultivariateStats, Statistics, StatsBase, LinearAlgebra
using Formatting, LaTeXStrings, Colors
using Measures, LsqFit, Makie, CairoMakie, DelimitedFiles, CurveFit

export checkwindowopti
export distribcv_multiow
export distribcv_multipc
export energyspec
export StcFct
export StcFctWithFit
export StcFctExponent




function adjustspless(file,NORD::Int,NCOL::Int,NROW::Int,LAGTOFIT)
    #file = readdlm(DATPATH,comment_char='#',skipstart=1)
    spl = file[2:NORD+1,2:end]
    f = Figure(backgroundcolor = RGBf(0.98, 0.98, 0.98),size = (1000, 700))
    gg = f[1, 1] = GridLayout()
    axs = [Axis(gg[row, col],xtickalign = 1.0,ytickalign = 1.0, aspect = nothing,xscale=log10,yscale=log10,xlabel=L"S_3(l)",ylabel=L"S_p(l)\times S_3(l)^{-\zeta_p/\zeta_3}") for row in 1:NROW, col in 1:NCOL]
    fitted = zeros(NORD,2)
    for ord=1:NORD
        fitted[ord,:] .= CurveFit.power_fit(spl[3,LAGTOFIT],spl[ord,LAGTOFIT])
    end
    
    for row in 1:NROW, col in 1:NCOL
        if (row-1)*NCOL+col<=NORD 
            Makie.scatter!(axs[row,col],spl[3,:],spl[(row-1)*NCOL+col,:].*spl[3,:].^(.-fitted[(row-1)*NCOL+col,2]./fitted[3,2]),marker=:xcross,color=:black)
        # lines!(axs[row,col],spl[3,:],fitted[1].*spl[3,:].^fitted[2])

            if row!=NROW || col!=1 
                Makie.hidedecorations!.(axs[row, col],ticks = false,ticklabels=false,grid=false)
            else
                Makie.hidedecorations!.(axs[row, col],ticks = false,label=false,ticklabels=false,grid=false)
            end
            Makie.text!(axs[row, col],0.5,0.5,text="p=$((row-1)*NCOL+col)",space = :relative)
        else 
            Makie.hidedecorations!.(axs[row, col])
        end
        

    end
    f

end



function checkwindowopti(sourcecube,cubewo,minimap,VELOCITYVECTOR,NBROW,NBCOL;limx=(minimum(VELOCITYVECTOR),maximum(VELOCITYVECTOR)))
    l = @layout [grid(NBROW,NBCOL)]
    p = plot(layout=l,legend = false, xformatter=_->"", yformatter=_->"", grid=:false,primary=false,link=:both )#,size=(1000,700),dpi=1000)
    xsize = size(sourcecube)[1]
    posx = []
    for ix = 1:NBROW*NBCOL
        append!(posx,rand(1:xsize))
        #yy = Array{Float64}(undef,floor(Int,size(VELOCITYVECTOR)[1]/minimap[posx[ix]]))
        #yy .= 0
        #xx = (1:minimap[posx[ix]]:size(VELOCITYVECTOR)[1]) .|> Int
        p = plot!(p[ix],VELOCITYVECTOR,sourcecube[posx[ix],:],color=:navy,xlim=limx,linewidth=0.5)
        #p = plot!(p[ix],VELOCITYVECTOR,cubewo[posx[ix],:],color=:red,title="$(posx[ix])_$(minimap[posx[ix]])",titlefontsize=5,xlim=limx)
        p = plot!(p[ix],VELOCITYVECTOR,cubewo[posx[ix],:],color=:red,title="$(posx[ix])",titlefontsize=5,xlim=limx,linewidth=0.5)
        #p = plot!(p[ix],VELOCITYVECTOR[xx],yy,seriestype=:scatter,markershape=:+,xlim=limx)
    end
    #p = plot!(p[1],showaxis=:hide,foreground_color_text=:white,titlefontsize=1)

    display(p)
end


function convpca(METRICPATH::String,FIRST::Int,BURNALL::Int;BURNMET=BURNALL,NMEAN=3)
    file = readdlm(METRICPATH,skipstart=2)
    allm = zeros(size(file)[1],5)
    allm[:,1] .= file[:,2]  #Moment 1
    allm[:,2] .= file[:,4]  #Moment 3
    allm[:,3] .= file[:,5]  #Moment 4
    allm[:,4] .= sqrt.(allm[:,1].^2 .+(allm[:,2]).^2 .+allm[:,3].^2 )  #Metric 
    allm[:,5] .= file[:,6]  #PC

    SMOO = zeros(convert(Int64,ceil(size(file)[1]/NMEAN)))
    NEWX = similar(SMOO)
    NEWX .= allm[1:NMEAN:end,5]
    range = FIRST:BURNALL
    f = Figure(backgroundcolor = RGBf(0.98, 0.98, 0.98),size = (1000, 700))

    gg = f[1, 1] = GridLayout()
    #axs = [Axis(gg[row,1],xtickalign = 1.0,ytickalign = 1.0, aspect = nothing,yscale=log10) for row in 1:4]
    #hidedecorations!.(axs, ticks = true)


    for row =1:4
        SMOO .= mean.(abs.(allm[1:NMEAN:end,row]))
        
        if row==1
            axs = Axis(gg[row,1], ylabel = "Moment 1",xtickalign = 1.0,ytickalign = 1.0,xminortickalign = 1.0, aspect = nothing,yscale=log10,xminorticksvisible = true, xminorgridvisible = true,yminorticksvisible = true, yminorgridvisible = true,xminorticks= IntervalsBetween(5))
            hidedecorations!.(axs,ticks = false,ticklabels=false,grid=false,label=false,minorticks=false,minorgrid = false)
        elseif row==2 || row==3
            axs = Axis(gg[row,1], ylabel = "Moment $(row+1)",xtickalign = 1.0,ytickalign = 1.0,xminortickalign = 1.0, aspect = nothing,yscale=log10,xminorticksvisible = true, xminorgridvisible = true,yminorticksvisible = true, yminorgridvisible = true,xminorticks= IntervalsBetween(5))
            hidedecorations!.(axs,ticks = false,ticklabels=false,grid=false,label=false,minorticks=false,minorgrid = false)
        else 
            axs = Axis(gg[row,1], ylabel = "Metric", xlabel="Principal Component",xtickalign = 1.0,xminortickalign = 1.0,ytickalign = 1.0, aspect = nothing,yscale=log10,xminorticksvisible = true,yminorticksvisible = true, yminorgridvisible = true, xminorgridvisible = true,xminorticks= IntervalsBetween(5))
            hidedecorations!.(axs,ticks = false,ticklabels=false,grid=false,label=false,minorticks=false,minorgrid = false)
        end
        Makie.vspan!(axs,FIRST,BURNALL,color=(:grey,0.5))
        Makie.scatter!(axs,allm[:,5],abs.(allm[:,row]),markersize=4,marker=:cross,color=:red)
        mea = mean(filter(x->x!=0,allm[FIRST:BURNALL,row]))
        disp = std(filter(x->x!=0,allm[FIRST:BURNALL,row]))
        Makie.hlines!(axs,mea+3*disp,linestyle=:dot,color=:black)
       
        Makie.scatter!(axs,NEWX,SMOO,markersize=10,marker=:star4,color=:black)
        pc = ceil(Int64,BURNMET/NMEAN)
        if pc>size(SMOO)[1]
            pc = size(SMOO)[1]
        end
        nopt = 0
        while (SMOO[pc]<=mea+3*disp) && pc<size(SMOO)[1] && pc>1
            nopt = pc*3-2
            pc -= 1
        end

        Makie.vlines!(axs,nopt,linestyle=:dash,color=:black)
        Makie.text!(axs,0.9,0.9,text="Nopt = $nopt",space = :relative)
        if row==4
            Makie.xlims!(axs,-1,BURNMET)
        else 
            Makie.xlims!(axs,-1,maximum(allm[:,5])+1)
        end
    end
    rowgap!(gg,1)

    return(f)
end 



function convintegrantspl(cvi,NCOL::Int,NROW::Int,WINDSIZE::Int)
    ptnb = convert(Int64,ceil(size(cvi)[1]/WINDSIZE)) 
    f = Figure(backgroundcolor = RGBf(0.98, 0.98, 0.98),size = (1000, 700))
    gg = f[1, 1] = GridLayout()
    axs = [Axis(gg[row, col],xtickalign = 1.0,ytickalign = 1.0, aspect = nothing) for row in 1:NROW, col in 1:NCOL]
    for row in 1:NROW, col in 1:NCOL
        mea = zeros(ptnb)
        mea .= -1000
        dis = similar(mea)
        ttx = similar(mea)
        tx,ty = Analysis.distrcvi(cvi)
        ty = ty.*abs.(tx).^((row-1)*NCOL+col)
        mini = minimum(ty)
        arr = 1
        for arr=1:ptnb
            if arr*WINDSIZE<=size(ty)[1] || (arr-1)*WINDSIZE+1<=size(ty)[1]
                mea[arr] = mean(ty[(arr-1)*WINDSIZE+1:arr*WINDSIZE])
                dis[arr] = std(ty[(arr-1)*WINDSIZE+1:arr*WINDSIZE])
                ttx[arr] = mean(tx[(arr-1)*WINDSIZE+1:arr*WINDSIZE]) 
            elseif arr*WINDSIZE<=size(ty)[1] #(arr-1)*WINDSIZE+1<=size(ty)[1]
                mea[arr] = mean(ty[(arr-1)*WINDSIZE+1:end])
                dis[arr] = std(ty[(arr-1)*WINDSIZE+1:end])
                ttx[arr] = mean(tx[(arr-1)*WINDSIZE+1:end])
            end
        end

        mea = replace(mea,NaN=>0)
        mea = replace(mea,-1000=>0)

        dis = replace(dis,NaN=>0)
        dis = replace(dis,-1000=>0)

        posi = findall(x->x!=0,mea)
        mea = mea[posi]
        dis = dis[posi]
        ttx = ttx[posi]

        Makie.scatter!(axs[row,col],ttx,mea,markersize=5,marker=:ltriangle,color=:blue)
        Makie.errorbars!(axs[row,col],ttx,mea,dis)

        Makie.text!(axs[row, col],0.9,0.9,text="p=$((row-1)*NCOL+col)",space=:relative)
        Makie.hidedecorations!.(axs[row, col],ticks = false,label=false,ticklabels=false)
        # if minimum(mea)!=0
        #     ylims!(axs[row, col],minimum(mea.-dis)-0.5*minimum(mea.-dis),maximum(mea.+dis)+1e-3*maximum(mea.+dis))
        # else
        #     ylims!(axs[row, col],minimum(filter(x->x!=0,mea.-dis))-0.5*minimum(filter(x->x!=0,mea.-dis)),maximum(mea.+dis)+1e-3*maximum(mea.+dis))
        # end

        
        # xlims!(axs[row, col],-10,10)
    end
    Makie.colgap!(gg,5)
    Makie.rowgap!(gg,5)
    return(f)
end



function convpca(METRICPATH::String;NMEAN=3)
    file = readdlm(METRICPATH,skipstart=2)
    allm = zeros(size(file)[1],5)
    allm[:,1] .= file[:,2]  #Moment 1
    allm[:,2] .= file[:,4]  #Moment 3
    allm[:,3] .= file[:,5]  #Moment 4
    allm[:,4] .= file[:,1]  #Metric 
    allm[:,5] .= file[:,6]  #PC

    SMOO = zeros(convert(Int64,ceil(size(file)[1]/NMEAN)))
    NEWX = similar(SMOO)
    NEWX .= allm[1:NMEAN:end,5]
    f = Figure(backgroundcolor = RGBf(0.98, 0.98, 0.98),size = (1000, 700))

    gg = f[1, 1] = GridLayout()

    for row =1:4
        SMOO .= mean.(abs.(allm[1:NMEAN:end,row]))
        
        if row==1
            axs = Axis(gg[row,1], ylabel = "Moment 1",xtickalign = 1.0,ytickalign = 1.0,xminortickalign = 1.0, aspect = nothing,yscale=log10,xminorticksvisible = true, xminorgridvisible = true,xminorticks= IntervalsBetween(5))
            hidedecorations!.(axs,ticks = false,ticklabels=false,grid=false,label=false,minorticks=false,minorgrid = false)
        elseif row==2 || row==3
            axs = Axis(gg[row,1], ylabel = "Moment $(row+1)",xtickalign = 1.0,ytickalign = 1.0,xminortickalign = 1.0, aspect = nothing,yscale=log10,xminorticksvisible = true, xminorgridvisible = true,xminorticks= IntervalsBetween(5))
            hidedecorations!.(axs,ticks = false,ticklabels=false,grid=false,label=false,minorticks=false,minorgrid = false)
        else 
            axs = Axis(gg[row,1], ylabel = "Metric", xlabel="Principal Component",xtickalign = 1.0,xminortickalign = 1.0,ytickalign = 1.0, aspect = nothing,yscale=log10,xminorticksvisible = true, xminorgridvisible = true,xminorticks= IntervalsBetween(5))
            hidedecorations!.(axs,ticks = false,ticklabels=false,grid=false,label=false,minorticks=false,minorgrid = false)
        end
        Makie.scatter!(axs,allm[:,5],abs.(allm[:,row]),markersize=4,marker=:cross,color=:red)
        Makie.scatter!(axs,NEWX,SMOO,markersize=10,marker=:star4,color=:black)

        Makie.xlims!(axs,-1,maximum(allm[:,5])+1)
        
    end
    rowgap!(gg,1)

    return(f)
end 



function distreigenimage(eigen,NCOL::Int,NROW::Int)
    #eigen = replace(eigen,-1000=>0)
    f = Figure(backgroundcolor = RGBf(0.98, 0.98, 0.98),size = (1000, 700))
    gg = f[1, 1] = GridLayout()
    axs = [Axis(gg[row, col],xtickalign = 1.0,ytickalign = 1.0, aspect = nothing,yscale=log10) for row in 1:NROW, col in 1:NCOL]
    #hidedecorations!.(axs, ticks = true)
    for row in 1:NROW, col in 1:NCOL
        tt = reshape(eigen[:,:,(row-1)*NCOL+col],size(eigen[:,:,(row-1)*NCOL+col])[1]*size(eigen[:,:,(row-1)*NCOL+col])[2])
        tx,ty = Analysis.calcdistr(tt)
        Makie.scatter!(axs[row, col],tx,ty,markersize=5,marker=:cross,color=:black)
        gau = exp.(.-(tx.-0).^2 ./2)./sqrt.(2pi)./1
        Makie.lines!(axs[row, col],tx,gau,linestyle=:dash)
        Makie.ylims!(axs[row, col],1e-3,1)
        Makie.xlims!(axs[row, col],-6,6)
        if row!=NROW || col!=1
            Makie.hidedecorations!.(axs[row, col],ticks = false)
        end
        Makie.text!(axs[row, col],-5.5,3e-1,text="$((row-1)*NCOL+col)")
    end
    Makie.hidedecorations!.(axs[NROW,1],ticks = false,label=false,ticklabels=false)
    Makie.colgap!(gg,1)
    Makie.rowgap!(gg,1)
    return(f)

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
    minimetr = minimum(metric)
    minipc   = xvector[findall(x->x==minimetr,metric)]
    metric = log10.(metric)
    p = plot!(p[5],xvector,metric,st=:scatter,shape=:cross,ms=1.5,c=:black,alpha=0.5,ylabel="Log(metric)",xlabel=L"PC number",tickfontsize=5,xaxis=(0:10:xvector[end]),minorgrid=true)
    p = plot!(p[5],[minipc,minipc],[log10.(minimetr),log10.(minimetr)],st=:scatter,shape=:cross,ms=1.5,c=:red,alpha=0.5)

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
    minimetr = minimum(metric)
    minipc   = xvector[findall(x->x==minimetr,metric)]

    metric = log10.(metric)
    p = plot!(p[6],xvector,metric,st=:scatter,shape=:cross,ms=1.5,c=:black,ylabel="Log(metric)",xlabel="Size of integration \n (number of velocity canal)",tickfontsize=5,minorgrid=true)
    p = plot!(p[6],[minipc,minipc],[log10.(minimetr),log10.(minimetr)],st=:scatter,shape=:cross,ms=1.5,c=:red,alpha=0.5)


    display(p)

end #distribcv_multiow




function distribmom_multipc(mom1,mom2,mom3,mom4,xvector)
    gr()
    l = @layout [grid(2,2)]
    p = plot(layout=l,legend = false,  grid=:true,link=:x,leftmargins=0.3cm,labelfontsize=7)#,primary=false)#,link=:y )#,size=(1000,700),dpi=1000)

    #ytick = collect((minimum(mom1):size(mom1)[1]/10*(maximum(mom1)-minimum(mom1))/size(mom1)[1]:maximum(mom1)))
    #yticklabel = [sci_not(x,2) for x in ytick]
    p = plot!(p[1],xvector,mom1,st=:scatter,shape=:cross,ms=1.5,ylabel=L"\mu",c=:black,alpha=0.5,tickfontsize=5,yformatter=:scientific,xaxis=(0:20:xvector[end]),yguidefontrotation=-90)
    p = plot!(p[2],xvector,mom2,st=:scatter,shape=:cross,ms=1.5,ylabel=L"\sigma",c=:black,alpha=0.5,tickfontsize=5,yformatter=:scientific,xaxis=(0:20:xvector[end]),yguidefontrotation=-90)#,yticks=(ytick,yticklabel))
    p = plot!(p[3],xvector,mom3,st=:scatter,shape=:cross,ms=1.5,ylabel=L"\gamma",c=:black,alpha=0.5,tickfontsize=5,xaxis=(0:20:xvector[end]),yguidefontrotation=-90)#,xlabel="Number of PC")#,yticks=(ytick,yticklabel))
    p = plot!(p[4],xvector,mom4,st=:scatter,shape=:cross,ms=1.5,ylabel=L"\kappa",c=:black,alpha=0.5,tickfontsize=5,xaxis=(0:20:xvector[end]),yguidefontrotation=-90)#,xlabel="Number of PC")#,yticks=(ytick,yticklabel))

    display(p)

end #distribcv_multipc




"""
   energyspec(powerspec,karr,imsize,PATHTOSAVE;fitted=true,SAVENAME="")

Plot an energy power spectra. karr is the frequencies array (pixel^-1). Also add a powerlaw fit by default.
"""
function energyspec(powerspec,karr,imsize,PATHTOSAVE;fitted=true,SAVENAME="")
    Np   = trunc(Int,imsize/2)
    mid = floor(Int,Np/5)
    gr()
    p = plot(karr[2:Np],powerspec[2:Np]./karr[2:Np],seriestype=:scatter,alpha=0.8,minorgrid=true,yaxis=:log,label="",xaxis=:log,xlabel="k (pixel⁻¹)",ylabel="P(k)")  
    display(p)
    if fitted==true
      model(x,xhi) = xhi[2].*x.^xhi[1]    
      
      FitFbm = LsqFit.curve_fit(model, karr[5:mid], powerspec[5:mid]./karr[5:mid], [-3.0,1])
      #FitFbm = LsqFit.curve_fit(model, karr[5:Np], powerspec[5:Np]./karr[5:Np], [-3.0,1])
      a = "$(trunc(FitFbm.param[2],sigdigits=3))"
      ind ="$(trunc(FitFbm.param[1],sigdigits=3))"
      plot!(karr[5:mid], FitFbm.param[2].*karr[5:mid].^FitFbm.param[1],label=latexstring("$(a)\\times k^{$(ind)}"))
      #FitFbm = curve_fit(model, karr[10:Np], mean(arr,dims=1)[10:Np]./karr[10:Np], [-3.0,10])
      FitFbm = LsqFit.curve_fit(model, karr[20:Np], powerspec[20:Np]./karr[20:Np], [-3.0,1])
      FitFbm = LsqFit.curve_fit(model, karr[40:Np], powerspec[40:Np]./karr[40:Np], [-3.0,1])
      a = "$(trunc(FitFbm.param[2],sigdigits=3))"
      ind ="$(trunc(FitFbm.param[1],sigdigits=3))"
      #plot!(karr[20:Np], FitFbm.param[2].*karr[20:Np].^FitFbm.param[1],label=latexstring("$(a)\\times k^{$(ind)}"))
      plot!(karr[40:Np], FitFbm.param[2].*karr[40:Np].^FitFbm.param[1],label=latexstring("$(a)\\times k^{$(ind)}"))
      display(p)
    end
    Plots.savefig("$PATHTOSAVE/Figures/powerspectra_$(SAVENAME).pdf")
end



function expo(file,NORD::Int,LAGTOFIT,LAG)
    #file = readdlm(DATPATH,comment_char='#',skipstart=1)
    spl = file[2:NORD+1,2:end]
    f = Figure(backgroundcolor = RGBf(0.98, 0.98, 0.98),size = (1000, 700))
    ax = Axis(f[1, 1],xlabel=L"p",ylabel=L"\zeta_p/\zeta_3")
    #axs = Axis(xtickalign = 1.0,ytickalign = 1.0, aspect = nothing,xscale=log10,yscale=log10,xlabel=L"S_3(l)",ylabel=L"S_p(l)\times S_3(l)^{-\zeta_p/\zeta_3}")
    fitted = zeros(NORD,2)
    ordd = zeros(NORD)
    res = similar(ordd)
    for ord=1:NORD
        fitted[ord,:] .= CurveFit.power_fit(spl[3,LAGTOFIT],spl[ord,LAGTOFIT])
        eps = zeros(size(spl[3,LAGTOFIT])[1])
        for ix=1:size(spl[3,LAGTOFIT])[1]
            eps[ix] = (spl[ord,ix])-(fitted[ord,1]*(spl[3,ix])^fitted[ord,2])
        end
        res[ord] = sqrt.(sum(eps.^2))./(sum(spl[3,LAGTOFIT]) .-mean(spl[3,LAGTOFIT]).^2)
        # Error based on what I understood quickly of : https://en.wikipedia.org/wiki/Simple_linear_regression
        ordd[ord] = ord
    end
    Makie.scatter!(ordd,fitted[:,2]./fitted[3,2],marker=:xcross,markersize=15,label="[$(LAG[LAGTOFIT[1]+1]):$(LAG[LAGTOFIT[2]+1])]pix")
    Makie.errorbars!(ordd,fitted[:,2]./fitted[3,2],res,whiskerwidth = 10)

    Makie.lines!(ordd,ordd./3,label="K41",color=:black)
    Makie.lines!(ordd,ordd./9 .+2 .*(1 .-(2 ./3).^(ordd./3)),label="SL94",color=:black,linestyle=:dot)
    Makie.lines!(ordd,ordd./9 .+1 .-(1/3).^(ordd./3),label="B02",color=:black,linestyle=:dash)
    return(f)
end

function expoadd(file,NORD::Int,LAGTOFIT,LAG,f)
    spl = file[2:NORD+1,2:end]
    fitted = zeros(NORD,2)
    ordd = zeros(NORD)
    res = similar(ordd)
    for ord=1:NORD
        fitted[ord,:] .= CurveFit.power_fit(spl[3,LAGTOFIT],spl[ord,LAGTOFIT])
        eps = zeros(size(spl[3,LAGTOFIT])[1])
        for ix=1:size(spl[3,LAGTOFIT])[1]
            eps[ix] = (spl[ord,ix])-(fitted[ord,1]*(spl[3,ix])^fitted[ord,2])
        end
        res[ord] = sqrt.(sum(eps.^2))./(sum(spl[3,LAGTOFIT]) .-mean(spl[3,LAGTOFIT]).^2)
        # Error based on what I understood quickly of : https://en.wikipedia.org/wiki/Simple_linear_regression
        ordd[ord] = ord
    end
    Makie.scatter!(ordd,fitted[:,2]./fitted[3,2],marker=:xcross,markersize=15,label="[$(LAG[LAGTOFIT[1]+1]):$(LAG[LAGTOFIT[2]+1])]pix")
    Makie.errorbars!(ordd,fitted[:,2]./fitted[3,2],res,whiskerwidth = 10)
    Makie.axislegend(position = :lt)

    #Makie.axislegend(position = :lt)
    return(f)
end

function integrantspl(CVICUBE,NCOL::Int,NROW::Int)

    f = Figure(backgroundcolor = RGBf(0.98, 0.98, 0.98),size = (1000, 700))
    gg = f[1, 1] = GridLayout()
    axs = [Axis(gg[row, col],xtickalign = 1.0,ytickalign = 1.0, aspect = nothing,yscale=log10) for row in 1:NROW, col in 1:NCOL]
    for row in 1:NROW, col in 1:NCOL
        tx,ty = Analysis.distrcvi(CVICUBE)
        ty = ty.*abs.(tx).^((row-1)*NCOL+col)
        Makie.scatter!(axs[row,col],tx,ty,markersize=5,marker=:ltriangle,color=:blue,alpha=0.5)
        Makie.text!(axs[row, col],0.9,0.9,text="p=$((row-1)*NCOL+col)",space = :relative)
        Makie.hidedecorations!.(axs[row, col],ticks = false,label=false,ticklabels=false)
        Makie.ylims!(axs[row, col],1e-3,1e2)
        Makie.xlims!(axs[row, col],-10,10)
    end
    Makie.colgap!(gg,5)
    Makie.rowgap!(gg,5)
    f
end





function mapeigenimage(eigen,NCOL::Int,NROW::Int,ZLIM)
    data = Dataprep.permcolrow(data)
    eigen = replace(eigen,-1000=>NaN)
    f = Figure(backgroundcolor = RGBf(0.98, 0.98, 0.98),size = (1000, 700))
    gg = f[1, 1] = GridLayout()
    axs = [Axis(gg[row, col],xtickalign = 1.0,ytickalign = 1.0,xticksmirrored = true, xminorticksvisible = true,yticksmirrored = true, yminorticksvisible = true) for row in 1:NROW, col in 1:NCOL]
    for row in 1:NROW, col in 1:NCOL
        Makie.heatmap!( axs[row, col], eigen[:,:,(row-1)*NCOL+col],colorrange=(-ZLIM,ZLIM))
        if row!=NROW || col!=1
            Makie.hidedecorations!.(axs[row, col],ticks = false)
        end
        Makie.text!(axs[row, col],0.8,0.8,text="$((row-1)*NCOL+col)",space = :relative)

    end
    Makie.hidedecorations!.(axs[NROW,1],ticks = false,label=false,ticklabels=false)
    Colorbar(f[1, 2],limits=(-ZLIM,ZLIM))
    #Makie.colgap!(gg,2)
    #Makie.rowgap!(gg,2)
    colsize!(f.layout,1, Aspect(1,1.0))
    resize_to_layout!(f)
    
    return(f)
end # mapeigenimage



"""
    pratio(M::PCA,ylog::Bool,pc,titl)

Plot the explained percentage of the data per principal component. If ylog is True,
the yaxis will be plotted in log.
Pc arguments are used to name the plot.
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
    p = Plots.plot([1:1:size(explained_percentage)[1]],explained_percentage*100, label="p",title=titl,titlefontsize=10)
    p = Plots.plot!([1:1:size(explained_percentage)[1]],explained_percentage*100, markercolor=:blue, seriestype=:scatter, markershape=:+,label=false,markersize=2)
    p = Plots.plot!([1:1:(size(explained_percentage)[1])],100 .-cumsum(explained_percentage).*100,linestyle=:solid, linewidth=2,yaxis=loga, label="1-CumulSum(p)",xlabel="Number of Principal Components",ylabel = "Variance percentage explained", framestyle=:box)
    p = Plots.plot!(xaxis=(0:10:pc),xrotation=45,ylims=(1e-4,2e2),yminorticks=10) #ylim=[minimum(explained_percentage*100),110],,yminorticks=10,yticks=5
    display(p)
end



function metric_PCA(metric,xvector,PATHTOSAVE)
    # plot metric PCA
    X = xvector
    Y = metric
    scatter(X, Y, title="pol_metricPCA", xlabel="PCs", ylabel="Metric", mc=:red, ms=2, ma=0.5)
    savefig(PATHTOSAVE)
end



function residusPCA(PCA,ORI,VELOCITYVECTOR,NCOL::Int,NROW::Int,NCANMAX::Int;file="")
    f = Figure(backgroundcolor = RGBf(0.98, 0.98, 0.98),size = (1000, 700))
    gg = f[1, 1] = GridLayout()
    axs = [Axis(gg[row, col],xtickalign = 1.0,ytickalign = 1.0, aspect = nothing) for row in 1:NROW, col in 1:NCOL]
    if file==""
        open("posispec.dat","w") do io
        end
    else
        temp = readdlm(file)
    end
    yma = maximum(ORI.-PCA)+maximum(ORI.-PCA)*0.5
    ymi = minimum(ORI.-PCA)-abs(minimum(ORI.-PCA))*0.5
    for row in 1:NROW, col in 1:NCOL
        if file==""
            posx = rand(1:size(ORI)[1],1)[1]
            posy = rand(1:size(ORI)[2],1)[1]
            while ORI[posx,posy,1]==-1000 
                posx = rand(1:size(ORI)[1],1)[1]
                posy = rand(1:size(ORI)[2],1)[1]
            end
            open("posispec.dat","a") do io
                towrite = [posx posy]
                writedlm(io,towrite)
            end
        else 
            posx = convert.(Int64,temp[(row-1)*NCOL+col,1])
            posy = convert.(Int64,temp[(row-1)*NCOL+col,2])
        end
        rms = std(ORI[posx,posy,1:NCANMAX])
        Makie.stairs!(axs[row,col],VELOCITYVECTOR,ORI[posx,posy,:].-PCA[posx,posy,:],color=:black)
        Makie.hlines!(axs[row,col],rms,color=:blue)
        Makie.hlines!(axs[row,col],-rms,color=:blue)
        if row!=NROW || col!=1
            Makie.hidedecorations!.(axs[row, col],ticks = false)
        end
        Makie.text!(axs[row, col],8,50,text="$((row-1)*NCOL+col)")

        Makie.xlims!(axs[row, col],VELOCITYVECTOR[1],VELOCITYVECTOR[end])
        Makie.ylims!(axs[row, col],ymi,yma)

    end
    Makie.hidedecorations!.(axs[NROW,1],ticks = false,label=false,ticklabels=false)

    colgap!(gg,2)
    rowgap!(gg,2)
    return(f)
end



function spectrePCASWO(PCA,SWO,ORI,VELOCITYVECTOR,NCOL::Int,NROW::Int;file="")
    f = Figure(backgroundcolor = RGBf(0.98, 0.98, 0.98),size = (1000, 700))
    gg = f[1, 1] = GridLayout()
    axs = [Axis(gg[row, col],xtickalign = 1.0,ytickalign = 1.0) for row in 1:NROW, col in 1:NCOL]
    if file==""
        open("posispec.dat","w") do io
        end
    else
        temp = readdlm(file)
    end
    yma = maximum(ORI)+maximum(ORI)*0.5
    ymi = minimum(ORI)-abs(minimum(ORI))*0.5
    for row in 1:NROW, col in 1:NCOL
        if file==""
            posx = rand(1:size(ORI)[1],1)[1]
            posy = rand(1:size(ORI)[2],1)[1]
            while ORI[posx,posy,1]==-1000 
                posx = rand(1:size(ORI)[1],1)[1]
                posy = rand(1:size(ORI)[2],1)[1]
            end
            open("posispec.dat","a") do io
                towrite = [posx posy]
                writedlm(io,towrite)
            end
        else 
            posx = convert.(Int64,temp[(row-1)*NCOL+col,1])
            posy = convert.(Int64,temp[(row-1)*NCOL+col,2])
        end
        v2 = size(ORI)[3]
        v1 = 1
        while SWO[posx,posy,v1]==0
            v1 += 1
        end
        while SWO[posx,posy,v2]==0
            v2 -= 1
        end
        Makie.stairs!(axs[row,col],VELOCITYVECTOR,ORI[posx,posy,:],color=:grey)
        Makie.stairs!(axs[row,col],VELOCITYVECTOR,PCA[posx,posy,:],color=:blue)
        Makie.vlines!(axs[row,col],VELOCITYVECTOR[v1],color=:orange)
        Makie.vlines!(axs[row,col],VELOCITYVECTOR[v2],color=:orange)
        if row!=NROW || col!=1
            Makie.hidedecorations!.(axs[row, col],ticks = false)
        end
        Makie.text!(axs[row, col],8,50,text="$((row-1)*NCOL+col)")
        Makie.xlims!(axs[row, col],VELOCITYVECTOR[1],VELOCITYVECTOR[end])
        Makie.ylims!(axs[row, col],ymi,yma)

    end
    Makie.hidedecorations!.(axs[NROW,1],ticks = false,label=false,ticklabels=false)
    colgap!(gg,2)
    rowgap!(gg,2)
    return(f)
end



function spless(file,NORD::Int,NCOL::Int,NROW::Int ; LAG=0)
    if LAG==0
       spl = file[2:NORD+1,2:end]

    else 
        spl = file[2:NORD+1,2:LAG]
    end
    f = Figure(backgroundcolor = RGBf(0.98, 0.98, 0.98),size = (1000, 700))
    gg = f[1, 1] = GridLayout()
    axs = [Axis(gg[row, col],xtickalign = 1.0,ytickalign = 1.0, aspect = nothing,xscale=log10,yscale=log10,xlabel=L"S_3(l)",ylabel=L"S_p(l)") for row in 1:NROW, col in 1:NCOL]
    for row in 1:NROW, col in 1:NCOL
        Makie.scatter!(axs[row,col],spl[3,:],spl[(row-1)*NCOL+col,:],marker=:xcross,color=:black)
        Makie.text!(axs[row, col],0.2,0.7,text="p=$((row-1)*NCOL+col)",space=:relative)
        if row!=NROW || col!=1
            Makie.hidedecorations!.(axs[row, col],ticks = false,ticklabels=false,grid=false)
        end
    end
    Makie.hidedecorations!.(axs[NROW,1],ticks = false,label=false,ticklabels=false,grid=false)
    return(f)
end





###############################"
# DEPRECATED



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
    cols = distinguishable_colors(length(StcFct), [RGB(1,1,1), RGB(0,0,0)], dropseed=true)
    if add==true
        p = plot!(ThirdOrderStcFct,[StcFct[ix,:] for ix=1:size(StcFct)[1]],seriestype=:scatter,titlefontsize=10,yaxis=:log,xaxis=:log,labels=labs,xlabel=L"S_3(l)",ylabel=L"S_p(l)",legend=:bottomright)#,title="Structure function Sp(l) plotted against S3(l) for $(DataNameTitle)",
    else add=false
        p = plot(ThirdOrderStcFct,StcFct[1,:], c = cols[1],seriestype=:scatter,titlefontsize=10,labels="p=$(labs[1])",xlabel=L"S_3(l)",ylabel=L"S_p(l)",legend=:bottomright)
        p = plot!(ThirdOrderStcFct[fitcan],zeta[1,2].*ThirdOrderStcFct[fitcan].^zeta[1,1] , c = cols[1],legendfontsize=7, label="",seriestype=:line)
        for ix=2:size(StcFct)[1]
            p = plot!(ThirdOrderStcFct,StcFct[ix,:], c = cols[ix],seriestype=:scatter,titlefontsize=10,labels="p=$(labs[ix])",markershape=Plots.supported_markers()[ix+2],yaxis=:log,xaxis=:log,) #title="Structure function Sp(l) plotted against S3(l) for $(DataNameTitle)",
            p = plot!(ThirdOrderStcFct[fitcan],zeta[ix,2].*ThirdOrderStcFct[fitcan].^zeta[ix,1],yaxis=:log,xaxis=:log, c = cols[ix],legendfontsize=7,seriestype=:line,label="")
        end#yaxis=:log,xaxis=:log,
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
        p = plot!(OrderP,zeta[:,1]./ThirdOrderZeta,ribbon=zeta[:,3],fillalpha=0.5,label=labs,markershape = markers) 
    else add=false
        p = plot(OrderP,zeta[:,1]./ThirdOrderZeta,ribbon=zeta[:,3],alpha=0.5,label=labs,markershape = markers,fillalpha=0.5)

        p = plot!(OrderP,OrderP./3,seriestype=:line,label="K41",ylims=ylim,xlims=xlim,xlabel="p",ylabel=L"\zeta(p)/\zeta(3)",titlefontsize=10,xaxis=OrderP,legend=:topleft) #,title="Power-law exponents of the structure function in function of the order p \n $(DataNameTitle)"
        
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



end #module Graphic