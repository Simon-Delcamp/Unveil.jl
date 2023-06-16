
module Graphic

using Plots,MultivariateStats, Statistics, StatsBase, LinearAlgebra
using Formatting, LaTeXStrings, Colors
using Measures


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
    p = plot!(p[5],xvector,metric,st=:scatter,shape=:cross,ms=1.5,c=:black,alpha=0.5,ylabel="Log(metric)",xlabel=L"Number\ of\ PC",tickfontsize=5,xaxis=(0:10:xvector[end]),minorgrid=true)
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
        p = plot!(OrderP,zeta[:,1]./ThirdOrderZeta,ribbon=zeta[:,3],fillalpha=0.5,label=labs,markershape = markers) 
    else add=false
        println("yo")
        p = plot(OrderP,zeta[:,1]./ThirdOrderZeta,ribbon=zeta[:,3],alpha=0.5,label=labs,markershape = markers,fillalpha=0.5)

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



end #module Graphic