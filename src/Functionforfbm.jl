###################################################################
# Contains functions to produce fBms maps, compute their power spectra and construct 3D cubes from the fBms 2D maps.
###################################################################

module Functionforfbm


include("Data_preparation.jl") #Read and write fits
include("Data_analysis.jl")

using .Data_preparation
using .Data_analysis
using FITSIO, MultivariateStats, StatsPlots, Plots,StatsBase
using FFTW, AbstractFFTs, Random,Distributions, LaTeXStrings
using LsqFit

export fbm2D
export power_spectra
export ppv_from_fbm




"""
    fbm2D(imsize,powerlaw ; disp = false)

Construct a 2D fbm map of size imsize*imsize, following a power spectrum law equal to powerlaw. Can plot the map (disp = true).
"""
function fbm2D(imsize,powerlaw ; disp = false)
    imsize%2==0 && (Np1 = imsize/2 |> Int16)
    imsize%2==1 && (Np1 = (imsize-1)/2 |> Int16)

    xvecfreq = rfftfreq(imsize)
    yvecfreq = fftfreq(imsize) #Vector(0.001:0.0035:+1)
    xxvec = zeros(Float64,imsize,trunc(Int64,imsize/2+1))
    yyvec = similar(xxvec)
    for ix =1:imsize
        xxvec[ix,:] = xvecfreq
    end
    for ix = 1:Np1+1
        yyvec[:,ix] = yvecfreq
    end
    rr = (xxvec.^2 + yyvec.^2).^0.5
    replace!(x->x!=0.0 ? x : NaN,rr)


    phik = rand(Float64,(imsize,Np1 + 1))*2pi #-pi
    phases = cos.(phik) + sin.(phik).*1im

    # Rescale amplitude of the phases
    #phases ./= sqrt.(sum(phases.^2) ./ (size(phases)[1]*size(phases)[2]))


    output = (rr.^(powerlaw/2.)).*phases

    output[1,1] = 0.0
    if imsize % 2 == 0
        output[2:Np1, 1] = conj(output[imsize:-1:Np1+2,1])
        #output[2:Np1, 2:size(output)[2]] .= conj(output[imsize:-1:Np1+2, size(output)[2]:-1:2])
        #output[2:Np1, size(output)[2]] = conj(output[imsize:-1:Np1+2, size(output)[2]])

        output[Np1+1, 1] = real(output[Np1+1, 1]) + 1im * 0.0
        output[Np1+1, size(output)[2]] = real(output[Np1+1,size(output)[2]]) + 1im * 0.0
        output[2:Np1,Np1+1] .= conj(output[imsize:-1:Np1+2,Np1+1])
   else
       output[2:Np1 + 1, 1] = conj(output[imsize:-1:Np1+2, 1])
       #output[2:Np1 + 1, size(output)[2]] = conj(output[imsize:-1:Np1, size(output)[2]])
   end
    # Zero freq components must have no imaginary part to be own conjugate
    #output[1,Np1] = real(output[1,Np1]) +1im*0.0
    #output[1, 2:Np1-1] .= conj(output[1,Np1-1:-1:Np1+1])
    output[1, 1] = real(output[1, 1]) + 1im*0.0
    output[1, size(output)[2]] = real(output[1,size(output)[2]]) + 1im * 0.0



    output = convert(Array{Complex{Float64},2},Data_preparation.permcolrow(output))
    newmap = irfft(output,imsize)
    if disp==true
        p = plot(newmap,yflip=:true,seriestype=:heatmap,aspect_ratio=:equal)
        display(p)
        sleep(1)
    end
    # fourier_img = output
    # println("Fourier image    :",sum(sum(abs.(fourier_img).^2)))
    #
    # fourier_recovered = fft(ifft(fourier_img))
    # println("Fourier recovered:",sum(sum(abs.(fourier_recovered).^2)))
    #
    # test = abs.(fft(newmap))
    return(output,newmap,xvecfreq)
end


# function fbm3D(imsize,powerlaw)
#
# end

"""
    power_spectra(arr,karr,imsize ; fitted=true)

Compute the power spectrum of arr, print and plot it. karr is the frequencies array (pixel^-1). Also add a powerlaw fit by default.
"""
function power_spectra(arr,karr,imsize ; fitted=true)
    Np   = trunc(Int,imsize/2)
    arr  = abs.(fft(arr)).^2
    #return(arr)
    gr()
    p = plot(karr[2:Np],mean(arr,dims=1)[2:Np]./karr[2:Np],seriestype=:scatter,alpha=0.8,minorgrid=true,yaxis=:log,label="",xaxis=:log,xlabel="k (pixel^{-1})")
    display(p)
    if fitted==true
      model(x,xhi) = xhi[2].*x.^xhi[1]    
      FitFbm = curve_fit(model, karr[2:Np], mean(arr,dims=1)[2:Np]./karr[2:Np], [-3.0,10])
      plot!(karr[2:Np], FitFbm.param[2].*karr[2:Np].^FitFbm.param[1],label="Fit $(trunc(FitFbm.param[2],sigdigits=3))*k^$(trunc(FitFbm.param[1],sigdigits=3))")
      display(p)
    end
    num   = log(mean(arr,dims=1)[2]/karr[2])-log(mean(arr,dims=1)[Np]/karr[Np])
    denum = log(karr[2]) - log(karr[Np])
    coeff = num/denum
    println(coeff)
    return(coeff)
end



"""
    ppv_from_fbm(imsize,powerlaw,xarray,sig,stdev,tmax)

Compute a ppv cube, with two gaussian components centered on the values obtained in 2D fbm maps minus the mean of the xarray (== velocity array). A noise with a standard deviation of stdev is added.
"""
function ppv_from_fbm(imsize,powerlaw,xarray,sig,stdev,tmax)
    xarray_inc = abs(xarray[2]-xarray[1])
    map        = fbm2D(imsize,powerlaw)[2].+mean(xarray)
    #wing       = fbm2D(imsize,powerlaw)[2].+mean(xarray)
    distrib    = Normal(0,stdev)
    cube       = zeros(Float64,imsize*imsize,size(xarray)[1])
    for ix=1:imsize*imsize
        wing = Data_analysis.generate_gaussian(xarray,real(map[ix])+(rand()*10+20)*xarray_inc,sig*(rand()+1))[1].*(tmax/(rand()*2+3))
        cube[ix,:] .= Data_analysis.generate_gaussian(xarray,real(map[ix]),sig)[1].*(tmax+rand()-0.5).+rand(distrib,size(xarray)[1]).+wing
    end
    cube = reshape(cube,imsize,imsize,size(xarray)[1])
    return(cube)
end



"""
    ppv_from_fbm2(imsize,powerlaw,xarray,sig,stdev,tmax)

Compute a ppv cube, with one gaussian components centered on the values obtained in 2D fbm maps minus the mean of the xarray (== velocity array). A noise with a standard deviation of stdev is added.
"""
function ppv_from_fbm2(imsize,powerlaw,xarray,sig,stdev,tmax,fitsname)
    xarray_inc = abs(xarray[2]-xarray[1])
    map        = fbm2D(imsize,powerlaw)[2].+mean(xarray)
    #wing       = fbm2D(imsize,powerlaw)[2].+mean(xarray)
    distrib    = Normal(0,stdev)
    cube       = zeros(Float64,imsize*imsize,size(xarray)[1])
    
    for ix=1:imsize*imsize
        wing = Data_analysis.generate_gaussian(xarray,real(map[ix])+(rand()*10+20)*xarray_inc,sig*(rand()+1))[1].*(tmax/(rand()*2+3))
        cube[ix,:] .= Data_analysis.generate_gaussian(xarray,real(map[ix]),sig)[1].*(tmax+rand()-0.5).+rand(distrib,size(xarray)[1]).+wing
    end
    cube = reshape(cube,imsize,imsize,size(xarray)[1])
    Data_preparation.write_fits("/home/delcamps/Data/Simulated/WithWings/fake_0.00.fits","/home/delcamps/Data/Simulated/fBms/$(fitsname).fits",cube,(size(cube)[1],size(cube)[2],size(cube)[3]),-10000,finished=true)
    return(cube)
end


# cube = ppv_from_fbm(256,-3,velocity_vector,0.5,0.02,4)
# p = plot(velocity_vector,cube[100,140,:])
# display(p)




"""
    ppv_from_fbm3(imsize,powerlaw,xarray,sig,stdev,tmax)

Compute a ppv cube, with two gaussian components of maxima equal to the values obtained in 2D fbm maps minus the mean of the xarray (== velocity array). A noise with a standard deviation of stdev is added.  Xarray is the velocity axis wanted (write it as : -20:01:20 for example) xinc is an array (3 of size) given the increments for each components. sig is an array (3 of size) given the dispersions. A random value will also be added.
"""
#function ppv_from_fbm3(imsize,powerlaw,sig,stdev,tmax,fitspath,fitsname,T,xarray,xinc; more=[""])  # DO NOT ADD A RANDOM VALUE IN SUPPLEMENT OF THE POSITION GIVEN BY THE FBM.
#    #xarray = -20:0.1:20
#    xarray_inc = abs(xarray[2]-xarray[1])
#    map        = fbm2D(imsize,powerlaw)[2].+xinc[3]
#    wingmap1        = fbm2D(imsize,powerlaw)[2].+xinc[1]
#    wingmap2        = fbm2D(imsize,powerlaw)[2].+xinc[2]
#    #wing       = fbm2D(imsize,powerlaw)[2].+mean(xarray)
#    distrib    = Normal(0,stdev)
#    cube       = zeros(Float64,imsize*imsize,size(xarray)[1])
#    for ix=1:imsize*imsize
#        #wing1 = Data_analysis.generate_gaussian(xarray,real(wingmap1[ix]),sig*2)[1].*tmax/6
#        #wing2 = Data_analysis.generate_gaussian(xarray,real(wingmap2[ix]),sig*1.5)[1].*tmax/1.05
#        #cube[ix,:] .= Data_analysis.generate_gaussian(xarray,real(map[ix]),sig)[1].*(tmax).+rand(distrib,size(xarray)[1]).+wing1.+wing2
#        ### NEXT THREE LINE ADDED 16-02
#        wing1 = Data_analysis.generate_gaussian(xarray,real(wingmap1[ix]),sig[1]+rand(1)[1]*2)[1].*T[1] #tmax/20
#        wing2 = Data_analysis.generate_gaussian(xarray,real(wingmap2[ix]),sig[2]+rand(1)[1]*2)[1].*T[2]#tmax/1.35
#        cube[ix,:] .= Data_analysis.generate_gaussian(xarray,real(map[ix]),sig[3]+rand(1)[1]*2)[1].*T[3].+rand(distrib,size(xarray)[1]).+wing1.+wing2
#        #cube[ix,:] .= Data_analysis.generate_gaussian(xarray,real(map[ix]),sig)[1].*(tmax).+rand(distrib,size(xarray)[1]).+wing1.+wing2
#        #p = plot(xarray,wing1)
#        #p = plot!(xarray,wing2)
#        #p = plot!(xarray,cube[ix,:])
#        #p = plot!(xarray,Data_analysis.generate_gaussian(xarray,real(map[ix]),sig)[1].*(tmax).+rand(distrib,size(xarray)[1]))
#        #display(p)
#        #return
#    end
#    
#    cube = reshape(cube,imsize,imsize,size(xarray)[1])
#    Data_preparation.write_fits("/home/delcamps/Data/Simulated/fBms_low/Construction/noisefree.fits","$(fitspath)/$(fitsname)",cube,(size(cube)[1],size(cube)[2],size(cube)[3]),-10000,finished=true,more=more)
#    return(cube)
#end


function ppv_from_fbm3(imsize,powerlaw,sig,stdev,tmax,fitspath,fitsname,T,xarray,xinc; more=[""])  # DO NOT ADD A RANDOM VALUE IN SUPPLEMENT OF THE POSITION GIVEN BY THE FBM.
    #xarray = -20:0.1:20
    xarray_inc = abs(xarray[2]-xarray[1])

    fbmmap1    = fbm2D(imsize,powerlaw)[2]
    fbmmap2    = fbm2D(imsize,powerlaw)[2]
    fbmmap3    = fbm2D(imsize,powerlaw)[2]

    fbmmap1    .= (fbmmap1.-minimum(fbmmap1))
    fbmmap2    .= (fbmmap2.-minimum(fbmmap2))
    fbmmap3    .= (fbmmap3.-minimum(fbmmap3))
    fbmmap1    .= fbmmap1./maximum(fbmmap1).*2.5 .+xinc[1]
    fbmmap2    .= fbmmap2./maximum(fbmmap2).*2.5 .+xinc[2]
    fbmmap3    .= fbmmap3./maximum(fbmmap3).*2.5 .+xinc[3]
    #fbmmap1    .= fbmmap1.*3 .+xinc[1]
    #fbmmap2    .= fbmmap2.*3 .+xinc[2]
    #fbmmap3    .= fbmmap3.*3 .+xinc[3]
    #fbmmap4    .= fbmmap4.*3 .+xinc[4]
    distrib    = Normal(0,stdev)
    cube       = zeros(Float64,imsize*imsize,size(xarray)[1])
    for ix=1:imsize*imsize
        wing1 = Data_analysis.generate_gaussian(xarray,real(fbmmap1[ix]),sig[1]+rand(1)[1]/15)[1].*T[1].*(rand(1)[1])/10#tmax/20
        wing2 = Data_analysis.generate_gaussian(xarray,real(fbmmap2[ix]),sig[2]+rand(1)[1]/15)[1].*T[2].*(rand(1)[1])/10#tmax/1.35
        #wing3 = Data_analysis.generate_gaussian(xarray,real(fbmmap3[ix]),sig[3]+rand(1)[1]/15)[1].*T[3].*(rand(1)[1])/10#tmax/1.35
        cube[ix,:] .= Data_analysis.generate_gaussian(xarray,real(fbmmap3[ix]),sig[3]+rand(1)[1]/15)[1].*T[3].*(rand(1)[1])/10 .+rand(distrib,size(xarray)[1]).+wing1.+wing2#.+wing3
        #cube[ix,:] .= Data_analysis.generate_gaussian(xarray,real(map[ix]),sig)[1].*(tmax).+rand(distrib,size(xarray)[1]).+wing1.+wing2
        #p = plot(xarray,wing1)
        #p = plot!(xarray,wing2)
        #p = plot!(xarray,cube[ix,:])
        #p = plot!(xarray,Data_analysis.generate_gaussian(xarray,real(map[ix]),sig)[1].*(tmax).+rand(distrib,size(xarray)[1]))
        #display(p)
        #return
    end
    
    cube = reshape(cube,imsize,imsize,size(xarray)[1])
    Data_preparation.write_fits("/home/delcamps/Data/Simulated/fBms_low/Construction/noisefree.fits","$(fitspath)/$(fitsname)",cube,(size(cube)[1],size(cube)[2],size(cube)[3]),-10000,finished=true,more=more)
    return(cube)
end

"""
    ppv_from_fbm4(imsize,powerlaw,xarray,sig,stdev,tmax)

Compute a ppv cube, with 4 gaussian components of maxima equal to the values obtained in 2D fbm maps minus the mean of the xarray (== velocity array). A noise with a standard deviation of stdev is added. Xarray is the velocity axis wanted (write it as : -20:01:20 for example) xinc is an array (3 of size) given the increments for each components. sig is an array (3 of size) given the dispersions. A random value will also be added.
"""
function ppv_from_fbm4(imsize,powerlaw,sig,stdev,tmax,fitspath,fitsname,T,xarray,xinc; more=[""])  # DO NOT ADD A RANDOM VALUE IN SUPPLEMENT OF THE POSITION GIVEN BY THE FBM.
    #xarray = -20:0.1:20
    xarray_inc = abs(xarray[2]-xarray[1])
    fbmmap1    = fbm2D(imsize,powerlaw)[2]
    fbmmap2    = fbm2D(imsize,powerlaw)[2]
    fbmmap3    = fbm2D(imsize,powerlaw)[2]
    fbmmap4    = fbm2D(imsize,powerlaw)[2]

    fbmmap1    .= (fbmmap1.-minimum(fbmmap1))
    fbmmap2    .= (fbmmap2.-minimum(fbmmap2))
    fbmmap3    .= (fbmmap3.-minimum(fbmmap3))
    fbmmap4    .= (fbmmap4.-minimum(fbmmap4))
    fbmmap1    .= fbmmap1./maximum(fbmmap1).*3.3 .+xinc[1]
    fbmmap2    .= fbmmap2./maximum(fbmmap2).*3.3 .+xinc[2]
    fbmmap3    .= fbmmap3./maximum(fbmmap3).*3.3 .+xinc[3]
    fbmmap4    .= fbmmap4./maximum(fbmmap4).*3.3 .+xinc[4]
    #fbmmap1    .= fbmmap1.*3 .+xinc[1]
    #fbmmap2    .= fbmmap2.*3 .+xinc[2]
    #fbmmap3    .= fbmmap3.*3 .+xinc[3]
    #fbmmap4    .= fbmmap4.*3 .+xinc[4]
    distrib    = Normal(0,stdev)
    cube       = zeros(Float64,imsize*imsize,size(xarray)[1])
    for ix=1:imsize*imsize
        wing1 = Data_analysis.generate_gaussian(xarray,real(fbmmap1[ix]),sig[1]+rand(1)[1]/15)[1].*T[1].*(rand(1)[1])/10#tmax/20
        wing2 = Data_analysis.generate_gaussian(xarray,real(fbmmap2[ix]),sig[2]+rand(1)[1]/15)[1].*T[2].*(rand(1)[1])/10#tmax/1.35
        wing3 = Data_analysis.generate_gaussian(xarray,real(fbmmap3[ix]),sig[3]+rand(1)[1]/15)[1].*T[3].*(rand(1)[1])/10#tmax/1.35
        cube[ix,:] .= Data_analysis.generate_gaussian(xarray,real(fbmmap4[ix]),sig[4]+rand(1)[1]/15)[1].*T[4].*(rand(1)[1])/10 .+rand(distrib,size(xarray)[1]).+wing1.+wing2.+wing3
        #cube[ix,:] .= Data_analysis.generate_gaussian(xarray,real(map[ix]),sig)[1].*(tmax).+rand(distrib,size(xarray)[1]).+wing1.+wing2
        #p = plot(xarray,wing1)
        #p = plot!(xarray,wing2)
        #p = plot!(xarray,cube[ix,:])
        #p = plot!(xarray,Data_analysis.generate_gaussian(xarray,real(map[ix]),sig)[1].*(tmax).+rand(distrib,size(xarray)[1]))
        #display(p)
        #return
    end
    
    cube = reshape(cube,imsize,imsize,size(xarray)[1])
    Data_preparation.write_fits("/home/delcamps/Data/Simulated/fBms_low/Construction/noisefree.fits","$(fitspath)/$(fitsname)",cube,(size(cube)[1],size(cube)[2],size(cube)[3]),-10000,finished=true,more=more)
    return(cube)
end



end
