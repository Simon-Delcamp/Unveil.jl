


# Presentation 
This code unveil positions of extreme intermittent dissipation of turbulence on PPV cubes obtained by observing molecular gas inside molecular clouds. It uses statistical analysis by computing Centroid Velocity Increments, based on Lis et al. 1996, Pety et al. 1999, Hily-Blant et al. 2007, Hily-Blant et al. 2008, Hily-Blant et al. 2009, and on the PhD thesis of Simon Delcamp (in prep). 

One of the original features of the Unveil.jl analysis code is its ability to estimate uncertainties in statistical markers. On the one hand, Unveil.jl provides two independent methods of data treatment: the Spectral Window Optimisation (SWO), inspired by the work of Pety+2003, and a data compression using Principal Component Analysis (PCA). These two methods were tested on numerical simulations in order to calibrate the parameters of use, providing strong recommendations on the SNR required for observations. The code's ability to measure the statistical properties of turbulence was tested on magnetohydrodynamic numerical simulations covering sub- to supersonic and sub- to superAlfvénic regimes.

The idea is to give quick and simple tools to interested researchers, written in a new, promising and well maintain programming language for him to evolve. 

The module **Unveil.jl** is the main one. It contains 8 functions, each one comunicating with others modules inside this package. To be able to use one function of **Unveil**, you should first produce a '.txt' variable file. This allow to use functions with high numbers or long inputs more easily and quickier. Examples of these .txt files can be produced with the function Unveil.prodvarfile().


# Installation
This package is not registered in the General registry of Julia, following the recommended practices (the functionalities of this package will mainly be used by a small group of researchers).

However, the package installation is still simple. While in the local directory where you want to install it, two commands are needed : inside Julia, type :

```
julia> using Pkg
julia> Pkg.add(url="https://gricad-gitlab.univ-grenoble-alpes.fr/delcamps/unveil")
```

<!---The second command will ask you to enter your Username and your Password. For the Username, just enter the one you are using to connect to gitlab (email adress probably). For the Password, you need a "Project Access Tokens". 

*Project Access Tokens* :
Go to the project directory, into Settings then Access Tokens. Here, enter a Token name, remove the expiration date, select a role, check the box "api", and click on create. A new project access token will be prompt above : **save it inside your local computer, it will nether be prompt again !**

After that, everything inside the gitlab repo will be saved in the directory you are working in. --->


# Rules

To be used correctly, few rules should be followed : 
1. All files should be **.fits**
2. For better results, it is highly recommanded to blank the maps edges, and/or the noisiest spectra. 
3. For a PPV (Position-Position-Velocity) cube, at one given position in PP space, if one velocity channel is blanked, it should be blanked on every others velocity channels. We call that *regular blanking*. If not, the code will throw an error. 
4. The fits header should contain these three entries : "BLANK", "BZERO" and "BSCALE". See [this list of keywords definition](https://heasarc.gsfc.nasa.gov/docs/fcg/standard_dict.html) of the fits file convention for more informations. If these entries don't exist, an error will be print.




# Running it


First, in julia, enter : 
```
julia> using Unveil
```

From that, every functions of the package can be used. The first function to use is :
```
julia> Unveil.prodvarfile()
```
It will create all .txt variable files needed to run the other functions. Each function has an associated .txt file where variables are indicated (e.g fits path, values of the lags, ...).  They are used by given these modified txt files as input. For example you want to produce a PCA, just enter :
```
julia> Unveil.pca('path-to-the-txt-file')
```


Numerous others functions of this package can be used. They are inside others modules and should be used accordingly. For example, you can change every "NaN" positions to a given one in a PPV cube : 
```
julia> cube = Unveil.Data_preparation.read_fits_ppv(path,vel_units)[1]   # Read the fits inside 'path'. Need the unit of the veloity dimension
julia> Unveil.Data_preparation.replace_nantoblank(cube,newvalue)   # Change all NaN inside 'cube' into 'newvalue'
```

For a list and description of every function, please see the doc : download the folder 'docs/build' and open the '.html' files inside it. 