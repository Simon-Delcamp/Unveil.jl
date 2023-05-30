# Presentation 
This code is built in order to unveil positions of extreme intermittent dissipation of turbulence on PPV cubes obtained by observing molecular gas inside molecular clouds. Statistical analysis is also performed by computing Centroid Velocity Increments, based on Lis et al. 1996, Pety et al. 1999, Hily-Blant et al. 2007, Hily-Blant et al. 2008, Hily-Blant et al. 2009, and on the PhD thesis of Simon Delcamp (in prep).

The idea is to give quick and simple tools to interested researchers, written in a new, promising and well maintain programming language for him to evolve. 

As can be see on Fig, the module **Unveil.jl** is the main one. It contains 8 functions, each one comunicating with others modules inside this package. To be able to use one function of **Unveil**, you should first produce a '.txt' variable file. This allow to use functions with high numbers or long inputs more easily and quickier. The files should be of the same structure as the one specified in folder 'varfiles' (one file for each Unveil function, with corresponding names). 

## Statistical analysis
The statistical analysis is based on the computation of the Centroid Velocity, and the Centroid Velocity Increment. On a PPV (Position-Position-Velocity) cube, the centroid velocity is writen as :

  $$  C \ \text{[km/s]} = \frac{\sum_{i=1}^{N_v} T(v_i) \times v_i \ \delta v}{\sum_{i=1}^{N_v} T(v_i) \ \delta v}$$

Then, we compute the Centroid Velocity Increment by using a lag noted $l$ :
$$\Delta C(l) = |C(\textbf{r}+l)-C(\textbf{r}) |$$
This computation is made on every pixels of the map. For a given pixel, we compute this for multiple angles around it : the lag is just a distance, not a direction. If only the positions for high intermittent dissipation of turbulence are needed, we do the average for a given pixel. Else, for the statistical analysis, every values are needed. 

However, the centroid velocity computations is very dependant to the noise of the spectra we are working with. Two methods of data pre-treatment are proposed with this code to limit the uncertainties :
1. A Principal Component Analysis (PCA)
2. A Spectral Window Optimisation (SWO)

**PCA** \
The reader should be comfortable with the notions of **Principal Components**, **Data Reconstruction with PCA** and **Matrix Projection**. For reference, see Kendall (1957), Jeffer (1967) and Delcamp PhD thesis (2023). In our code, the idea is to use PCA to enhance the signal-to-noise ratio of our spectra without removing the signals. This procedure is based on the fact that each new PC added in the projection matrix will add a part of the variance of the data. Thus, the signal count for the highest variance, and the noise for the lowest variance. Then, the optimisation of this code ask to decide how many PC will be used to reconstruct the data. 

To do so, we compute the first 4 moments orders of the data projected on each PC. We also compute a metric :
  $$  m_\text{PCA} =  \sqrt{\left(\frac{\mu_\text{PCA}}{\delta v}\right)^2 + \left(\frac{\sigma_\text{PCA}}{\delta v}\right)^2+\gamma_\text{PCA}^2+(\kappa_\text{PCA}-3)^2}$$
When this metric is the lowest, then the distribution is the closest to a Gaussian. If the noise is correctly homogeneous in the cube, this metric should be a Gaussian when the corresponding PC is reproducing only the noise.  Thus we are looking for the lowest value of this metric.



**SWO** \
The idea behind this method is to reduce the size of the spectra by removing the noise channels. This method compute intensity integration of each spectra over a given number of channels (to smooth the intensity). For each new bigger channel, we compute the differences with the precedent one. At the end, we are going from noise onyl to signal part of a spectra when this result is close to $$\sigma_T^2 dv N_\text{T}$$  
with $\sigma_\text{T}$ being the dispersion of the noise, $dv$ the step for integration and $N_\text{T}$ the number of noise channels. However, this computation depends strongly on the step used for the integrations.

It can be demonstrate that if we compute the intensity integration of the original spectra and of the spectra reproduced by this method, if the method is working perfectly, their differences should be equal to an integration of a noise part of the spectra only. Thus, the first 4 moments orders should caracterise a Gaussian with a dispersion equal to the noise dispersion. Then, the given metric is the smallest possible when the method works the best :

$$   m_\text{OLS} =  \sqrt{\left(\frac{\mu_I^{N_T}}{\delta vdI}\right)^2 + \left(\frac{\sigma_I^{N_T}-\sigma_\text{b}}{\delta vdI}\right)^2+\gamma_I^2+(\kappa_I-3)^2} $$

# Where is the doc ?
The doc is in a form of html pages. You can find them in the folder 'docs/build'. Clone them and open them.

# Installation
This package is not registered in the General registry of Julia following the recommended practices. It is clear that the functionalities of this package will mainly by used by a small group of researchers. 

However, the package installation is still simple. While in the local directory where you want to install it, two commands are needed : inside Julia, type :

```
julia> using Pkg
julia> Pkg.add(url="https://gricad-gitlab.univ-grenoble-alpes.fr/delcamps/unveil")
```

The second command will ask you to enter your Username and your Password. For the Username, just enter the one you are using to connect to gitlab (email adress probably). For the Password, you need a "Project Access Tokens". 

*Project Access Tokens* :
Go to the project directory, into Settings then Access Tokens. Here, enter a Token name, remove the expiration date, select a role, check the box "api", and click on create. A new project access token will be prompt above : **save it inside your local computer, it will nether be prompt again !**

After that, everything inside the gitlab repo will be saved in the directory you are working in. 


# Recommendations and rules

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

From that, every functions of the package can be used. The first step to use one of these functions is to modify the '.txt' files correctly (see [Presentation](#Presentation)). Then, if for example you want to produce a PCA, just enter :
```
julia> Unveil.pca('path-to-the-txt-file')
```


Numerous others functions of this package can be used inside others modules. For example, you can change every "NaN" positions to a given one in a PPV cube : 
```
julia> cube = Unveil.Data_preparation.read_fits_ppv(path,vel_units)[1]   # Read the fits inside 'path'. Need the unit of the veloity dimension
julia> Unveil.Data_preparation.replace_nantoblank(cube,newvalue)   # Change all NaN inside 'cube' into 'newvalue'
```

For a list and description of every function, please see the doc : download the folder 'docs/build' and open the '.html' files inside it. 