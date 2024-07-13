###################################################################
# Script ading a mask on a datacube. 
###################################################################
#
include("../src/Data_preparation.jl") # Read and write fits


using .Data_preparation

println("")
println(" Path to the variable file ? (txt file containing all the informations relevant to read and work on the data)") 
VARFILEPATH = "../varfiles/addmask.txt"#"/home/delcamps/Prog/Centroid-Velocities-Increment/var_file/Large_data/Blagrave_HI_DF/stc_fct.txt" #readline()

FITSPATH,PATHTOSAVE,UNITVELOCITY,PHYSUNITS,BLANK,NBMASK = read_var_files(VARFILEPATH)
NBMASK==0 && (error("Change the value of the 'add_mask' parameter on your txt file to >0 ."))    

cube,VELOCITYVECTOR,DATADIMENSION,VELOCITYINCREMENT,HEAD = Data_preparation.read_fits_ppv(FITSPATH,UNITVELOCITY ; check=false)
cube = Data_preparation.replace_nantomissing(cube)

# Prepare directories where plots and data will be saved.
(isdir("$(PATHTOSAVE)/Plots/"))==0  && mkdir("$(PATHTOSAVE)/Plots/")
(isdir("$(PATHTOSAVE)/Data/"))==0   && mkdir("$(PATHTOSAVE)/Data/")

# Convert pixel position into physical units
XVEC,YVEC,DELTAXVEC,DELTAYVEC = Data_preparation.pixtocoord(HEAD,PHYSUNITS)


NBMASK!=0 && (cubemasked = Data_preparation.addmask(cube,(0,5),NBMASK,DATADIMENSION,DELTAXVEC,DELTAYVEC,BLANK,PATHTOSAVE))
cubemasked = reshape(cubemasked,DATADIMENSION)
cubemasked = Data_preparation.replace_missingtoblank(cubemasked,BLANK)
println("Saving Fits")
Data_preparation.write_fits("$(FITSPATH)","polaris2008_masked","$(PATHTOSAVE)/Data/",cubemasked,DATADIMENSION,BLANK)
println("Data masked saved in $(PATHTOSAVE)/Data/polaris2008_masked_NumberOfFilesWithTheSameNameAsPrefixe.fits as a fits.")
