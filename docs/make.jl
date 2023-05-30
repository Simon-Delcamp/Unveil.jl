push!(LOAD_PATH,"../src/")
using Pkg
pkg"activate .."


using Documenter, Unveil
using Data_preparation,Data_analysis,Graphic,Functionforcvi,Functionforpca,Spectralwindowopti,Functionforfbm,Structure_fct
#using CVI_light

push!(LOAD_PATH,"../src/")
makedocs(sitename="Unveil.jl Documentation",
        modules = [Unveil],
        authors = "Simon Delcamp",
        pages = [
            "Quick guide" => "Quickguide.md",
            "Main functions of Unveil" => "Unveil.md",
            "Functions used by Unveil and more"=> "Others.md",
            "Detailed API" => "Api.md"
         ],
         #format = Documenter.LaTeX()
         format = Documenter.HTML(prettyurls = false)
)
deploydocs(
    repo = "https://gricad-gitlab.univ-grenoble-alpes.fr/delcamps/unveil",
)
# Documenter can also automatically deploy documentation to gh-pages.
# See "Hosting Documentation" and deploydocs() in the Documenter manual
# for more information.
#deploydocs(
#    repo = "github.com/Simon-Delcamp/Unveil.jl.git",
#    devbranch = "main"
#)