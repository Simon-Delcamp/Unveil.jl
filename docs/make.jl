push!(LOAD_PATH,"../src/")
using Pkg
pkg"activate .."


using Documenter, Unveil
using Dataprep,Analysis,Graphic,CVI,PCA,SWO,Structure_functions
#using CVI_light

push!(LOAD_PATH,"../src/")
makedocs(sitename="Unveil.jl Documentation",
        modules = [Unveil],
        authors = "Simon Delcamp",
        pages = [
            "Quick guide" => "Quickguide.md",
            "Main functions of Unveil" => "Unveil.md",
            "Functions used by Unveil"=> "Others.md",
            "Detailed API" => "Api.md"
         ],
         #format = Documenter.LaTeX()
         format = Documenter.HTML(prettyurls = false)
)
#deploydocs(
#    repo = "gricad-gitlab.univ-grenoble-alpes.fr/delcamps/unveil",
#    branch = "main",
#    devbranch = "main"
#)
# Documenter can also automatically deploy documentation to gh-pages.
# See "Hosting Documentation" and deploydocs() in the Documenter manual
# for more information.
#deploydocs(
#    repo = "github.com/Simon-Delcamp/Unveil.jl.git",
#    devbranch = "main"
#)