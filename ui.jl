heading("Welcome in the module Unveil.jl !")
 row([
        cell(
            class="st-module",
            [
                h6("Global paramters used in almost every modules"),
                Html.form(action = "/", method="POST", ["Enter the path to fits source : ",
                    input(type="text",id="fitsourcepath", placeholder="Fits source path", name="fitsourcepath")]),
                Html.form(action = "/", method="POST", ["Enter the name of the fits source : ",
                    input(type="text",id="fitsname", placeholder="Fits source name", name="fitsname")]),      
                Html.form(action = "/", method="POST", ["Enter the path for futur saved data : ",
                    input(type="text",id="pathtosave", placeholder="Path for saved data", name="pathtosave")]),  
                Html.form(action = "/", method="POST", ["Enter generic name for futur saved data: ",
                    input(type="text",id="savename", placeholder="Name for saved data: ", name="savename")]),  
                Html.form(action = "/", method="POST", ["Enter indice of first channels with noise: ",
                    input(type="number",placeholder="First noise channel", name="noiseu")]),  
                Html.form(action = "/", method="POST", ["Enter indice of second channels with noise: ",
                    input(type="number",placeholder="Second noise channel", name="noises")]),
                Html.form(action = "/", method="POST", ["Blank value in fits file and to be used : ",
                    input(type="number",placeholder="Blank", name="blank")]),
                Html.form(action = "/", method="POST", ["Overwrite files with same names ? ",
                    input(type="text",placeholder="false/true", name="overwrite")]),
                Html.form(action = "/", method="POST", ["Units of the velocity axis (see header, m/s or km/s): ",
                    input(type="txt",placeholder="Units of velocity", name="velunit")])
                

            ]
        )
    #   cell(class="col-md-12", [
    #                            uploader( multiple = true,
    #                                     accept = ".fits",
    #                                     maxfilesize = 1024*1024*1,
    #                                     maxfiles = 3,
    #                                     autoupload = true,
    #                                     hideuploadbtn = true,
    #                                     label = "Upload datasets",
    #                                     nothumbnails = true,
    #                                     style="max-width: 95%; width: 95%; margin: 0 auto;",

    #                                     @on("rejected", :rejected),
    #                                     @on("uploaded", :uploaded)
    #                                    )

    #                           ])
     ])
row([
    cell(
        class="st-module",
        [
            h6("Convergence of the PCA"),
            Html.form(action = "/", method="POST", [ "Maximum number of PCs for the metric : "
                input(type="text",placeholder="PCmax", name="pcmax")]),  
            Html.form(action = "/", method="POST", [
                input(type="button",value="Run ConvPCA",name="goconvpca")]),  
        ]
        )
])
row([
     cell(
          class="st-module",
          [
           h6("PCA"),
           Html.form(action = "/", method="POST", [ "Number of PCs to use : "
                input(type="number",placeholder="Number of PCs", name="npc")]),  
            Html.form(action = "/", method="POST", [
                input(type="button",value="Run PCA",name="gopca")]),  
          ]
         )
     cell(
          class="st-module",
          [
           h6("SWO"),
           Html.form(action = "/", method="POST", [" Plot random spectra treated with SWO ? "
                input(type="texte",placeholder="YES/NO", name="example")]),  
            Html.form(action = "/", method="POST", [
                input(type="button",value="Run SWO",name="goswo")]), 
          ]
         )
    ])
row([
    cell(
        class="st-module",
        [
         h6("CV"),
         Html.form(action = "/", method="POST", [ "Path + name of file to compute CV (PPV, .fits). Can be a reconstructed (PCA) or treated (SWO) file  : "
            input(type="text",placeholder="Fits file to compute CV", name="cvfits")]),  
         Html.form(action = "/", method="POST", [ "Intensity threshold : "
            input(type="number",placeholder="Intensity threshold", name="intresh")]),  
         Html.form(action = "/", method="POST", [
            input(type="button",value="Run CV",name="gocv")]),  
        ]
       )
   cell(
        class="st-module",
        [
         h6("CVI"),
         Html.form(action = "/", method="POST", [ "Path + name of CV file (PP, .fits)  : "
            input(type="text",placeholder="CV fits file", name="cvifits")]),  
         Html.form(action = "/", method="POST", ["Values of lags : "
            input(type="texte",placeholder="lags", name="lags")]),  
        Html.form(action = "/", method="POST", [
            input(type="button",value="Run CVI",name="gocvi")]),  
        ]
       )
    ])
row([
    cell(
        class="st-module",
        [
            h6("Structure functions"),
            Html.form(action = "/", method="POST", [ "Path + name of CVI file (.fits). Should contains all rotations.  : "
                input(type="text",placeholder="CVI fits file", name="cvifits")]),  
            Html.form(action = "/", method="POST", [ "Orders to compute : "
                input(type="text",placeholder="Orders", name="orders")]),  
            Html.form(action = "/", method="POST", [
                input(type="button",value="Run Spl",name="gospl")]),  
        ]
        )
])