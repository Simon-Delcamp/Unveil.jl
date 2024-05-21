module ui
using GenieFramework
using GenieFramework.Genie.Requests: postpayload
@genietools

@app begin
    @in npc = 0
    @in togo = false
    @out aa = 0
    @in gopca = false
  #  @out outtext = ""
    @onchange npc begin
        aa = npc #calc_mean(gen_numbers(N))
    end
    @onbutton gopca begin
        #npc = parse(Int, postpayload(:npc))
        println("PCA running...")
        #Unveil.gopca(npc)
        println("$aa")
        cell(
        p("RECONSTRUCTED DATA SAVED IN {{aa}}"))
     #   outtext = "DATA SAVED IN ..."
    end
end

function fct()
    cell(
          class="st-module",
          [
           h6("PCA"),
           Html.form(action = "/", method="POST", [ "Number of PCs to use : "
                input(type="number",placeholder="Number of PCs", @bind(:npc), name="npc",:npc)]),  
           btn("Run PCA",color="white", textcolor="black",@click("gopca=true"), loading=:gopca,[
                tooltip(contentclass="bg-indigo", contentstyle="font-size: 16px", 
                style="offset: 10px 10px")],type="submit")
            ])
    cell(



    )
end

@page("/", fct)

# @app begin
#     @in npc = 0
#     @onbutton gopca begin
#         println("PCA running...")
#         #Unveil.gopca(npc)
#         println("$npc")
#     end

# end



#@page("/", fct)
end