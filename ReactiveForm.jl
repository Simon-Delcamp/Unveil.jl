module ReactiveForm
using GenieFramework
#using .Main.App.StatisticAnalysis
@genietools

@app begin
    @in N = ""
    @out m = 0.0
    @onchange N begin
        m = calc_mean(gen_numbers(N))
    end
end

function glob()
    cell([
        textfield("Path to the source fits", :fitsource)
    ])
end

@page("/reactive", glob)




end