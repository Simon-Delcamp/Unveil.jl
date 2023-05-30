#cp allinone_forbash.txt temp.txt
#sed -i 's+fitstopath+/home/delcamps/Data/Simulated/fBms/0-10/fbm_0-10.fits+g' temp.txt
#sed -i 's+pathtosave+/home/delcamps/Data/Simulated/fBms/DifDistr/+g' temp.txt
#sed -i 's+directoryname+0-10+' temp.txt
#sed -i 's+firstnoisecanal+140+' temp.txt
#sed -i 's+lastnoisecanal+201+' temp.txt
#mv temp.txt allinone.txt

for fb in 2 4 6 8; do
    echo "$fb"
    for npc in 15; do
        echo "$npc"

        cp pca_forbash.txt temp.txt
        sed -i "s+fitspath+/home/delcamps/Data/Simulated/fBms/0-0${fb}/fbm_0-0${fb}.fits+" temp.txt
        sed -i "s+pathtosave+/home/delcamps/Data/Simulated/fBms/0-0${fb}/+" temp.txt
        sed -i "s+nbpc+${npc}+" temp.txt
        mv temp.txt pca.txt
        cd ../scripts
        julia pca.jl
        cd ../varfiles 

        cp cvi_forbash.txt temp.txt
        sed -i "s+fitspath+/home/delcamps/Data/Simulated/fBms/0-0${fb}/Data/DataReconstructed_${npc}PC.fits+" temp.txt
        sed -i "s+pathtosave+/home/delcamps/Data/Simulated/fBms/0-0${fb}/+" temp.txt
        sed -i "s+nbpc+${npc}+" temp.txt
        mv temp.txt cvi.txt
        cd ../scripts
        julia cvi.jl
        cd ../varfiles 
    done
done
for npc in 15; do
    echo "$npc"

    cp pca_forbash.txt temp.txt
    sed -i "s+fitspath+/home/delcamps/Data/Simulated/fBms/0-10/fbm_0-10.fits+" temp.txt
    sed -i "s+pathtosave+/home/delcamps/Data/Simulated/fBms/0-10/+" temp.txt
    sed -i "s+nbpc+${npc}+" temp.txt
    mv temp.txt pca.txt
    cd ../scripts
    julia pca.jl

    cd ../varfiles 
    cp cvi_forbash.txt temp.txt
    sed -i "s+fitspath+/home/delcamps/Data/Simulated/fBms/0-10/Data/DataReconstructed_${npc}PC.fits+" temp.txt
    sed -i "s+pathtosave+/home/delcamps/Data/Simulated/fBms/0-10/+" temp.txt
    sed -i "s+nbpc+${npc}+" temp.txt
    mv temp.txt cvi.txt
    cd ../scripts
    julia cvi.jl
    cd ../varfiles 

done


