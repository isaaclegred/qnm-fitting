# Choose a program to download files, current options curl or wget
export downloader="wget"
# Choose which simulation you want to use from the catalog
export Simulation_Number=0305
mkdir "SXSDATA""$Simulation_Number"
cd "SXSDATA""$Simulation_Number"
# Choose which refinement levels you want
for i in 0 1 2 3 4 5 6
do
    mkdir "Lev"$i
    cd "Lev"$i
    if [ $downloader == "wget" ]
    then
        wget https://zenodo.org/record/3301877/files/SXS:BBH:$Simulation_Number/Lev$i/rhOverM_Asymptotic_GeometricUnits_CoM.h5?download=1 -O rhOverM_Asymptotic_GeometricUnits_CoM.h5 -o download_info.txt
    cd ..
    fi
done
