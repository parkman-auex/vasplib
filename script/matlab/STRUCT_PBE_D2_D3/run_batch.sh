(echo 102;echo 2;echo 0.04)|vaspkit >vaspkit.log
WORKDIRR=`pwd`

for Func in PBE D2 D3
do
    cp POSCAR ${Func}/
    cp KPOINTS ${Func}/
    cp POTCAR ${Func}/
    echo $Func > Func
    cp Func ${Func}/
    cp NP ${Func}/
    cp vasprun ${Func}/
    cd ${Func}
    echo ${Func} `date`
    matlab  -nodesktop -nodisplay -nosplash < ${WORKDIRR}/EvsV.m >log 2>err &
    wait
    cp EvsV_${Func}.mat ..
    cd ${WORKDIRR}
done

# PLOT
matlab  -nodesktop -nodisplay -nosplash < plot_EvsV.m>log 2>err&
wait
