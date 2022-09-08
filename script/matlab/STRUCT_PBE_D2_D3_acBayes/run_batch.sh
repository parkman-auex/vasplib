(echo 102;echo 2;echo 0.04)|vaspkit
cp POSCAR PBE
cp POSCAR D2
cp POSCAR D3
cp KPOINTS PBE
cp KPOINTS D2
cp KPOINTS D3
cp POTCAR PBE
cp POTCAR D2
cp POTCAR D3

cd PBE
echo PBE
matlab  -nodesktop -nodisplay -nosplash <EvsV.m >log 2>err &
wait
cp EvsV_PBE.mat ..
cd ..
cd D2
echo D2
matlab  -nodesktop -nodisplay -nosplash <EvsV.m >log 2>err &
wait
cp EvsV_D2.mat ..
cd ..
cd D3
echo D3
matlab  -nodesktop -nodisplay -nosplash <EvsV.m >log 2>err &
wait
cp EvsV_D3.mat ..
cd ..

matlab  -nodesktop -nodisplay -nosplash < plot_EvsV.m>log 2>err&
wait
