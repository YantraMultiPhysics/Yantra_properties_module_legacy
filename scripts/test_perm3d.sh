export YantraPath="../../yantra-dev"
export WorkingDir="."
python3 ../perm3d.py --directions="all" --input_image="../test_images/im3d.npy" --dx=1e-9 --P0=2.5e5 --DP=2.5e5 --export_vtk=1 --tol=1e-4
