export YantraPath="../../yantra-dev"
export WorkingDir="."
python3 ../diff2d.py --directions="all" --input_image="../test_images/im2d.npy" --dx=1e-6 --D0=1e-9 --export_vtk=1 --tol=1e-4
