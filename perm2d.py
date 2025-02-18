"""
computes water permeability from image
"""
#%%Import relvant modules
import os
import sys
import warnings
YantraPath = os.environ['YantraPath']  #defined as environmental variable
WorkingDir = os.environ['WorkingDir'] #defined as environmental variable
sys.path.append(YantraPath)
import yantra
import numpy  as np
try:
    from skimage import io
except ImportError:
    warnings.warn("skimage not found cannot import image of .tif or .tiff format")
    
#%%Argument parser to get input
import argparse
parser = argparse.ArgumentParser(description='Argument parser for permeability solver')
parser.add_argument('--dx',type=float,default=1.,help="voxel size of image in physical units")
parser.add_argument('--input_image',type=str,default='',help="path for input image")
parser.add_argument('--directions',type = str,default='all',help="x,y or all to get permeability in x,y or all directions respectively")
parser.add_argument('--save_yantra_state',type=int,default=0,help="whether to save yantra state or not")
parser.add_argument('--export_vtk',type=int,default=1,help="whether to export velocity to vtk or not")
parser.add_argument('--max_iters',type=int,default=10000,help="max number of iterations allowed to reach steady state")
parser.add_argument('--DP',type=float,default=0.001,help="volumetric force for permeability computation")
parser.add_argument('--P0',type=float,default=1,help="reference pressure for exact comptuation of pressure")    
parser.add_argument('--tol',type=float,default=1e-4,help="tolerence to reach steady state")
parser.add_argument('--array_job',type=int,default=0,help="tolerence to reach steady state")
parser.add_argument('--array_jobID',type=int,default=0,help="tolerence to reach steady state")
parser.add_argument('--array_job_input',type=str,default='',help="tolerence to reach steady state")  
args = parser.parse_args()
dx_phy=args.dx
directions = args.directions
input_image=args.input_image
save_yantra_state = args.save_yantra_state
export_vtk = args.export_vtk
max_iters = args.max_iters
DP = args.DP
Pref = args.P0
tol = args.tol
array_job = args.array_job
array_jobID = args.array_jobID
array_job_input = args.array_job_input
if array_job:
    with open(array_job_input,'r')  as f:
        inputs = f.readlines()
    input_image,dx_phy =inputs[array_jobID].split(',')
    dx_phy = float(dx_phy)
#%% read image
imext = input_image.split('/')[-1].split('.')[1]
if imext=='npy':
    im = np.load(input_image)
elif imext=='tif' or imext=='tiff':
    im = io.imread(input_image)
else:
    raise ValueError('image data not understood')
imname = input_image.split('/')[-1].split('.')[0]
#transform input image to match with yantra data structure 
im = im[:,::-1]
im = im.T
#%% permeability simulation function
def get_permeability(direction):
    print("preparing to compute permeability in %s direction"%direction)
    global Fv,tol
    if direction == 'x':
        domain = yantra.Domain2D.use_img(im,dx=dx_phy,periodicity={'x':0,'y':0})
        bc_params={}
        dp_params = {"rhof":1000,
                     "visc":8.917e-7,
                     "p0":Pref}
        bc_params['left']= ['p',[Pref+DP/2,0]]
        bc_params['right']= ['p',[Pref-DP/2,0]]
        bc_params['top']=['u',[0,0]]
        bc_params['bottom']=['u',[0,0]]
    if direction == 'y':
        domain = yantra.Domain2D.use_img(im,dx=dx_phy,periodicity={'x':0,'y':0})
        bc_params={}
        dp_params = {"rhof":1000,
                     "visc":8.917e-7,
                     "p0":Pref}
        bc_params['top']= ['p',[Pref+DP/2,0]]
        bc_params['bottom']= ['p',[Pref-DP/2,0]]
        bc_params['left']=['u',[0,0,0]]
        bc_params['right']=['u',[0,0,0]]
    sp_params={'lattice':'D2Q9','collision_model':'trt','magic_para':3./16.,'tauref':1}
    ns = yantra.IncompressibleNavierStokes(domain,dp_params,bc_params,sp_params)
    print (ns._bc)
    print(ns.convfactors)
    ns.steady_state(tol=tol,check_iters=10,err_cal_method='l2norm',max_iters=max_iters)
    if direction == 'x': uavg=np.mean(ns.u[0])
    if direction == 'y': uavg=np.mean(ns.u[1])
    if direction == 'x': DX=np.mean(ns.x[-1]-ns.x[0])
    if direction == 'y': DX=np.mean(ns.y[-1]-ns.y[0])
    kphy = uavg*ns.rhof*8.917e-7/(DP/DX)   
    hydCond = kphy * 9.81/8.917e-7
    return kphy,hydCond, ns
#%%Get permeabilites 
timer = yantra.timer()
kx = 0
ky = 0
Kx = 0
Ky = 0
ns_x = None
ns_y = None
if directions =='all':
   kx,Kx, ns_x=  get_permeability('x')
   ky,Ky, ns_y=  get_permeability('y')
else:
    k,K,ns = get_permeability(directions)
    if directions == 'x':   kx, Kx,ns_x = k,K,ns
    if directions == 'y':   ky, Ky,ns_y = k,K,ns
#%%write outputs
k = ['int. permeability (m^2)',kx,ky]
k = [str(i) for i in k ]
K =['hydraullic conductivity (m/s)',Kx,Ky]
K = [str(i) for i in K ]
xyz = ['','x','y']
with open('%s.out'%imname,'w') as f:
    f.write(','.join(xyz) + '\n')
    f.write(','.join(k) + '\n')
    f.write(','.join(K))
if ns_x is not None:
    if  save_yantra_state:yantra.save(ns_x,'%s.x'%imname)
    if  export_vtk: yantra.export_to_vtk(ns_x,'%s.x'%imname,['nodetype','u','p'])
if ns_y is not None:
    if  save_yantra_state:yantra.save(ns_y,'%s.y'%imname)
    if  export_vtk: yantra.export_to_vtk(ns_y,'%s.y'%imname,['nodetype','u','p'])
        
print(timer)