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
parser.add_argument('--directions',type = str,default='all',help="x,y or all to get permeability in x,y or all directiodiff respectively")
parser.add_argument('--save_yantra_state',type=int,default=0,help="whether to save yantra state or not")
parser.add_argument('--export_vtk',type=int,default=1,help="whether to export velocity to vtk or not")
parser.add_argument('--max_iters',type=int,default=10000,help="max number of iteratiodiff allowed to reach steady state")
parser.add_argument('--cin',type=float,default=1,help="Inlet concentration")
parser.add_argument('--cout',type=float,default=0,help="Outlet concentration")    
parser.add_argument('--D0',type=float,default=0,help="Diffusion coefficent")    
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
cin = args.cin
cout = args.cout
D0 = args.D0
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
#tradiffform input image to match with yantra data structure 
im = im[:,::-1]
im = im.T
#%% permeability simulation function
def get_diffusivity(direction):
    print("preparing to compute diffusivity in %s direction"%direction)
    global Fv,tol
    if direction == 'x':
        domain = yantra.Domain2D.use_img(im,dx=dx_phy,periodicity={'x':0,'y':0})
        bc_params={}
        x0 = domain.x[0]
        xl = domain.x[-1]
        dcdl = (cout-cin)/(xl-x0)
        x,y = domain.meshgrid()
        c = cin + dcdl*(x-x0)
        dp_params = {"c":c,
                     "D":D0,
                    }
        bc_params['left']= ['c',cin]
        bc_params['right']= ['c',cout]
        bc_params['top']=['flux',0]
        bc_params['bottom']=['flux',0]
    if direction == 'y':
        domain = yantra.Domain2D.use_img(im,dx=dx_phy,periodicity={'x':0,'y':0})
        bc_params={}
        y0 = domain.y[0]
        yl = domain.y[-1]
        dcdl = (cout-cin)/(yl-y0)
        x,y = domain.meshgrid()
        c = cin + dcdl*(y-y0)
        dp_params = {"c":c,
                     "D":D0,
                    }
        bc_params['top']= ['c',cin]
        bc_params['bottom']= ['c',cout]
        bc_params['left']=['flux',0]
        bc_params['right']=['flux',0]
    sp_params={'collision_model':'trt','magic_para':1./4.,'tauref':1}
    diff = yantra.Diffusion(domain,dp_params,bc_params,sp_params)
    diff.steady_state(tol=tol,check_iters=10,err_cal_method='l2norm',max_iters=max_iters)
    if direction == 'x': fluxavg=np.mean(diff.flux[0])
    if direction == 'y': fluxavg=np.mean(diff.flux[1])
    if direction == 'x': DX=np.mean(diff.x[-1]-diff.x[0])
    if direction == 'y': DX=np.mean(diff.y[-1]-diff.y[0])
    DC = cin -cout
    De = fluxavg/(DC/DX)   
    return De, diff
#%%Get permeabilites 
timer = yantra.timer()
Dex = 0
Dey = 0
diff_x = None
diff_y = None
if directions =='all':
   Dex, diff_x=  get_diffusivity('x')
   Dey, diff_y=  get_diffusivity('y')
else:
    De,diff = get_diffusivity(directions)
    if directions == 'x':   Dex, diff_x = De,diff
    if directions == 'y':   Dey, diff_y = De,diff
#%%write outputs
De = ['Diffusivity (m^2/s)',Dex,Dey]
De = [str(i) for i in De ]
xyz = ['','x','y']
with open('%s.out'%imname,'w') as f:
    f.write(','.join(xyz) + '\n')
    f.write(','.join(De))
if diff_x is not None:
    if  save_yantra_state:yantra.save(diff_x,'%s.x'%imname)
    if  export_vtk: yantra.export_to_vtk(diff_x,'%s.x'%imname,['nodetype','flux','c'])
if diff_y is not None:
    if  save_yantra_state:yantra.save(diff_y,'%s.y'%imname)
    if  export_vtk: yantra.export_to_vtk(diff_y,'%s.y'%imname,['nodetype','flux','c'])
        
print(timer)