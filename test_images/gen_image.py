import sys
sys.path.append('../../yantra-dev')
import yantra
import numpy as np
#%%2d image
dom = yantra.Domain2D((0,0),(80,80),dx=1,grid_type='midway')
dom.draw_circle((0,0),20)
dom.draw_circle((80,0),20)
dom.draw_circle((0,80),20)
dom.draw_circle((80,80),20)
dom.draw_circle((40,40),20)
yantra.export_to_vtk(dom,'2dimg',['nodetype'])
np.save('im2d',dom.nodetype)
#%%3d image
dom = yantra.Domain3D((0,0,0),(80,80,80),dx=1,grid_type='midway')
dom.draw_sphere((0,0,0),35)
dom.draw_sphere((80,0,0),35)
dom.draw_sphere((0,80,0),35)
dom.draw_sphere((0,0,80),35)
dom.draw_sphere((80,80,0),35)
dom.draw_sphere((80,0,80),35)
dom.draw_sphere((0,80,80),35)
dom.draw_sphere((80,80,80),35)
dom.draw_sphere((40,40,40),35)
yantra.export_to_vtk(dom,'3dimg',['nodetype'])
np.save('im3d',dom.nodetype)
print (1-np.mean(dom.nodetype==1))