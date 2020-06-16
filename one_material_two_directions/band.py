import math
import sys
import os
import meep as mp
import matplotlib.pyplot as plt
from meep import mpb

# fast-tun parameter for debugging
# N = 10
# num_bands = 3
# resolution = 128
# interp = int(sys.argv[1])

# real job
N = 10
num_bands = 180
resolution = 128
interp = int(sys.argv[1])

# honeycomb:
k_points = [mp.Vector3(-1/2,0)
           ,mp.Vector3( 1/2,0)]
#triangular
# k_points = [mp.Vector3(),          # Gamma
#             mp.Vector3(0.5,0),     # M
#             mp.Vector3(1/3, 1/3),  # K
#             mp.Vector3(0,0)] 
k_points = mp.interpolate(interp, k_points)

print("k_points",k_points)

eps_back_1 = 4.3**2
eps_cen_1 = 3.7**2

eps_back_2 = 3.9**2
eps_cen_2 = 4.5**2

# percentage = 0.1*int(sys.argv[3]) # degree of transistion from amorphous to crystaline

# argv[1]: radius of airhole
# argv[2]: radius of PCM

# triangular lattice with background InSb and middle PCM, air holes around
geometry = []

# geometry += [mp.Block(center=mp.Vector3(0,-N*3**0.5/4), material=mp.Medium(epsilon=eps_InSb_amor+(eps_InSb_crys-eps_InSb_amor)*0),
#                      size=mp.Vector3(1,N*3**0.5/2+2))]

# vertices = [mp.Vector3(0,-0.5/3**0.5),
#             mp.Vector3(1/2, -1/3**0.5),
#             mp.Vector3(1/2, N*3**0.5/2-1/3**0.5),
#             mp.Vector3(0, N*3**0.5/2-.5/3**0.5),
#             mp.Vector3(-1/2,N*3**0.5/2-1/3**0.5),
#             mp.Vector3(-1/2,-1/3**0.5)]

# geometry += [mp.Prism(vertices, height = mp.inf, material=mp.Medium(epsilon=eps_InSb_amor+(eps_InSb_crys-eps_InSb_amor)*1))]

# r_air = 0.001*int(sys.argv[1])
r_air = 1/6+1e-8
# r_cen = 0.001*int(sys.argv[2])
# r_cen = 1/3-r_air

vertices_hex = [mp.Vector3(     0,  2/(3**0.5)*r_air),
                mp.Vector3( r_air,  1/3**0.5*r_air),
                mp.Vector3( r_air, -1/3**0.5*r_air),
                mp.Vector3(     0, -2/3**0.5*r_air),
                mp.Vector3(-r_air, -1/3**0.5*r_air),
                mp.Vector3(-r_air,  1/3**0.5*r_air)]

vertices_l = [mp.Vector3(     0,  0),
              mp.Vector3( r_air,  1/3**0.5*r_air),
              mp.Vector3(     0,  2/(3**0.5)*r_air),
              mp.Vector3(-r_air,  1/3**0.5*r_air)]

vertices_r = [mp.Vector3(     0,  0),
              mp.Vector3( r_air,  -1/3**0.5*r_air),
              mp.Vector3(     0,  -2/(3**0.5)*r_air),
              mp.Vector3(-r_air,  -1/3**0.5*r_air)]

vertices_ul = [mp.Vector3(     0,  0),
               mp.Vector3(     0,  2/(3**0.5)*r_air),
               mp.Vector3( r_air,  1/3**0.5*r_air),
               mp.Vector3( r_air,  -1/3**0.5*r_air)]

vertices_lr = [mp.Vector3(     0,  0),
               mp.Vector3(     0,  -2/(3**0.5)*r_air),
               mp.Vector3(-r_air,  -1/3**0.5*r_air),
               mp.Vector3(-r_air,  1/3**0.5*r_air)]

vertices_ur = [mp.Vector3(     0,  0),
               mp.Vector3(     0,  -2/(3**0.5)*r_air),
               mp.Vector3( r_air,  -1/3**0.5*r_air),
               mp.Vector3( r_air,  1/3**0.5*r_air)]

vertices_ll = [mp.Vector3(     0,  0),
               mp.Vector3(     0,  2/(3**0.5)*r_air),
               mp.Vector3( -r_air, 1/3**0.5*r_air),
               mp.Vector3( -r_air, -1/3**0.5*r_air)]

for i in range(2):
    for j in range(N):
        center = (-1/2+i+1/2*(j%2),(-N+j)*3**0.5/2)
        geometry += [mp.Prism(vertices_hex, center = mp.Vector3(center[0],center[1]), height = mp.inf, material=mp.Medium(epsilon=eps_cen_1)),
                     mp.Prism(vertices_l, center = mp.Vector3(center[0],center[1]-3**.5/3+1/3**0.5*r_air), height = mp.inf, material=mp.Medium(epsilon=eps_back_1)),
                     mp.Prism(vertices_r, center = mp.Vector3(center[0],center[1]+3**.5/3-1/3**0.5*r_air), height = mp.inf, material=mp.Medium(epsilon=eps_back_1)),
                     mp.Prism(vertices_ul, center = mp.Vector3(center[0]-1/2+1/3**0.5*r_air*3**0.5/2,center[1]-3**.5/6+1/3**0.5*r_air/2), height = mp.inf, material=mp.Medium(epsilon=eps_back_1)),
                     mp.Prism(vertices_ur, center = mp.Vector3(center[0]-1/2+1/3**0.5*r_air*3**0.5/2,center[1]+3**.5/6-1/3**0.5*r_air/2), height = mp.inf, material=mp.Medium(epsilon=eps_back_1)),
                     mp.Prism(vertices_ll, center = mp.Vector3(center[0]+1/2-1/3**0.5*r_air*3**0.5/2,center[1]-3**.5/6+1/3**0.5*r_air/2), height = mp.inf, material=mp.Medium(epsilon=eps_back_1)),
                     mp.Prism(vertices_lr, center = mp.Vector3(center[0]+1/2-1/3**0.5*r_air*3**0.5/2,center[1]+3**.5/6-1/3**0.5*r_air/2), height = mp.inf, material=mp.Medium(epsilon=eps_back_1))]
for i in range(2):
    for j in range(N):
        center = (-1/2+i+1/2*(j%2),(j)*3**0.5/2)
        geometry += [mp.Prism(vertices_hex, center = mp.Vector3(center[0],center[1]), height = mp.inf, material=mp.Medium(epsilon=eps_cen_2)),
                     mp.Prism(vertices_l, center = mp.Vector3(center[0],center[1]-3**.5/3+1/3**0.5*r_air), height = mp.inf, material=mp.Medium(epsilon=eps_back_2)),
                     mp.Prism(vertices_r, center = mp.Vector3(center[0],center[1]+3**.5/3-1/3**0.5*r_air), height = mp.inf, material=mp.Medium(epsilon=eps_back_2)),
                     mp.Prism(vertices_ul, center = mp.Vector3(center[0]-1/2+1/3**0.5*r_air*3**0.5/2,center[1]-3**.5/6+1/3**0.5*r_air/2), height = mp.inf, material=mp.Medium(epsilon=eps_back_2)),
                     mp.Prism(vertices_ur, center = mp.Vector3(center[0]-1/2+1/3**0.5*r_air*3**0.5/2,center[1]+3**.5/6-1/3**0.5*r_air/2), height = mp.inf, material=mp.Medium(epsilon=eps_back_2)),
                     mp.Prism(vertices_ll, center = mp.Vector3(center[0]+1/2-1/3**0.5*r_air*3**0.5/2,center[1]-3**.5/6+1/3**0.5*r_air/2), height = mp.inf, material=mp.Medium(epsilon=eps_back_2)),
                     mp.Prism(vertices_lr, center = mp.Vector3(center[0]+1/2-1/3**0.5*r_air*3**0.5/2,center[1]+3**.5/6-1/3**0.5*r_air/2), height = mp.inf, material=mp.Medium(epsilon=eps_back_2))]


geometry_lattice = mp.Lattice(size=mp.Vector3(1, N*3**0.5),
                                 basis1=mp.Vector3(1,0),
                                 basis2=mp.Vector3(0,N*3**0.5))

ms = mpb.ModeSolver(num_bands=num_bands,
                    # target_freq=0.7,
                    tolerance = 1e-5, 
                    k_points=k_points,
                    geometry_lattice = geometry_lattice,
                    geometry=geometry,
                    resolution=resolution
                   )
# kxm = 0.5
# kym = 0.5
# interp = 4;

# kxm = 2
# kym = 2
# interp = 49
# N = interp+2
# step = kym/(interp+1)
# ms.k_points = []
# for i in range(N):
# 	ms.k_points += mp.interpolate(interp,[mp.Vector3(0,step*i),               # Gamma
#               mp.Vector3(kxm,step*i)])

# print("Triangle lattice of rods: TE bands(Hz)")
# ms.run_te(mpb.output_hfield_z)


ms.run_te()
# ms.run_tm()


md = mpb.MPBData(rectify=True, periods=1, resolution=resolution)
eps = ms.get_epsilon()
converted_eps = md.convert(eps)
plt.figure()
plt.imshow(converted_eps, interpolation='spline36', cmap='binary')
plt.colorbar()
# plt.axis('off')
try:
    this_str = './figures/'+'correct_n_cen_4100'+'n_diff_200'+'_reso_'+str(resolution)+'_interp_'+str(interp)
    os.mkdir(this_str)
except:
    print("error")
plt.savefig('./figures/'+'correct_n_cen_4100'+'n_diff_200'+'_reso_'+str(resolution)+'_interp_'+str(interp)+'/eps_'+str(resolution)+'.png', dpi=300)
# plt.imshow(eps.T, interpolation='spline36', cmap='binary')
# plt.axis('off')
# plt.show()
