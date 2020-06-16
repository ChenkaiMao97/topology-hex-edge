import matplotlib.pyplot as plt
import numpy as np
import sys

f_te = open("bandte.dat")
# f_tm = open("bandtm.dat")

# percentage = 0.1*int(sys.argv[3])

x_axis = np.linspace(-1/2,1/2,int(sys.argv[1])+2)

te_bands = [[]]*(int(sys.argv[1])+2)
f_te.readline()

counter_k = 0
counter_j = 0

for k in range(int(sys.argv[1])+2):
	x = f_te.readline()
	print(x)
	x = x[:-1] # get rid of \n
	te_bands[k]=[float(i) for i in x.split(",")[6:]]
# print(te_bands) 
plt.figure()

plt.plot(x_axis,te_bands,'r.', markersize=1)
plt.savefig('./te.png', dpi=300)

plt.show()


# tm_bands = [[]]*(int(sys.argv[3])+2)
# f_tm.readline()

# for k in range(int(sys.argv[3])+2):
# 		x = f_tm.readline()
# 		x = x[:-1] # get rid of \n
# 		tm_bands[k]=[float(i) for i in x.split(",")[6:]]

# # print(tm_bands) 
# plt.figure()
# plt.plot(x_axis,tm_bands,'b.', markersize=2)
# plt.show()
# plt.savefig('./figures/donuts/'+'r_pcm_'+str(0.001*int(sys.argv[1]))[:5]+'r_air_hole_'+str(0.001*int(sys.argv[2]))[:5]+'/tm.png',
# 			dpi=300)
# for x in f_tm:
#   	print(x)

