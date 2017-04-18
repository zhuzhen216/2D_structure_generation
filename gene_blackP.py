
# coding: utf-8

# In[97]:

import numpy as np
import math
f_position = open('position.txt','w')
f_coord = open('coord.txt','w')
f_ANI = open('xyz.ANI','w')
f_neib = open('neib.txt','w')
# leave this to final
# input the number of cells to repeat
na1=int(input('Enter number of unit cells in a1 direction:'))
na2=int(input('Enter number of unit cells in a2 direction:'))
#
# pick the type of structures
#
sum_j=int(input('Please pick the sum of neighbore, 1 or 2 or 3:'))


# In[74]:

# comment this at the end
#na1 = 2
#na2 = 2
#sum_j = 2


# In[75]:

#
# atomic sites in primitive cell based on the primitive cell of black P
# this can be changed
#
ntau = 4


# In[76]:

# open output files
#f_neib = open('neib.txt','w')
#f_position = open('position.txt','w')
#f_coord = open('coord.txt','w')
#f_ANI = open('xyz.ANI','w')


# In[77]:

# atomic index
pho = 'P'


# In[78]:

#
# number of atoms in the supercell
#
nat = ntau*na1*na2


# In[79]:

lattice_vector = np.array([5.0*na1,3.3*na2,1.0])
#lattice_vector[0]=5.0*na1
#lattice_vector[1]=3.3*na2
#lattice_vector[2]=1.0


# In[80]:

#for i in range(3):
#    print(lattic_vector[i])


# In[81]:

#
# the fractional coordinate of black P
#
coor_unit = np.array([[0.0,0.16667,0.5,0.66667],[0.0,0.5,0.5,0.0]])
# or
# define a tuple list
#coor_unit = [(0,0),(0,0.5),(0.16667,0.5),(0.66667,0)]
#coor_unit[2][0]
#coor_unit2 = [1,1,-1,-1]
#list(zip(coor_unit[:][0],coor_unit2))


# In[82]:

# assign fraction coordinate to each atomic site in the supercell
coor = np.zeros((3,4*na1*na2))
for i in range(na1):
    for j in range(na2):
        atm_count=4*i*na2+j
        coor[0][atm_count]=coor_unit[0][0]/na1+i/na1
        coor[1][atm_count]=coor_unit[1][0]/na2+j/na2
        coor[0][atm_count+na2]=coor_unit[0][1]/na1+i/na1
        coor[1][atm_count+na2]=coor_unit[1][1]/na2+j/na2
        coor[0][atm_count+2*na2]=coor_unit[0][2]/na1+i/na1
        coor[1][atm_count+2*na2]=coor_unit[1][2]/na2+j/na2
        coor[0][atm_count+3*na2]=coor_unit[0][3]/na1+i/na1
        coor[1][atm_count+3*na2]=coor_unit[1][3]/na2+j/na2
#print(coor)


# In[83]:

coor_xyz = 1.0*np.zeros((3,4*na1*na2))
for i in range(nat):
    coor_xyz[0][i]=lattice_vector[0]*coor[0][i]
    coor_xyz[1][i]=lattice_vector[1]*coor[1][i]
    coor_xyz[2][i]=lattice_vector[2]*coor[2][i]
#coor_xyz


# In[84]:

# generate neighbor matrix
neib_dist = np.zeros((nat,nat))
# calculate the distance
i=0
while i < nat:
    j = i
    while j < nat:
        dis0=math.sqrt((coor_xyz[0][i]-coor_xyz[0][j])**2+(coor_xyz[1][i]-coor_xyz[1][j])**2+(coor_xyz[2][i]-coor_xyz[2][j])**2)
        dis1=math.sqrt((coor_xyz[0][i]+lattice_vector[0]-coor_xyz[0][j])**2+(coor_xyz[1][i]-coor_xyz[1][j])**2+(coor_xyz[2][i]-coor_xyz[2][j])**2)
        dis2=math.sqrt((coor_xyz[0][i]-coor_xyz[0][j])**2+(coor_xyz[1][i]+lattice_vector[1]-coor_xyz[1][j])**2+(coor_xyz[2][i]-coor_xyz[2][j])**2)
        dis3=math.sqrt((coor_xyz[0][i]+lattice_vector[0]-coor_xyz[0][j])**2+(coor_xyz[1][i]+lattice_vector[1]-coor_xyz[1][j])**2+(coor_xyz[2][i]-coor_xyz[2][j])**2)
        neib_dist[i][j]=min(dis0,dis1,dis2,dis3)
        neib_dist[j][i]=neib_dist[i][j]
        j = j + 1
    i = i + 1
#print(neib_dist)


# In[85]:

# each atom has three nearest neighbors
# neib_map would store the nearest neighbore info
# np.argsort would sort the array and return the index
neib_map = np.zeros((nat,3))
for i in range(nat):
    neib_map[i][0]=int(np.argsort(neib_dist[:][i])[1])
    neib_map[i][1]=int(np.argsort(neib_dist[:][i])[2])
    neib_map[i][2]=int(np.argsort(neib_dist[:][i])[3]) 
    output_value=(neib_map[i][0],neib_map[i][1],neib_map[i][2])
    f_neib.write(str(output_value))
#f_neib.close()
#neib_map


# In[86]:

#generate lattice gas model to the superstructure
ntot=1
for i in range(nat):
    ntot = ntot*2
ntot


# In[87]:

zat=np.zeros(nat)
zat


# In[95]:

###############################################
# this function will convert a number to binary
# in this program, I use binary to represent the
# up-and-down of the structures
# for 16 atoms, there are 2^16 structures
###############################################
def convert_bin(array_size,input_value):
    zat=np.zeros(array_size)
    k = input_value
    for j in range(array_size-1):
        zat[j] = k % 2
        k = k//2
    return zat
#convert_bin(16,0)


# In[98]:

#
#
# m, n, l are three neighboring atoms of j
#
#zat=np.zeros(nat)
for i in range(ntot):
    #zat[nat-1]=0
    k=i
#    for j in range(nat-1):
#        zat[j]=k%2
#        k=k//2
    zat=convert_bin(nat,i)
    add_zat = 0
    #print(zat)
    for j in range(nat):
        add_zat = add_zat + zat[j]
    if add_zat != nat/2:
        continue
    
    countj = 0
    for j in range(nat):
        z_neib=0
        m=int(neib_map[j][0])
        n=int(neib_map[j][1])
        l=int(neib_map[j][2])
        z_neib=int(zat[j]+zat[m]+zat[n]+zat[l])
        type(z_neib)
        if sum_j==2:
            if z_neib!=2: #==0 or z_neib==3 or z_neib==1 or z_neib==4:
                break
        elif sum_j==1 or sum_j==3:
            if z_neib==0 or z_neib==2 or z_neib==4:
                break
        countj = countj + 1
    #print(countj)
    if countj == nat:
        f_position.write(str(nat))
        f_position.write('\n\n')
        f_coord.write(str(nat))
        f_coord.write('\n\n')
        f_ANI.write(str(nat))
        f_ANI.write('\n\n')
        for j in range(nat):
            f_position.write(str(j+1)+ ' ' + str(zat[j])+'\n')
            f_ANI.write(pho+'\t'+str(coor_xyz[0][j])+'\t'+str(coor_xyz[1][j])+'\t'+str(zat[j])+'\n')
            f_position.write('\n')


# In[99]:

f_position.close()
f_ANI.close()
f_position.close()

# 2D-structure-generation
