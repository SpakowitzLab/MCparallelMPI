#  ----------------------------------------
#
#    Parses output file and makes pdb
#    Writen by Quin MacPherson
#    Started: 11/11/16
#
#   Example usage:
#   for i in `seq 1 31`; do python r2pdb.py r30v$i > snap0$i.pdb; done
#-------------------------------------------

nboundary=100
import math
import sys
if len(sys.argv) != 2:
      sys.exit('one arguement please')
fname= sys.argv[1]
#fname='full/r109v1'


# ---------------------
#    Read from file
# --------------------

nbeadsPerPoly=50
maxN=50000
lbox=25.0
totalPointsToDraw=60000

X=[];
Y=[];
Z=[];
AB=[];
METH=[];
n=0;
skip=1; count=0
with open(fname) as f:
   for line in f:
       temp=line.split()
       count=count+1;
       if count%skip != 0:
           continue
       #if count>maxN:
       #    continue
       x=float(temp[0])%lbox
       if x>5:
           continue
       X.append(float(temp[0])%lbox)  
       Y.append(float(temp[1])%lbox) 
       Z.append(float(temp[2])%lbox)
       if len(temp)==3:
           continue  
       AB.append(int(temp[3]))
       if len(temp)==4:
           continue  
       METH.append(int(temp[4]))

while len(X)<totalPointsToDraw:
    X.append(0.0)
    Y.append(0.0)
    Z.append(0.0)
    AB.append(7)

nbeads=len(X)

# -------------------------
#  Set atom types by sequence
# -------------------------
atomType=[]
for n in range(0,nbeads):
    #if X[n]<23 or X[n] >27:
    #    atomType.append('  A4')
    #    continue
    if AB[n]==1:
        atomType.append('  A2')
    elif AB[n]==7:
        atomType.append('  A7')
    else:
        atomType.append('  A1')

# -------------------------
#  Color in blocks
# -------------------------
#atomType=[]
#nAtomTypes=9
#for n in range(0,nbeads):
#    #if Z[n]<18 or Z[n] >22:
#    #    atomType.append('  A4')
#    #    continue
#    typ=int(round(nAtomTypes*float(n)/float(nbeads)+0.5))
#    typ=min(typ,nAtomTypes)
#    typ=max(typ,1)
#    atomType.append('  A%1d'%(typ))

# -------------------------
#  Set atom types by sequence
# -------------------------
#atomType=[]
#for n in range(0,nbeads):
#    #if X[n]<23 or X[n] >27:
#    #    atomType.append('  A4')
#    #    continue
#    if METH[n]==1:
#        atomType.append('  A2')
#    else:
#        atomType.append('  A1')


# -------------------
#  Write to pdb format
# ---------------------

chemName='Leviathan!'
if nbeads > 99999:
    raise 'too many beads for pdb format'

resname='WLC'
if len(resname) != 3:
    raise 'resname must be 3 letters'
print('HET    %s  A   1   %5d     Pseudo atom representation of DNA'
      %(resname,nbeads+nboundary))
print('HETNAM     %s %s'%(resname,chemName))
print('FORMUL  1   %s       '%(resname))

Ntot=0;
for n in range(0, nbeads):
    atomName=atomType[n]
    if len(atomName) != 4:
        raise 'Atom name needs to be 4 characters'
    x=X[n]
    y=Y[n]
    z=Z[n]
    print('ATOM  %5d %s %s          %8.3f%8.3f%8.3f  1.00  1.00           C'
            %(n+1,atomName,resname,x,y,z))
    #if n%nbeadsPerPoly != 0:
        #if(atomType[n-1] != '  A4' and atomType[n] !='  A4'):
        #    print('CONECT%5d%5d'%(n, n+1))
        #print('CONECT%5d%5d'%(n, n+1))
    Ntot=Ntot+1;


# -----------------
#  Draw confinement
# -----------------

# ----------------
#  Circle
# ----------------
#R=25
#for N in range(0,nboundary):
#     atomName='  A3'
#     x=25;
#     y=25+25*math.sin(2*3.141592*N/nboundary)
#     z=25+25*math.cos(2*3.141592*N/nboundary)
#     Ntot=Ntot+1;
#     print ('ATOM  %5d %s %s          %8.3f%8.3f%8.3f  1.00  1.00           C'
#           %(Ntot,atomName,resname,x,y,z))
#     if N != 0:
#         print('CONECT%5d%5d'%(Ntot,Ntot-1))
#print('CONECT%5d%5d'%(nbeads+1,Ntot))

# -------------------------
#  elongnated
#-------------------------
#R=10
#Length=20
#for N in range(0,nboundary):
#    atomName='  P1'
#    if N/nboundary < 0.5:
#        x=R+Length+R*math.sin(2*3.14159*N/nboundary)
#    else:
#        x=R+R*math.sin(2*3.14159*N/nboundary)
#    y=R+R*math.cos(2*3.141592*N/nboundary)
#    z=R
#    Ntot=Ntot+1;
#    print ('ATOM  %5d %s %s          %8.3f%8.3f%8.3f  1.00  1.00           C'
#          %(Ntot,atomName,resname,x,y,z))
#    if N != 0:
#        print('CONECT%5d%5d'%(Ntot,Ntot-1))
#print('CONECT%5d%5d'%(nbeads+1,Ntot))



