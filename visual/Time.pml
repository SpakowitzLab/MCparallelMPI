delete *
reset

#movie.load snap*.pdb, snap
load time.pdb,snap 

set connect_mode,1
#rotate y,+70,state=0,camera=0
#rotate x,+20,state=0,camera=0
snap,mode=hist,gradient=bgr,minimum=0,maximum=1,nbins=40,sat=1.,value=1.


hide all
show spheres, (name A1)
color blue,(name A1)
alter (name A1),vdw=0.2
show sticks, (name A1)

show spheres, (name A2)
color red,(name A2)
alter (name A2),vdw=0.2
show sticks, (name A2)


show (name P1)
color black, (name P1)
set stick_radius=0.1
show sticks, (name P1)

#run bar.py
#spectrumbar black

bg_color white
reset 
zoom center,15.5

mpng snap
#png snap060.png, dpi = 1000, width=1200, height=1200, ray=1
#set ray_trace_frames = 1
