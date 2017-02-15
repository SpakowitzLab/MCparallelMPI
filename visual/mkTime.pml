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
color br0,(name A1)
alter (name A1),vdw=0.1

show spheres, (name A2)
color r4,(name A2)
alter (name A2),vdw=0.1

show spheres, (name A3)
color r9,(name A3)
alter (name A3),vdw=0.1

show spheres, (name A4)
color brightorange,(name A4)
alter (name A4),vdw=0.1

show spheres, (name A5)
color brown,(name A5)
alter (name A5),vdw=0.1

show spheres, (name A6)
color hotpink,(name A6)
alter (name A6),vdw=0.1

show spheres, (name A7)
color deepsalmon,(name A7)
alter (name A7),vdw=0.1

show spheres, (name A8)
color green,(name A8)
alter (name A8),vdw=0.1

show spheres, (name A9)
color forest,(name A9)
alter (name A9),vdw=0.1



show (name P1)
color black, (name P1)
set stick_radius=0.05
show sticks, (name P1)

#run bar.py
#spectrumbar black

bg_color white
reset 
zoom center,15.5

mpng snap
#png snap060.png, dpi = 1000, width=1200, height=1200, ray=1
#set ray_trace_frames = 1
