# Bacterio
mol load pdb PROTEIN.pdb

while {[molinfo top get numreps] > 0} {mol delrep 0 top}
mol representation NewCartoon
mol color ColorID 8
mol selection all
mol material Opaque
mol addrep top


mol load pdb movie.pdb

graphics 0 line {-1000 5 -35} {1000 5 -35} width 2
graphics 0 line {-1000 5 -65} {1000 5 -65} width 2
graphics 0 line {-1000 -5 -35} {1000 -5 -35} width 2
graphics 0 line {-1000 -5 -65} {1000 -5 -65} width 2


while {[molinfo top get numreps] > 0} {mol delrep 0 top}
mol representation NewCartoon
mol color ColorID 1
mol selection all
mol material Opaque
mol addrep top


axes location off
display projection orthographic
display cuedensity 0
color Display Background white

user add key q {rotate x by 90}
user add key w {rotate x by -90}
puts "ahora anda lindo todo!!!"
