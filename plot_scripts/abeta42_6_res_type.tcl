# Bacterio
#mol load pdb 2xov.pdb

mol load pdb movie.pdb

graphics 0 line {-1000 45 15} {1000 45 15} width 2
graphics 0 line {-1000 45 -15} {1000 45 -15} width 2
graphics 0 line {-1000 -45 15} {1000 -45 15} width 2
graphics 0 line {-1000 -45 -15} {1000 -45 -15} width 2

# graphics 0 line {15 -1000 15} {15 1000 15} width 2
# graphics 0 line {-15 -1000 -15} {-15 1000 -15} width 2


while {[molinfo top get numreps] > 0} {mol delrep 0 top}


mol representation NewCartoon
mol material AOChalky
#mol modmaterial 0 0 AOChalky


axes location off
display projection orthographic
display cuedensity 0
color Display Background white

mol delrep 0 0
user add key q {rotate x by 90}
user add key w {rotate x by -90}



mol addrep 0
mol modcolor 0 0 ResType
mol modstyle 0 0 NewCartoon 0.300000 10.000000 4.100000 0
mol smoothrep 0 0 10


puts "ahora anda lindo todo!!!"

animate style Once
animate goto 0
display resetview
rotate x by 90.000000
rotate x by 90.000000
rotate x by 90.000000
scale by 1.200000

animate speed 1.000000
color Labels Bonds black

#animate speed 0.168831
#animate forward
#animate goto 450
#source "~/Downloads/take_picture.tcl"
#take_picture
#animate goto 200
#render Tachyon frame200.dat "'/Applications/VMD 1.9.2.app/Contents/vmd/tachyon_MACOSXX86'" -aasamples 12 %s -format TARGA -o frame200.tga -res 2000 2000
#animate goto 450
#render Tachyon frame450.dat "'/Applications/VMD 1.9.2.app/Contents/vmd/tachyon_MACOSXX86'" -aasamples 12 %s -format TARGA -o frame450.tga -res 2000 2000
#exit

# pbc set {272 272 272} -all
pbc set {100 100 100} -all
pbc box -center origin