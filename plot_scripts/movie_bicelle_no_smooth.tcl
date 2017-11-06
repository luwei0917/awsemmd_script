# Bacterio
#mol load pdb 2xov.pdb

mol load pdb movie.pdb

graphics 0 line {-1000 15 15} {1000 15 15} width 2
graphics 0 line {-1000 15 -15} {1000 15 -15} width 2
graphics 0 line {-1000 -15 15} {1000 -15 15} width 2
graphics 0 line {-1000 -15 -15} {1000 -15 -15} width 2


while {[molinfo top get numreps] > 0} {mol delrep 0 top}
mol representation NewCartoon
mol material AOChalky
#mol modmaterial 0 0 AOChalky
mol color ColorID 1
mol selection all
mol addrep top

axes location off
display projection orthographic
display cuedensity 0
color Display Background white

user add key q {rotate x by 90}
user add key w {rotate x by -90}

puts "ahora anda lindo todo!!!"

animate style Once
animate goto 0
display resetview
rotate x by 90.000000
rotate x by 90.000000
rotate x by 90.000000
scale by 1.200000

animate speed 1.000000
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
