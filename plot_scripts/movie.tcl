mol load pdb movie.pdb

while {[molinfo top get numreps] > 0} {mol delrep 0 top}

mol representation NewCartoon
mol material AOChalky
#mol modmaterial 0 0 AOChalky


while {[molinfo top get numreps] > 0} {mol delrep 0 top}
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

translate by 0.0 1.0 0.0
scale by 0.833000
scale by 0.9
animate speed 1.0
mol smoothrep 0 0 10

animate goto 99
display resetview
animate goto 1

#animate forward
#animate goto 450
#source "~/Downloads/take_picture.tcl"
#take_picture
#animate goto 200
#render Tachyon frame200.dat "'/Applications/VMD 1.9.2.app/Contents/vmd/tachyon_MACOSXX86'" -aasamples 12 %s -format TARGA -o frame200.tga -res 2000 2000
#animate goto 450
#render Tachyon frame450.dat "'/Applications/VMD 1.9.2.app/Contents/vmd/tachyon_MACOSXX86'" -aasamples 12 %s -format TARGA -o frame450.tga -res 2000 2000
#exit
