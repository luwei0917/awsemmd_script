# Bacterio
#mol load pdb 2xov.pdb

mol load pdb movie.pdb




while {[molinfo top get numreps] > 0} {mol delrep 0 top}

graphics 0 line {-1000 15 14.4} {1000 15 14.4} width 2
graphics 0 line {-1000 15 -14.4} {1000 15 -14.4} width 2
graphics 0 line {-1000 -15 14.4} {1000 -15 14.4} width 2
graphics 0 line {-1000 -15 -14.4} {1000 -15 -14.4} width 2

mol representation NewCartoon
mol material AOChalky
#mol modmaterial 0 0 AOChalky



mol addrep top
mol color ColorID 1
mol selection resid 1 to 25
mol modrep 1 top

mol addrep top
mol color ColorID 3
mol selection resid 26 to 56
mol modrep 2 top

mol addrep top
mol color ColorID 4
mol selection resid 57 to 79
mol modrep 3 top


mol addrep top
mol color ColorID 7
mol selection resid 80 to 106
mol modrep 4 top

mol addrep top
mol color ColorID 10
mol selection resid 107 to 133
mol modrep 5 top

mol addrep top
mol color ColorID 0
mol selection resid 134 to 158
mol modrep 6 top

mol addrep top
mol color ColorID 11
mol selection resid 159 to 181
mol modrep 7 top
mol addrep top

axes location off
display projection orthographic
display cuedensity 0
color Display Background white

mol delrep 0 0
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
animate speed 0.168831
#animate forward
#animate goto 450
#source "~/Downloads/take_picture.tcl"
#take_picture
#animate goto 200
#render Tachyon frame200.dat "'/Applications/VMD 1.9.2.app/Contents/vmd/tachyon_MACOSXX86'" -aasamples 12 %s -format TARGA -o frame200.tga -res 2000 2000
animate goto 1000
render Tachyon frame1000.dat "'/Applications/VMD 1.9.2.app/Contents/vmd/tachyon_MACOSXX86'" -aasamples 12 %s -format TARGA -o frame1000.tga -res 2000 2000
exit
