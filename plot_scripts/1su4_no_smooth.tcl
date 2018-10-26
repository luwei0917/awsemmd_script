# Bacterio
#mol load pdb 2xov.pdb

mol load pdb movie.pdb
mol addfile {movie.psf} type {psf} first 0 last -1 step 1 waitfor 1 0
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



mol addrep top
mol color ColorID 1
mol selection resid 1 to 100
mol modrep 1 top

mol addrep top
mol color ColorID 3
mol selection resid 101 to 200
mol modrep 2 top

mol addrep top
mol color ColorID 4
mol selection resid 201 to 300
mol modrep 3 top


mol addrep top
mol color ColorID 7
mol selection resid 301 to 400
mol modrep 4 top

mol addrep top
mol color ColorID 10
mol selection resid 401 to 500
mol modrep 5 top

mol addrep top
mol color ColorID 0
mol selection resid 501 to 600
mol modrep 6 top

mol addrep top
mol color ColorID 11
mol selection resid 601 to 700
mol modrep 7 top
mol addrep top

mol color ColorID 12
mol selection resid 701 to 800
mol modrep 8 top
mol addrep top

mol color ColorID 13
mol selection resid 801 to 900
mol modrep 9 top
mol addrep top

mol color ColorID 14
mol selection resid 901 to 1000
mol modrep 10 top
mol addrep top

# mol color ColorID 15
# mol color ColorID 16

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
# pbc set {100 100 100} -all
# pbc box -center origin