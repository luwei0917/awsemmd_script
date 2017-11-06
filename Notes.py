# Useful codes
# awk '{print $NF}' all_wham.dat > e_total
# tr " " "\n"
# sed 1d
# sort -u -k 3
# sed -e 's/+T//'

# modify line
# awk '$5=-$5' data

# xargs - build and execute command lines from standard input
# sort    -g, --general-numeric-sort  compare according to general numerical value
# ls */wham.dat | sort -g | xargs cat > d2.dat
