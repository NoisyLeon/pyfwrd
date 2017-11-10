#!/bin/csh

if ( $#argv != 5 )then
	  echo "############# hey, problem in running Minos_bran!!!!!"
	    echo $#argv
	      exit
endif

set cmplst = ( R T S )
set jcmp = $argv[1]
set moddir = $argv[2]
set modnm = $argv[3]
set w1 = $argv[4] #mHz, freq_min
set w2 = $argv[5] #freq_max
set cmp = ${cmplst[$jcmp]}

time /home/jixi7887/progs/jy/Mineos/Mineos-Linux64-1_0_2/minos_bran << EOF
$moddir/$modnm.txt
${modnm}_$cmp
none
1.0e-10 1
$jcmp
2 8000 $w1 $w2 0 0
EOF

