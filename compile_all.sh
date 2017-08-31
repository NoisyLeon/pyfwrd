cdir=`pwd`
echo $cdir
cd $cdir/tdisp96_src
./compile_tdisp96.sh

cd $cdir/tregn96_src
./compile_tregn96.sh

cd $cdir/tlegn96_src
./compile_tlegn96.sh

