#!/bin/bash
#
bindir="/home/youyir/Cascadia/ProjectQ/src/prediction/"
bin="${bindir}/Q-test-yamauchi"
pargen="/home/youyir/Cascadia/ProjectQ/src/prediction/1DQ-Pargen.py"
home=`pwd`

echo $home
#=== compile code 
cd $bindir
echo "$bindir"
sh compile.sh

cd $home
echo $home

echo "generate Parameter files"
python $pargen -f h2000wet121 > dirlist

dirs=`cat dirlist`

for dir in ${dirs[@]}
do
    echo "$dir"
    cd $dir
    echo "$bin < PAR.${dir}.YT2016"
    $bin < PAR.${dir}.YT2016
    cd $home
done

rm  dirlist
