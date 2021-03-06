#!/bin/bash
path="~/sundials"
installPath="InstallSundials"
echo "Trying to download and install SUNDIALS/CVODE solver for you."
echo ""
read -p "Where to install?(Default: $path) :" path
path=${path:-~/sundials}
echo $path

echo "Installing SUNDIALS to $path"
echo "Installing examples to $path/example" 
echo ""

rm -rf $installPath
rm -rf sundials-master

echo "download from github: \n wget https://codeload.github.com/LLNL/sundials/zip/master"
wget -c https://codeload.github.com/LLNL/sundials/zip/master -O sundials-master.zip

echo "unzip sundials-master.zip"
unzip sundials-master.zip

mkdir $installPath
echo "cd $installPath"
cd $installPath
rm -rf builddir

echo "makedir builddir"
mkdir builddir
echo "cd builddir"
cd builddir/

#ccmake ../../sundials/
echo "start configure...."
cmake -DCMAKE_INSTALL_PREFIX=$path \
 -DEXAMPLES_INSTALL_PATH="$path/example" \
 -DBUILD_CVODE=ON \
 -DBUILD_CVODES=OFF \
 -DBUILD_ARKODE=OFF \
 -DBUILD_CVODES=OFF \
 -DBUILD_IDA=OFF \
 -DBUILD_IDAS=OFF \
 -DBUILD_KINSOL=OFF \
 -DOPENMP_ENABLE=OFF \
 ../../sundials-master

echo "make"
make

echo "make install"
make install
