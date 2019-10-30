set -e

cd build/
cmake ..
make
cd $OLDPWD

