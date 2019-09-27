# setup.sh
CC="gcc"   \
CXX="g++"   \
CFLAGS="-I/usr/include/"    \
LDFLAGS="-L/usr/lib64"   \
    python setup.py build_ext --inplace