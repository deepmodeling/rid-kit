# Installation
**You may encounter the error when compiling the plumed**
```
undefined reference to `std::__throw_bad_array_new_length()@glibcxx_3.4.29'
```
The reason for this error is that the tensorflow_cc library from conda requires GLIBCXX 3.4.29, but the highest version of GLICBCXX associated with the docker image is 3.4.28 (located in /usr/lib/x86_64-linux-gnu), so we need to upgrade the libstdc++6 library
```
DEBIAN_FRONTEND=noninteractive apt-get install -y software-properties-common
add-apt-repository -y ppa:ubuntu-toolchain-r/test
apt-get update
apt-get install -y libstdc++6
```

**You may encounter the error when compiling gromacs**
```
/opt/conda/lib/libcurl.so.4: no version information available (required by /usr/bin/cmake)
...
/opt/conda/lib/libtinfo.so.6: no version information available (required by /bin/bash)
```
The reason for this error is that there is version conflict of libcurl.so.4 and libtinfo.so.6 between the conda libraries and the libraries associated with ubuntu, we can relink the two libraries
```
rm ${CONDA_PREFIX}/lib/libtinfo.so.6
ln -s /usr/lib/x86_64-linux-gnu/libtinfo.so.6 ${CONDA_PREFIX}/lib/libtinfo.so.6
rm ${CONDA_PREFIX}/lib/libcurl.so.4
ln -s /usr/lib/x86_64-linux-gnu/libcurl.so.4 ${CONDA_PREFIX}/lib/libcurl.so.4
```