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
Or in slurm machine you might encounter this error in conda enviroment
```
/data1/ddwang/software/anaconda3/envs/rid-drop/lib/libgrpc.so.20: undefined reference to `std::__throw_bad_array_new_length()@GLIBCXX_3.4.29'
/data1/ddwang/software/anaconda3/envs/rid-drop/lib/libtensorflow_framework.so.2: undefined reference to `std::__cxx11::basic_ostringstream<char, std::char_traits<char>, std::allocator<char> >::basic_ostringstream()@GLIBCXX_3.4.26'
/data1/ddwang/software/anaconda3/envs/rid-drop/lib/libtensorflow_cc.so.2: undefined reference to `std::__cxx11::basic_stringstream<char, std::char_traits<char>, std::allocator<char> >::basic_stringstream()@GLIBCXX_3.4.26'
```
Try to install libgcc in conda
```
conda install libgcc
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