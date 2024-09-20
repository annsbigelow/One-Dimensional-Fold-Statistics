# 1D fold stats


## Setup of configuration file
To compile the code it is necessary to create a common configuration file
called **config.mk** in the parent directory. Several templates are provided in
the **config** directory. To use, copy one of the templates into the parent
directory. From the main directory, on a Linux computer, type
```Shell
cp config/config.mk.linux ../config.mk
```
On a Mac with libraries installed via [MacPorts](http://www.macports.org), type
```Shell
cp config/config.mk.mac_mp ../config.mk
```
An additional configuration file for Windows will be added later. After setting up
this file, the code can be compiled by typing
```Shell
make
```
This will build several object files and executable files.
