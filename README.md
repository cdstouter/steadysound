# Readme

## What is it

SteadySound is a small, feature-packed utility for adjusting the volume throughout audio files. It's basically an RMS meter, AGC (automatic gain control), and limiter combined, with a whole bunch of options. More documentation is coming soon, but you can see a list of the program options along with short descriptions by running `steadysound --help`.

## How to compile

Right now there's no build system, but it's very easy to compile. The development version of the following libraries are required:

* libsndfile
* Boost.Program_options

### Linux

On Ubuntu, install the dependencies like this:

    sudo apt install libsndfile1-dev libboost-program-options-dev


Other Linux distributions should have similar packages.

Compile the program:

    g++ *.cpp -Wall -lsndfile -lboost_program_options -o steadysound

### Windows

I use [MSYS2](http://www.msys2.org/) to compile on Windows. It's a great way to set up a MinGW environment, and it includes bash, the pacman package management system, and a whole bunch of ready to use packages. After installing it, you'll need launch the MSYS2 MSYS shell and install the required packages:

    pacman -S base-devel mingw-w64-i686-gcc mingw-w64-x86_64-gcc mingw-w64-i686-libsndfile mingw-w64-x86_64-libsndfile mingw-w64-i686-boost mingw-w64-x86_64-boost

Running that command will install everything you need to compile both 32-bit and 64-bit binaries. To build, launch either the MSYS2 MinGW 32-bit or 64-bit shell, depending on which architecture you want to build, and run this command:

    g++ *.cpp -Wall -lsndfile -lboost_program_options-mt -o steadysound.exe

Once it's compiled, you will probably have to copy some DLLs from \[MSYS2 installation directory\]/mingw\[32|64\]/bin/ to get it to run.

## License

This program is released under the MIT License. See the LICENSE file for more information. It also uses a few software libraries that are under different licenses:

* libsndfile - GNU LGPL
* boost\_program\_options - Boost Software License

For the details of these licenses, see the LICENSE-3RD-PARTY file.