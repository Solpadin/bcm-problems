### Installation


      $ sudo apt-get install libparmetis-dev
      $ sudo apt-get install gfortran
      $ sudo apt-get install g++
      $ sudo apt-get install openmpi-bin

### Run

       lhead@head:/srv/foamrun$ git clone https://github.com/Solpadin/bcm-problems.git tmp$(date +%Y-%m-%dT%H:%M:%S)
       lhead@head:/srv/foamrun$ cd tmp* && make
