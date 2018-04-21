
Gurobi
Download gurobi here
http://www.gurobi.com/downloads/user/gurobi-optimizer

# from http://abelsiqueira.github.io/blog/installing-gurobi-7-on-linux/

tar -zxf gurobi7.0.2_linux64.tar.gz
gedit $HOME/.bashrc

input the following at the end of the file:

 export GUROBI_HOME="/FULL/PATH/gurobi702/linux64"
 export PATH="${PATH}:${GUROBI_HOME}/bin"
 export LD_LIBRARY_PATH="${LD_LIBRARY_PATH}:${GUROBI_HOME}/lib"
 export GRB_LICENSE_FILE="$HOME/gurobi.lic"


Notice the /FULL/PATH part. This is where you downloaded the package, but it has to be the full path there. You can use $HOME to reference your home folder. For instance, your Desktop is in $HOME/Desktop. My packages folder is $HOME/packages. Also notice the GRB_LICENSE_FILE variable. If you intend to put your license somewhere else, you’ll need to change that part.

After that is done you need to either (a) close and reopen the terminal or (b) input source $HOME/.bashrc.



Go to http://www.gurobi.com
get an account
become an academic user
get academic license key here

http://www.gurobi.com/downloads/user/licenses/free-academic

Run grbgetkey as instructed there.



# get ets_fiber_assigner (aka libETS)


pip install pybind11
git clone https://github.com/Subaru-PFS/ets_fiber_assigner.git
cd ets_fiber_assigner
python setup.py install

# get ics_cobraOps (aka collision detection code)
git clone https://github.com/Subaru-PFS/ics_cobraOps


