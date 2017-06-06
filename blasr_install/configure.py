from subprocess import call
import os


current_path = os.getcwd()
cmd = "export LD_LIBRARY_PATH=%s/blasr/libcpp/alignment:%s/blasr/libcpp/hdf:%s/blasr/libcpp/pbdata:%s/hdf5/hdf5-1.8.16-linux-centos6-x86_64-gcc447-shared/lib/:" %(current_path, current_path, current_path, current_path)
path = "%s/blasr/libcpp/alignment:%s/blasr/libcpp/hdf:%s/blasr/libcpp/pbdata:%s/hdf5/hdf5-1.8.16-linux-centos6-x86_64-gcc447-shared/lib/" %(current_path, current_path, current_path, current_path)
#call(cmd, shell = True)
os.environ["LD_LIBRARY_PATH"] += os.pathsep + path

cmd = "export"
call(cmd, shell = True)

cmd = "./blasr/blasr --version"
call(cmd, shell = True)
