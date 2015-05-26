## Folder hierarchy

I usually have a root directory that contains all simulation data, let's call it root. 

The folder hierarchy goes like root/ N / type / temp where :  

* N is the number of repeat unit is the polymer (actual number of particles is 2N)
* type is **cir** for circular or **lin** for linear systems
* temp is temperature in ?.?? format

*loop-create.sh* is a bash script that creates blank folders for a given N. It should be called from the $root directory. It takes a single integer parameter which is N.


## Submission scripts

The SGE submission scripts are named xsubmit-machine and they take 3 parameters in this order

* type (cir,lin)
* N (integer)
* temperature (float)

*loop-submit.sh* is a bash script to sumbit a bulk of jobs. It takes 2parameters which are

* submission script to be used (full path)
* N
