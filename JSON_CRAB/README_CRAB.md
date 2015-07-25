# Running Code on Grid by using CRAB2
## Reference Twiki
-------------
https://twiki.cern.ch/twiki/bin/view/CMSPublic/WorkBookCRAB2Tutorial

## Config files 
----------------
### I) crab.cfg
To run a single job, we only need this file. It has following three blocks: 
#### i) [CMSSW]
We need to provide any two inputs out of these three:  
- Number of Events  
- Number of Jobs  
- Number of Events/Job  
A job will crash if number of jobs exceeds than 500. For example, if  
number of events = 5,000,000 
number of events/job = 8,000 
Then total number of jobs are 625 will result in a crash. 
- pset is the path to python file that runs our code.
- Input file location 
- Output file name 
 
#### ii) [USER]
- return_data
- copy_data 
- ui_working_dir asks path to input directory where jobs are created.  

#### iii) [GRID]
It has couple of standard definitions to use.

### II) multicrab.cfg
To run multiple jobs this file is needed along with above file. This file  
calls any block of crab.cfg and defines/redefines it.

## General Commands for crab/multicrab
--------------------------------
First set CMSSW environment in the working directory:  
```
cmsenv
```
Then source CRAB environment:  
```
source /afs/cern.ch/cms/ccs/wm/scripts/Crab/crab.sh
```
And create jobs in the directory defined in .cfg file:  
```
multicrab -create
```
Now submit the job:  
```
multicrab -submit -c /afs/cern.ch/work/a/asaddiqu/public/BH_MC
```
Status of all the jobs can be checked:
```
multicrab -status -c /afs/cern.ch/work/a/asaddiqu/public/BH_MC
```
Status of a single job can also be checked:  
```
crab -status -c /afs/cern.ch/work/a/asaddiqu/public/BH_MC/QCD_1000_inf
```
Output of a single data set can be retreived:
```
crab -getoutput -c /afs/cern.ch/work/a/asaddiqu/public/BH_MC/QCD_1000_inf
```
Note:  
The above output will be stored in directory:  
/afs/cern.ch/work/a/asaddiqu/public/BH_MC/QCD_1000_inf/res
