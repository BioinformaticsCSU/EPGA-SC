Latest Version
==============
Please see the latest version of EPGA-SC:https://github.com/BioinformaticsCSU/EPGA-SC


License
=======

Copyright (C) 2017 Jianxin Wang(jxwang@mail.csu.edu.cn), Xingyu Liao(liaoxingyu@csu.edu.cn)

This program is free software; you can redistribute it and/or
modify it under the terms of the GNU General Public License
as published by the Free Software Foundation; either version 3
of the License, or (at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program; if not, see <http://www.gnu.org/licenses/>.

Jianxin Wang(jxwang@mail.csu.edu.cn), Xingyu Liao(liaoxingyu@csu.edu.cn)
School of Information Science and Engineering
Central South University
ChangSha
CHINA, 410083


Installation and operation of EPGA-SC 
==================================

### Dependencies

When running EPGA-SC from GitHub source the following tools are required:

* [jdk 1.8.0 or above]
* [GNU C++ 4.6.3 or above] 
* [perl 5.6.0 or above] 
* [python 2.7.14 or above]
 
### Install EPGA-SC

EPGA-SC is written Java and therefore will require a machine with jdk pre-installed.

Create a main directory (eg:EPGA-SC). Copy all source code to this directory.

Run command line: javac epga_sc.java 

### Run EPGA-SC.
	
    Run command line:  
	java -Xmx256G epga_sc -R 100 -c 13 -k 21 -K 31 -t 64 -q1 lib1.1.fastq -q2 lib1.2.fastq -o /home/epga-sc/finalResults/ [options] 
	
	[options]
	
	     -Xmx256G <This parameter is only used when running large data sets (For example: the data set size exceeds 5Gb)>
	     -R <The average length of read>
	     -c <The k-mer size used in coverage estimate(Default value: 13)>
	     -k <The k-mer size used in low depth reads assembly(Default value: 21)>
	     -K <The k-mer size used in normal depth reads assembly(Default value: 31)>
	     -t <The number of threads(Default value: 64)>
	     -i <The value of insertsizes>
	     -s <The standard devision of insertsizes(Default value: 0.1*insertsizes)>
	     -l <The threshold of low depth k-mer(Default value: 0.1)>
	     -r <The threshold of the number of low depth k-mer in a low depth reads(Default value: 0.01)>
	     -q1 <The first FASTQ file>
	     -q2 <The second FASTQ file>
         -o  <The path used to save final assemblies>		 
	
	If the system prompts "operation not permitted" ,we need to run the following commands to modify the permissions of EPGA-SC folder.
	
	chmod -R 777  EPGA-SC-master
	

### Output.
    
	(1)The final assemblies will be stored in the path specified by '-o'.
    
        /home/.../epga-sc/finalResults/Contigs.fa
		
        /home/.../epga-sc/finalResults/Scaffolds.fa
