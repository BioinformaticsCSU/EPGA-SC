Latest Version
==============
Please see the latest version of EPGA-SC:https://github.com/BioinformaticsCSU/EPGA-SC


License
=======

Copyright (C) 2019 Jianxin Wang(jxwang@mail.csu.edu.cn), Xingyu Liao(liaoxingyu@csu.edu.cn)

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

Create a main directory (eg:EPGA-SC-master). Copy all source code to this directory.

Run command line: javac epga_sc.java 

### Run EPGA-SC.

[Usage]:

 java [options] epga_sc [main arguments]

[options]:

 * -Xmx300G : This parameter is only used when working with large data sets.

[Main arguments]:
 
 * -R  <int>: The average length of read (100).
 * -c  <int>: The k-mer size is used in coverage estimation (13bp).
 * -k  <int>: The k-mer size is used in low depth reads assembly (21bp).
 * -K  <int>: The k-mer size is used in normal depth reads assembly (31bp).
 * -t  <int>: The number of threads (64).
 * -i  <int>: The average value of insertsizes.
 * -s  <double>: The standard devision of insertsizes (0.1*insertsizes).
 * -l  <double>: The threshold of filtering low depth k-mer (0.1).
 * -r  <double>: The threshold of the number of low depth k-mers in a low depth read (0.1).
 * -q1 <string>: The file with left reads for paired-end reads.
 * -q2 <string>: The file with right reads for paired-end reads.
 * -o  <string>: The path used to save the final assemblies.
		 
    If the system prompts "operation not permitted", we need to run the following command to modify the permissions of the EPGA-SC folder.
	
    chmod -R 777  EPGA-SC-master
	

### Output.
    
	The final assemblies will be stored in the path specified by '-o'. (For example: -o /home/.../EPGA-SC-master/finalResults/ , the final assemblies will be stored in the following path.)
    
        /home/.../EPGA-SC-master/finalResults/Contigs.fa
		
        /home/.../EPGA-SC-master/finalResults/Scaffolds.fa
