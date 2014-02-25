/*************************************************************************\

  Copyright 2012 The University of North Carolina at Chapel Hill.
  All Rights Reserved.

  Permission to use, copy, modify and distribute this software and its
  documentation for educational, research and non-profit purposes, without
   fee, and without a written agreement is hereby granted, provided that the
  above copyright notice and the following three paragraphs appear in all
  copies.

  IN NO EVENT SHALL THE UNIVERSITY OF NORTH CAROLINA AT CHAPEL HILL BE
  LIABLE TO ANY PARTY FOR DIRECT, INDIRECT, SPECIAL, INCIDENTAL, OR
  CONSEQUENTIAL DAMAGES, INCLUDING LOST PROFITS, ARISING OUT OF THE
  USE OF THIS SOFTWARE AND ITS DOCUMENTATION, EVEN IF THE UNIVERSITY
  OF NORTH CAROLINA HAVE BEEN ADVISED OF THE POSSIBILITY OF SUCH
  DAMAGES.

  THE UNIVERSITY OF NORTH CAROLINA SPECIFICALLY DISCLAIM ANY
  WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF
  MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE.  THE SOFTWARE
  PROVIDED HEREUNDER IS ON AN "AS IS" BASIS, AND THE UNIVERSITY OF
  NORTH CAROLINA HAS NO OBLIGATIONS TO PROVIDE MAINTENANCE, SUPPORT,
  UPDATES, ENHANCEMENTS, OR MODIFICATIONS.

  The authors may be contacted via:

  US Mail:       GAMMA Research Group at UNC
                       Department of Computer Science
                       Sitterson Hall, CB #3175
                       University of N. Carolina
                       Chapel Hill, NC 27599-3175

  Phone:        (919)962-1749

  Authors:	Tang, Min tang_m@zju.edu.cn
		Manocha, Dinesh   dm@cs.unc.edu

  Publication:	a) Min Tang, Dinesh Manocha, Ruofeng Tong, MCCD: Multi-Core
                   collision detection between deformable models using
                   front-based decomposition, Graphical Models,
                   Vol. 72, No. 2, pp.7-23, 2010. 
		b) Min Tang, Dinesh Manocha, Ruofeng Tong. Multi-Core collision 
		   detection between deformable models. 2009 SIAM/ACM Joint 
		   Conference on Geometric and Physical Modeling, 2009, pp.355-360.
                     
  EMail:	geom@cs.unc.edu; tang_m@zju.edu.cn


\**************************************************************************/


This is a parallel continuous collision detection package optimized for multi-core CPUs.

------------------------------------------------------------------------------
1. Installation Instructions:
------------------------------------------------------------------------------

Version 1.0 of MCCD is developed by Visual Studio 2005, so currently it
 can only be used on Windows platform.
 

For rebuilding or running the demos:

Open mccd.sln, then build and run it.

To run the demos, please download data files from following links:
download http://gamma.cs.unc.edu/SR/benchmarks/cloth_ball.plys.zip (100 MB, 
104,915,287 bytes) and put the extracted ply files into d:\data\cloth_ball.plys\

There are two release versions to run.
Release-Serial: a serial implementation only uses a single thread.
Release-Parallel: a parallel implementation uses all the available cores.

After build, run mccd.exe for testing.

------------------------------------------------------------------------------
2. Who will need it?
------------------------------------------------------------------------------
I assume the reader is familiar with continuous collision detection.
MCCD is just designed for performing parallel continuous collision detection
between deformable model. It checks both inter-object & intra-object collisions.
MCCD (http://gamma.cs.unc.edu/PCD/) is optimized for multi-core platforms by 
using front-based decomposition. The idea has also be used for GPU-based CCD
acceleration (http://gamma.cs.unc.edu/CSTREAMS/).

Features:
Front-based decomposition

Please refer following papers for details:

a) Min Tang, Dinesh Manocha, Ruofeng Tong, MCCD: Multi-Core collision detection 
 between deformable models using front-based decomposition, Graphical Models,
 Vol. 72, No. 2, pp.7-23, 2010. 

b) Min Tang, Dinesh Manocha, Ruofeng Tong. Multi-Core collision detection between
 deformable models. 2009 SIAM/ACM Joint Conference on Geometric and Physical Modeling,
 2009, pp.355-360.

------------------------------------------------------------------------------
3. Bug report
------------------------------------------------------------------------------

We would be interested in knowing more about your application as well as any
bugs you may encounter in the collision detection library. You can
report them by sending e-mail to geom@cs.unc.edu or tang_m@zju.edu.cn.

------------------------------------------------------------------------------
4. Limitations
------------------------------------------------------------------------------

a. The current version supports up to 16-cores. It can be easily extended to more 
cores with some modificaitons (We have tested the algorithm on 24-core and 32-core 
machines.)

b. X64 version will be faster.

Release 1.0: 2012/3/22


Update at  2013/7/10:

1. Fixed the compile error for Visual Studio 2010.

2.  To set the number of threads, you need to change the environment variable OMP_NUM_THREADS, e.g.:

set OMP_NUM_THREADS[=num]

3. To run other benchmarks, you can download other benchmakrs from UNC Dynamic Scene Benchmarks (http://gamma.cs.unc.edu/DYNAMICB/).

For example,  download http://gamma.cs.unc.edu/SR/benchmarks/balls16_.plys.zip and put the extracted ply files into d:\data\balls16_.plys\, you will also need to modifiy the params.txt:

d:\data\balls16_.plys\balls16_
75
