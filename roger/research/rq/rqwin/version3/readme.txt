Quantile Regression function for Splus for Windows
--------------------------------------------------


The following files are distributed with this README file:

  rq.s:    	S-PLUS functions for computing quantile regressions. 
  rq.obj:  	Splus object
  rqhelp.txt:	Help file (ASCII)

Copy the files to a directory in your hard drive, say, 

  c:\programs

Then, start S-PLUS for Windows, and type

  source("c:\\programs\\rq.s")

At the beginning of each session, type

  dyn.load("c:\\programs\\rq.obj") 

and you will then be able to use rq. 

Notes: The rq function was modified from its UNIX implementation 
by commenting one of its lines. Do not try to use the UNIX version
when running S-PLUS for Windows.




Please, report any problems to Roger Koenker 
(roger@ysidro.econ.uiuc.edu). January 9, 1996



