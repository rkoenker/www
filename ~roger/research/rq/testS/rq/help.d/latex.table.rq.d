.if \n(.g \{\ 
.warn 0
.fp 5 B
.ftr C Bn\}
.BG
.FN  latex.table.rq 
.TL
Latex Table of Quantile Regression Results 
.DN
This function converts the output of table.rq to a latex table.
.CS
latex.table.rq(object, caption="caption goes here.", digits=3, file="a")
.RA
.AG object
an object of class table.rq
.OA
.AG caption
caption for the table
.AG digits
significant digits for table entries
.AG file
name of the latex file
.nr XX 0
.ft 1
.fi
.in .75i
.sp
.ti 0
.ta .75i 1i 1.25i 1.5i 1.75i 2i
VALUE:\t\c
.br
NULL
.SE
generates table in  specified file
.DT
Uses Frank Harrell's function latex.table.
.SA
table.rq, latex.table
.EX

.KW

.SH Documentation Class
function
.WR
