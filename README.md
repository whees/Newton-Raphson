# Newton-Raphson

note to users:\
if you want to enter your own roots,\
you need to edit the variable 'roots'\
in the file roots.py;\
'roots' must be a list\
of all roots to be used;\
example: roots = [1,quat(0,-1,0,0),quat(0,1,0,0)]\
is a polynomial\
with roots at 1, -i, and i;\
 \

Also, for an odd degree polynomial,\
at least one root must be real;\
Any root with a non-zero imaginary part\
must be accompanied by its complex conjugate,\
so a polynomial with only imaginary roots must be of even degree;\
\



 
I recommend you keep the images small\
(around 200x200)\
to avoid taking forever\
 \
but for larger images\
lower the number of iterations\
in the 'newtons' function
