# MIM
Clone this repo with:

git clone github.com/hrgrimsl/Psi_MIM/edit/master

You will need a working Psi4 installation; binary and source should both work.  Instructions available here:

http://www.psicode.org/psi4manual/master/build_obtaining.html#how-to-obtain-psi4-start-with-find-the-code-quiz-end-in-top-level-psi4-dir

To test if it is working, call:

psi4 psicarbonyl out.dat

The optimization should run to completion.  

The file 'psicarbonyl' will serve as a template for all of your MIM jobs.  Simply mirror the input structure shown in that file.  (For any line with a "/" in it, the left references the smaller fragments and the right of the "/" references the larger fragments, e.g. assigning methods as 'ccsd/scf' will perform ccsd on your smaller fragments and scf on your large ones.  Please note that at this time, asymmetric bases are not supported.)  You will also need a file described by 'name'.cml to describe bonds and what-not.  The geometry can be empty or incorrect, but you need the bond order shell.  Making new files can be done with Avogadro or whatever.

An excellent explanation of the theory of MIM is available at:

https://pubs.acs.org/doi/pdf/10.1021/ct200033b

Please submit bug reports, questions, or miscellaneous vitriol to:

hrgrimsl@vt.edu
