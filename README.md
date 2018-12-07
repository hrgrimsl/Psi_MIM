# MIM
Clone this repo with:

git clone https://github.com/hrgrimsl/Psi_MIM.git

then from inside the Psi_MIM directory, use:

python setup.py build_ext --inplace

You will need a working Psi4 installation; binary and source should both work.  Instructions available here:

http://www.psicode.org/psi4manual/master/build_obtaining.html#how-to-obtain-psi4-start-with-find-the-code-quiz-end-in-top-level-psi4-dir

To test if it is working, call:

psi4 psicarbonyl out.dat

The optimization should run to completion.  

The files psinv_diamond and nvdiamond.cml will serve as a template for all of your MIM jobs.  Simply mirror the input structure shown in them file.  (For any line with a "/" in it, the left references the smaller fragments and the right of the "/" references the larger fragments, e.g. assigning methods as 'ccsd/scf' will perform ccsd on your smaller fragments and scf on your large ones.  Please note that at this time, asymmetric bases are not supported.)

The 'envir' variable in the Psi4 input file examples will be 'local' or 'cluster' depending on whether you're using an HPCC cluster or local environment.  I have no reason to think that this code will work on your cluster.  

An excellent explanation of the theory of MIM is available at:

https://pubs.acs.org/doi/pdf/10.1021/ct200033b

Please submit bug reports, questions, or miscellaneous vitriol to:

hrgrimsl@vt.edu
