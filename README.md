[Source](https://sites.google.com/a/lbl.gov/cs267-spring-2017/home/homework1?tmpl=/system/app/templates/print/&showPrintDialog=1 "Permalink to Assignment 1: Optimize Matrix Multiplication")

# Assignment 1: Optimize Matrix Multiplication

## Problem statement

Your task is to optimize [matrix multiplication][1] (matmul) code to run fast on a single processor core of [NERSC's Edison][2] cluster.

We consider a special case of matmul:

_C_ := _C_ \+ _A_*_B_

where _A_, _B_, and _C_ are _n_ x _n_ matrices. This can be performed using 2_n_3 floating point operations (_n_3 adds, _n_3 multiplies), as in the following pseudocode:
    
    
      for i = 1 to n
    
    
        for j = 1 to n
    
    
          for k = 1 to n
    
    
            C(i,j) = C(i,j) + A(i,k) * B(k,j)
    
    
          end
    
    
        end
    
    
      end

## Remote XSEDE/Moodle Students: Please Read

Dear Remote Students, we are thrilled to be a part of your parallel computing learning experience and to share these resources with you! To avoid confusion, please note that the assignment instructions, deadlines, and other assignment details posted here were designed for the local UC Berkeley students. You should check with your local instruction team about submission, deadlines, job-running details, etc. and utilize Moodle for questions. With that in mind, the problem statement, source code, and references should still help you get started (just beware of institution-specific instructions). Best of luck and we hope you enjoy the assignment!

## Instructions

* You will work in assigned teams (TBA) for this assignment.
* Your submission should be a gzipped tar archive, formatted (for Team 4) like: team04_hw1.tgz. It should contain: 
    * dgemm-blocked.c, a C-language source file containing your implementation of the routine: void square_dgemm(int, double*, double*, double*);
    * described in pseudocode above. We provide an example dgemm-blocked.c, below.
    * Makefile, only if you modified it. If you modified it, make sure it still correctly builds the provided benchmark.c, which we will use to grade your submission.
    * (e.g. for Team 4) team04_hw1.pdf, your write-up.
    * **Please do use these formats and naming conventions.** Not following these instructions leads to more busy work for the GSI's, which makes the GSI's sad...
* [This link][3] tells you how to use tar to make a .tgz file. 
* Submit your .tgz through bCourses.
* Your write-up should contain: 
    * the names of the people in your group (and each member's contribution),
    * the optimizations used or attempted,
    * the results of those optimizations,
    * the reason for any odd behavior (e.g., dips) in performance, and
    * how the performance changed when running your optimized code on a _different_ machine.
* For the last requirement, you may run your implementation on [another NERSC machine][4], on your laptop/cellphone, etc.
* Please carefully read the notes for implementation details. Stay tuned to [Piazza][5] ([signup][6]) for updates and clarifications, as well as discussion.
* If you are new to optimizing numerical codes, we recommend reading the papers in the references section.

### Notes:

* Your grade will mostly depend on two factors: 
    * performance sustained by your codes on Edison,
    * explanations of the performance features you observed (including what didn't work)
* There are other formulations of matmul (e.g., [Strassen][7]) that are mathematically equivalent, but perform asymptotically fewer computations - we will not grade submissions that do fewer computations than the two _n _cubed algorithm. This is actually an optional part of HW1.
* Your code must use [double-precision][8] to represent real numbers. A common reference for double-precision matrix multiplication is the [dgemm][9] (double-precision general matrix-matrix multiply) routine in the level-3 [BLAS][10]. We will compare your implementation with the tuned dgemm implementation in the vendor-provided BLAS library - on Edison (Cray XC30), we will compare with the Cray LibSci implementation of dgemm. Note that dgemm has a more general interface than square_dgemm - an optional part of HW1 encourages you to explore this richer tuning space.
* You may use [any compiler available][11]. We recommend starting with the GNU C compiler (gcc). If you use a compiler other than gcc, you will have to change the provided Makefile, which uses gcc-specific flags. Note that the default compilers, every time you open a new terminal, are PGI - you will have to type "module swap PrgEnv-pgi PrgEnv-gnu" to switch to, eg, GNU compilers. You can type "module list" to see which compiler wrapper you have loaded.
* Besides compiler intrinsic functions and built-ins, your code (dgemm-blocked.c) must only call into the C standard library.
* You may not use compiler flags that automatically detect dgemm kernels and replace them with BLAS calls, i.e. Intel's [-matmul][12] flag.
* You should try to use your compiler's automatic vectorizer before manually vectorizing. 
    * GNU C provides [many][13] extensions, which include intrinsics for vector (SIMD) instructions and data alignment. (Other compilers may have different interfaces.)
    * Ideally your compiler injects the appropriate intrinsics into your code automatically (eg, automatic vectorization and/or automatic data alignment). [gcc's auto-vectorizer][14] reports diagnostics that may help you identify if manual vectorization is required.
    * To manually vectorize, you must add compiler intrinsics to your code.
    * Consult your compiler's documentation.
* You may assume that A and B do not alias C; however, A and B may alias each other. It is semantically correct to qualify C (the last argument to square_dgemm) with the C99 restrict keyword. There is a lot online about restrict and pointer-aliasing - [this][15] is a good article to start with.
* The matrices are all stored in [column-major order][16], i.e. _Ci,j_ == C(i,j) == C[(i-1)+(j-1)*n], for i=1:n, where n is the number of rows in C. Note that we use 1-based indexing when using mathematical symbols and MATLAB index notation (parentheses), and 0-based indexing when using C index notation (square brackets).
* We will check correctness by the following componentwise error bound: |square_dgemm(n,A,B,0) - A*B| < eps*n*|A|*|B|. 
* where eps := 2-52 = 2.2 * 10-16 is the [machine epsilon][17].
* The [target processor][18] is a 12-core Intel Ivy Bridge running at 2.4 GHz, yielding 4 double-precision (ie, 64-bit) flops per pipeline * 2 pipelines * 2.4 GHz = 19.2 Gflops/s.

![][19]

### Optional:

These parts are not graded. You should be satisfied with your square_dgemm results and write-up before beginning an optional part.

* Implement Strassen matmul. Consider switching over to the three-nested-loops algorithm when the recursive subproblems are small enough.
* Support the dgemm interface (ie, rectangular matrices, transposing, scalar multiples).
* Try float (single-precision). This means you can use 8-way SIMD parallelism on Edison.
* Try complex numbers (single- and double-precision) - note that complex numbers are part of C99 and [supported in gcc][20]. [This forum thread][21] gives advice on vectorizing complex multiplication with the conventional approach - but note that there are [other algorithms][22] for this operation.
* Optimize your matmul for the case when the inputs are symmetric. Consider [conventional][23] and [packed][24] symmetric storage.
* Try optimizing matrix multiply for the [Cori Intel Xeon Phi Nodes][25] (the newest addition to NERSC).  The difference in architectures may make it very interesting for you.  You will need to write your own job scripts etc.  A good place to start is the [Getting Started page][26].  To add some incentive, **we are offering extra credit for this optional part this year. ** 

## Source Files

We provide two simple implementations for you to start with: a naive three-loop implementation similar to the pseudocode above, and a more cache-efficient blocked implementation.

The necessary files are in [cs267_hw1_2017.tar][27] posted for local students on bCourses. (Remote students, please check with your instructors and instruction site.) Included in the tar are the following: 

dgemm-naive.c

A naive implementation of matrix multiply using three nested loops,

dgemm-blocked.c

A simple blocked implementation of matrix multiply,

dgemm-blas.c

A wrapper for the vendor's optimized BLAS implementation of matrix multiply (default: Cray LibSci),

benchmark.c

The driver program that measures the runtime and verifies the correctness by comparing with the vendor's implementation,

Makefile

A simple makefile to build the executables,

job-blas, job-blocked, job-naive

_**Example**_ scripts*** **to run the executables on Edison compute nodes. For example, you can type "sbatch job-blas" to benchmark the BLAS version. 

The documentation for Edison's programming environment can be found below.

### Documentation:

You are also welcome to learn from the source code of state-of-art BLAS implementations such as [GotoBLAS][28] and [ATLAS][29]. However, you should not reuse those codes in your submission.

## Getting Running

Edison works with a batch submission system. Most of you have probably never worked with such a system, so welcome to the wonderful new world. The code you download should work out of the box. Below is an example of how you use batch submission. However, the NERSC documentation is great and is your friend. Once you work through this example, we strongly recommend you get more detail from NERSC. They have a nice example of [submitting your first job][30] that you should also work through. Details on how to modify the batch scripts can be found [here][31].

* SSH into edison and "wget" the above tarball of source code into some directory. "cd" into that directory and untar the source.  

`    ``    ``ssh -l  edison.nersc.gov  
`  
`        wget --no-check-certificate https://people.eecs.berkeley.edu/~qijing.huang/cs267_2017/hw1/cs267_hw1_2017.tar  
`  
`        tar xvf cs267_hw1_2017.tar`

* Compile the code with the makefile that we give you by typing "make" in the directory. This yields binaries that the computer can execute.
* Note that because we are writing sequential code for this assignment, the code will run on the login nodes (the node that you ssh into). You can run the naive code by typing "./benchmark-naive". Note that while this is feasible, running code on the login nodes is highly frowned upon and if you abuse it NERSC will get in touch with you. We also will not test your code on the login nodes, so you want to be benchmarking things on the compute nodes.
* To submit a job to be processed on the compute nodes, you need to submit a batch job script. These are easy to write, but we gave you one anyway *****, because we like you. To submit your job, you must choose a queue. There is more information online, but we'll test your code on the regular queue and you can debug in the debug queue. We chose the debug queue for you in the batch job script we give you, but use the regular queue to get real data for your final report. Submit the job-naive file with "sbatch job-blas". This will spit out something like "Submitted batch job JOBID" where JOBID is the id of your job.
* You can check on how your job is doing by querying the system with "squeue -u MYUSERNAME". That will spit out some status while you impatiently wait for your job to finish (it's so exciting!).
* When your job finishes running on a compute node, it will write an output (standard out) file and a error file (standard error). It's all very cool
* Note - lots of scientists all over are using the computers, so don't do bad things like run lots of code on the login nodes or query squeue constantly. This slows things down for other people and NERSC will get mad at you, which will make it hard for you to do your homework.
* Finally, please bring questions to office hours! We want you spending time on learning how to optimize code and understand the hardware.

### Important:

This documentation has been updated to reflect recent changes in the NERSC batching system. We have done our best to update all references, but beware artifacts may exist in outside sources and even these documents. Specifically, if you come across "qsub" instructions in NERSC-related resources, follow instead the "sbatch" instructions current on the NERSC site. If you are confused, please email us or stop by our office hours.

***** We emphasize these are example scripts because for these as well as all other assignment scripts we provide, you may need to adjust the number of requested nodes and cores and amount of time according to your needs (your allocation and the total class allocation is limited). To understand how you are charged, [READ THIS][32] alongside the given scripts. For testing (1) try running and debugging on your laptop first, (2) try running with the minimum resources you need for testing and debugging, (3) once your code is fully debugged, use the amount of resources you need to collect the final results for the assignment. This will become more important for later assignments, but it is good to get in the habit now.

## References

* Goto, K., and van de Geijn, R. A. 2008. Anatomy of High-Performance Matrix Multiplication, _ACM Transactions on Mathematical Software 34_, 3, Article 12. 
* (Note: explains the design decisions for the GotoBLAS dgemm implementation, which also apply to your code.)
* Chellappa, S., Franchetti, F., and PÃƒÂ¼schel, M. 2008. [How To Write Fast Numerical Code: A Small Introduction][33], _Lecture Notes in Computer Science 5235_, 196-259.
* (Note: how to write C code for modern compilers and memory hierarchies, so that it runs fast. Recommended reading, especially for newcomers to code optimization.)
* Bilmes, _et al._ [The PHiPAC (Portable High Performance ANSI C) Page for BLAS3 Compatible Fast Matrix Matrix Multiply][34].
* (Note: PHiPAC is a code-generating autotuner for matmul that started as a submission for this HW in a previous semester of CS267. Also see [ATLAS][29]; both are good examples if you are considering code generation strategies.)
* Lam, M. S., Rothberg, E. E, and Wolf, M. E. 1991. The Cache Performance and Optimization of Blocked Algorithms, _ASPLOS'91_, 63-74.
* (Note: clearly explains cache blocking, supported by with performance models.)
* Hints for HW1 and SIMD example pptx and pdf

 | 

[1]: http://www.google.com/url?q=http%3A%2F%2Fen.wikipedia.org%2Fwiki%2FMatrix_multiplication&sa=D&sntz=1&usg=AFQjCNF7P3ya4RJWT8a_cBRQhQhtHFLj8g
[2]: http://www.google.com/url?q=http%3A%2F%2Fwww.nersc.gov%2Fusers%2Fcomputational-systems%2Fedison%2F&sa=D&sntz=1&usg=AFQjCNEofn66_HYEJza5Ib-DUvj4CLA0MA
[3]: http://www.google.com/url?q=http%3A%2F%2Fkb.iu.edu%2Fdata%2Facfi.html&sa=D&sntz=1&usg=AFQjCNGG4xSgOGLH7zgAmIv_EPYedRbdJQ
[4]: http://www.google.com/url?q=http%3A%2F%2Fwww.nersc.gov%2Fusers%2Fcomputational-systems%2F&sa=D&sntz=1&usg=AFQjCNFpczBKfea_3UZRkBl9Gj-f39xLKw
[5]: https://www.google.com/url?q=https%3A%2F%2Fpiazza.com%2Fberkeley%2Fspring2017%2Fcs267%2Fhome&sa=D&sntz=1&usg=AFQjCNGgn8LIfjCRCK7tDZhYAf4fiOWKFg
[6]: https://www.google.com/url?q=https%3A%2F%2Fpiazza.com%2Fberkeley%2Fspring2016%2Fcs267&sa=D&sntz=1&usg=AFQjCNHBwSYeAdgW5aLhTv4fR_OI3DG8wA
[7]: http://www.google.com/url?q=http%3A%2F%2Fen.wikipedia.org%2Fwiki%2FStrassen_algorithm&sa=D&sntz=1&usg=AFQjCNHnZS7hPd_keDRJqilM_VuzvNLlOQ
[8]: http://www.google.com/url?q=http%3A%2F%2Fen.wikipedia.org%2Fwiki%2FDouble-precision_floating-point_format&sa=D&sntz=1&usg=AFQjCNHW3ExAuxbFdmsFQjZI0yxdmvJ2Ig
[9]: http://www.google.com/url?q=http%3A%2F%2Fwww.netlib.org%2Fblas%2Fdgemm.f&sa=D&sntz=1&usg=AFQjCNF55_TvoFch3H7sc14RC4Sspp1yag
[10]: http://www.google.com/url?q=http%3A%2F%2Fnetlib.org%2Fblas%2Findex.html&sa=D&sntz=1&usg=AFQjCNH8LcjkYZmN7pVINDPLy9qNWbl4zg
[11]: http://www.google.com/url?q=http%3A%2F%2Fwww.nersc.gov%2Fusers%2Fsoftware%2Fcompilers%2F&sa=D&sntz=1&usg=AFQjCNGLnwCqKY9K8aNPj-9_URHOTzkUZQ
[12]: http://www.google.com/url?q=http%3A%2F%2Fsoftware.intel.com%2Fsites%2Fproducts%2Fdocumentation%2Fhpc%2Fcomposerxe%2Fen-us%2F2011Update%2Ffortran%2Fwin%2Fcopts%2Fcommon_options%2Foption_opt_matmul.htm&sa=D&sntz=1&usg=AFQjCNFnGmJtz7xFqbFOK7G_gwbPtncOXQ
[13]: http://www.google.com/url?q=http%3A%2F%2Fgcc.gnu.org%2Fonlinedocs%2Fgcc%2FC-Extensions.html&sa=D&sntz=1&usg=AFQjCNFTuXG1DC_s0Bo6GJNgduA6AVKIyw
[14]: http://www.google.com/url?q=http%3A%2F%2Fgcc.gnu.org%2Fprojects%2Ftree-ssa%2Fvectorization.html&sa=D&sntz=1&usg=AFQjCNFfihLxGT75s7zqaJ0lrt-QM7GsUg
[15]: http://www.google.com/url?q=http%3A%2F%2Fcellperformance.beyond3d.com%2Farticles%2F2006%2F05%2Fdemystifying-the-restrict-keyword.html&sa=D&sntz=1&usg=AFQjCNGTYyQUZZowLWc4Uh-lGqPH3DUJJw
[16]: http://www.google.com/url?q=http%3A%2F%2Fen.wikipedia.org%2Fwiki%2FRow-major_order&sa=D&sntz=1&usg=AFQjCNEwA_zeB5feZYAGG-LiZRk4yYWKDQ
[17]: http://www.google.com/url?q=http%3A%2F%2Fen.wikipedia.org%2Fwiki%2FMachine_epsilon&sa=D&sntz=1&usg=AFQjCNHv1pzh_w99W-9rnsCMKx-x0C9bdQ
[18]: http://www.google.com/url?q=http%3A%2F%2Fwww.nersc.gov%2Fusers%2Fcomputational-systems%2Fedison%2Fconfiguration%2F&sa=D&sntz=1&usg=AFQjCNF2VCGxuccbHChss-ip_icRwPzUzw
[19]: https://ce8b96b5-a-46f4763e-s-sites.googlegroups.com/a/lbl.gov/cs267-spring-2017/home/homework1/mm_cs267.png?attachauth=ANoY7cqN4J3DcnJM92NINtcq6lQ4Rn7NnFZfu8AaokwQzktF40WQ_m290ZQqqzsOzSmCt9L_-X0Tvo-okuLBbYSAH293V-CP84In-04WWLCHlvI-TYKOWKaJnPsvgYkG599Qba2xHNY6EUb3aPdURyZ8buUv5rAdpFPQipowGx4UUUgh1hY_lNhoNV4o224PwY1AWvURij0kQ5Dp3JNrDPQR97c0VFrKn9FPqUrV_rmAiCdnBa12mq4%3D&attredirects=0
[20]: http://www.google.com/url?q=http%3A%2F%2Fgcc.gnu.org%2Fonlinedocs%2Fgcc%2FComplex.html&sa=D&sntz=1&usg=AFQjCNEuK_RzagaRPIhs5_oaLVzjxcyzUg
[21]: http://www.google.com/url?q=http%3A%2F%2Fstackoverflow.com%2Fquestions%2F3211346%2Fcomplex-mul-and-div-using-sse-instructions&sa=D&sntz=1&usg=AFQjCNFlVz3VkBTwoDXWqQN5RtVLzz1h7Q
[22]: http://www.google.com/url?q=http%3A%2F%2Fen.wikipedia.org%2Fwiki%2FMultiplication_algorithm%23Gauss.27s_complex_multiplication_algorithm&sa=D&sntz=1&usg=AFQjCNGxM1aGMfdzDlu_DuDBpQH0Es2GzA
[23]: http://www.google.com/url?q=http%3A%2F%2Fwww.netlib.org%2Flapack%2Flug%2Fnode122.html&sa=D&sntz=1&usg=AFQjCNGV9G9tzNpn3DqC3whVIPjPhZR4fA
[24]: http://www.google.com/url?q=http%3A%2F%2Fwww.netlib.org%2Flapack%2Flug%2Fnode123.html&sa=D&sntz=1&usg=AFQjCNGkU-RXjyZKWpnxmQZihGfjkB3_cg
[25]: http://www.nersc.gov/users/computational-systems/cori/cori-intel-xeon-phi-nodes/
[26]: http://www.nersc.gov/users/computational-systems/cori/getting-started/
[27]: https://www.google.com/url?q=https%3A%2F%2Fpeople.eecs.berkeley.edu%2F~qijing.huang%2Fcs267_2017%2Fhw1%2Fcs267_hw1_2017.tar&sa=D&sntz=1&usg=AFQjCNGUNrUH_0vYsG6ApO5-ToqnYv44Gw
[28]: http://www.google.com/url?q=http%3A%2F%2Fwww.tacc.utexas.edu%2Ftacc-projects%2Fgotoblas2&sa=D&sntz=1&usg=AFQjCNHVWHR0sUbhsWS83J0UwyW5Ut7o7w
[29]: http://www.google.com/url?q=http%3A%2F%2Fmath-atlas.sourceforge.net%2F&sa=D&sntz=1&usg=AFQjCNG0a6Gljk-YMhGkgnGpvKB3SxJcuQ
[30]: http://www.google.com/url?q=http%3A%2F%2Fwww.nersc.gov%2Fusers%2Fcomputational-systems%2Fedison%2Fgetting-started%2Fyour-first-program%2F&sa=D&sntz=1&usg=AFQjCNGoU6m0S6j2tnZ_VanwO_-SjvXtKg
[31]: http://www.google.com/url?q=http%3A%2F%2Fwww.nersc.gov%2Fusers%2Fcomputational-systems%2Fedison%2Frunning-jobs%2Fbatch-jobs%2F&sa=D&sntz=1&usg=AFQjCNGWdXnFyD1MqSOJnxEqs4py2CPJTw
[32]: http://www.google.com/url?q=http%3A%2F%2Fwww.nersc.gov%2Fusers%2Faccounts%2Fuser-accounts%2Fhow-usage-is-charged%2F&sa=D&sntz=1&usg=AFQjCNET3VLJ7VzGbeh-wEm7c4XEHirjxA
[33]: http://www.google.com/url?q=http%3A%2F%2Fspiral.ece.cmu.edu%3A8080%2Fpub-spiral%2Fabstract.jsp%3Fid%3D100&sa=D&sntz=1&usg=AFQjCNF78ztoHz6K4vBo4yTSfg4G5-wEcg
[34]: http://www.google.com/url?q=http%3A%2F%2Fwww.icsi.berkeley.edu%2F~bilmes%2Fphipac%2F&sa=D&sntz=1&usg=AFQjCNEHDeCYWhbZF_cpADrkJ5NatuPJaw

  
