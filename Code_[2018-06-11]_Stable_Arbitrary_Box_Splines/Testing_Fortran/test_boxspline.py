import numpy as np

# It's happy hour time, we need to leave, this is out of
# hand. (Handle?, grey hair, tanned skin)


# Execute a blocking command with a subprocess, on completion provide
# the return code, stdout as string, and stderr as string. This should
# work across both Python2.7 and Python3.x as well as cross-platform.
#  INPUT:
#   command -- A list of strings or string (space separated) describing
#              a standard command as would be given to subprocess.Popen
# 
#  OUTPUT:
#   return_code -- Straight from the subprocess returncode
#   stdout      -- A list of strings that are the lines of stdout 
#                  produced by the execution of <command>
#   stderr      -- A list of strings that are the lines of stderr
#                  produced by the execution of <command>
def run(command, **popen_kwargs):
    import sys, subprocess
    # For Python 2.x the encoding is a string by default
    # For Python 3.6 and later the encoding can be given as an arguemnt
    if sys.version_info >= (3,6):
        popen_kwargs.update( dict(encoding="UTF-8") )
    proc = subprocess.Popen(command, stdout=subprocess.PIPE,
                            stderr=subprocess.PIPE, **popen_kwargs)
    stdout, stderr = proc.communicate()
    # Before python 3.6, the encoding had to be handled after collection
    if ((sys.version_info >= (3,)) and (sys.version_info[0] == 3)):
        if (type(stdout) != str): stdout = str(stdout, encoding="UTF-8")
        if (type(stderr) != str): stderr = str(stderr, encoding="UTF-8")
    # Remove windows specific characters and split by the new line
    if stdout: stdout = stdout.replace("\r","").split("\n")
    else:      stdout = ""
    if stderr: stderr = stderr.replace("\r","").split("\n")
    else:      stderr = ""
    # Return the exit code, standard out, and standard error
    return proc.returncode, stdout, stderr


# ========================================================
#                   Generate Test Cases     
# ========================================================

elements          = ["TenP", "Diag", "ZP"]
numbers_of_points = ["1", "2K", "4K"]
multiplicities    = [1, 2, 3, 4, 5]

for mult in multiplicities:
    for element in elements:
        for num_points in numbers_of_points:
            if (mult > 2) and (element == "ZP"): continue
            if (mult > 3) and (element == "Diag"): continue
            name = "_{num_points}_{element}_{mult}.csv".format(**{
                "num_points":num_points, "element":element, "mult":mult})
            # print(f"Starting {name}", end="  ...   ")
            # ====================================================================
            if element == "TenP":
                # Bi-linear / bi-quadratic / bi-cubic function
                dvecs = np.array([[1.,0.],
                                  [0.,1.]], order="F")
                mults = np.array([1 ,1 ], order="F", dtype=np.int32) * mult
            if element == "Diag":
                # Not sure if this is named, but it's another test.
                dvecs = np.array([[1.,0.,1.],
                                  [0.,1.,1.]], order="F")
                mults = np.array([1 ,1 ,1 ], order="F", dtype=np.int32) * mult
            if element == "ZP":
                # Zwart-Powell element
                dvecs = np.array([[1.,0., 1., 1.],
                                  [0.,1.,-1., 1.]], order="F")
                mults = np.array([1 ,1 ,1 ,1 ], order="F", dtype=np.int32) * mult
            # ====================================================================

            if num_points == "1":
                eval_pts = np.asarray(np.asarray([np.sum(dvecs*mult,axis=1)/2]).T, order='F', dtype=np.float64)
            if num_points == "2K":
                eval_pts = np.meshgrid(list(range(50)), list(range(50)))
                eval_pts = np.vstack((eval_pts[0].flatten(), eval_pts[1].flatten())).T
                eval_pts = np.asarray(eval_pts.T, order='F', dtype=np.float64)
            if num_points == "4K":
                padding = .05
                # Get the x min and max
                x_min_max = [sum(np.where(mults*dvecs[0,:] < 0, mults*dvecs[0,:], 0)), 
                             sum(np.where(mults*dvecs[0,:] > 0, mults*dvecs[0,:], 0))]
                x_min_max[0] -= (x_min_max[1] - x_min_max[0]) * padding
                x_min_max[1] += (x_min_max[1] - x_min_max[0]) * padding
                # Get the y min and max
                y_min_max = [sum(np.where(mults*dvecs[1,:] < 0, mults*dvecs[1,:], 0)), 
                             sum(np.where(mults*dvecs[1,:] > 0, mults*dvecs[1,:], 0))]
                y_min_max[0] -= (y_min_max[1] - y_min_max[0]) * padding
                y_min_max[1] += (y_min_max[1] - y_min_max[0]) * padding
                plot_points = int(4000**(0.5) + 0.5)

                x_vals = (np.linspace(*x_min_max, num=plot_points),)
                x_vals += (np.linspace(*y_min_max, num=plot_points),)
                x_vals = tuple(x.flatten() for x in np.meshgrid(*x_vals))
                x_vals = np.vstack(x_vals).T
                eval_pts = np.asarray(x_vals.T, order='F')

            # ====================================================================

            np.savetxt("dvecs"    +name, dvecs.T,    delimiter=",", header=",".join(map(str,dvecs.shape)),    comments="")
            np.savetxt("mults"    +name, mults,      delimiter=",", header=",".join(map(str,mults.shape)),    comments="", fmt="%d")
            np.savetxt("eval_pts" +name, eval_pts.T, delimiter=",", header=",".join(map(str,eval_pts.shape)), comments="")
            # print("done.")


# ========================================================
#                    Generate Test Code     
# ========================================================


program_body = """
PROGRAM TEST_BOXSPLINE
  ! Define local variables for testing the box-spline code.
  INTEGER :: DIM, NUM_DVECS, NUM_PTS
  REAL :: START, FINISH
  ! Define local variables for reading direction vectors,
  ! multiplicities, and evaluation points.
  INTEGER, PARAMETER :: R8=SELECTED_REAL_KIND(13)
  REAL(KIND=R8), DIMENSION(:,:), ALLOCATABLE :: UNIQUE_DVECS
  INTEGER,       DIMENSION(:),   ALLOCATABLE :: DVEC_MULTS
  REAL(KIND=R8), DIMENSION(:,:), ALLOCATABLE :: EVAL_PTS
  REAL(KIND=R8), DIMENSION(:),   ALLOCATABLE :: BOX_EVALS
  INTEGER :: ERROR
  ! Define an interface for the BOXSPLEV subroutine.
  INTERFACE 
    SUBROUTINE BOXSPLEV(UNIQUE_DVECS, DVEC_MULTS, EVAL_PTS, BOX_EVALS, ERROR)
      INTEGER, PARAMETER :: R8=SELECTED_REAL_KIND(13)
      REAL(KIND=R8), INTENT(IN),  DIMENSION(:,:) :: UNIQUE_DVECS
      INTEGER,       INTENT(IN),  DIMENSION(:)   :: DVEC_MULTS
      REAL(KIND=R8), INTENT(IN),  DIMENSION(:,:) :: EVAL_PTS
      REAL(KIND=R8), INTENT(OUT), DIMENSION(:)   :: BOX_EVALS
      INTEGER,       INTENT(OUT)                 :: ERROR
    END SUBROUTINE BOXSPLEV
  END INTERFACE
  ! Body code for running all tests.

  {read_and_call}

END PROGRAM TEST_BOXSPLINE
"""

call_code_to_format = """
  ! __________________________________________________________________
  ! 
  ! {element} element, {num_points} points, multiplicity {mult}
  ! 
  ! Read direction vectors
  OPEN(UNIT=10,FILE='dvecs{name}')
  READ(10,*) DIM, NUM_DVECS
  ALLOCATE(UNIQUE_DVECS(1:DIM,1:NUM_DVECS))
  READ(10,*) UNIQUE_DVECS
  CLOSE(10)
  ! Read multiplicities  
  OPEN(UNIT=10,FILE='mults{name}')
  READ(10,*) NUM_DVECS
  ALLOCATE(DVEC_MULTS(1:NUM_DVECS))
  READ(10,*) DVEC_MULTS
  CLOSE(10)
  ! Read evaluation points
  OPEN(UNIT=10,FILE='eval_pts{name}')
  READ(10,*) DIM,NUM_PTS
  ALLOCATE(EVAL_PTS(1:DIM,1:NUM_PTS))
  READ(10,*) EVAL_PTS
  CLOSE(10)
  ! Allocate evaluation result storage
  ALLOCATE(BOX_EVALS(1:NUM_PTS))
  ! Evaluate box-spline
  CALL CPU_TIME(START)
  CALL BOXSPLEV(UNIQUE_DVECS, DVEC_MULTS, EVAL_PTS, BOX_EVALS, ERROR)
  CALL CPU_TIME(FINISH)
  PRINT *, '{mult},{element},{num_points},', ERROR, ',', FINISH - START
  ! Free memory
  DEALLOCATE(UNIQUE_DVECS, DVEC_MULTS, EVAL_PTS, BOX_EVALS)
  ! ------------------------------------------------------------------
"""

all_calls = ""

# Generate all of the calls to the boxspline code that read various files.
for mult in multiplicities:
    for element in elements:
        for num_points in numbers_of_points:
            if (mult > 2) and (element == "ZP"): continue
            if (mult > 3) and (element == "Diag"): continue
            name = "_{0}_{1}_{2}.csv".format(
                num_points, element, mult)
            # Generate the calling code for the box spline program.
            all_calls += call_code_to_format.format(name=name,
                mult=mult, element=element, num_points=num_points)

with open("test_boxspline.f90", "w") as f:
    f.write(program_body.format(read_and_call=all_calls))



comp_comm = "{compiler} -lblas -llapack {opt} -o test real_precision.f90 test_boxspline.f90 {version}"
clean_comm = "rm *.o *.mod test"
run_comm = "./test"

compilers = [
#    "gfortran",
    "/home/f/ltw/bin/ftn95.sun",
    "ifort",
]

versions = [
    # "bs-dynamic.f90",
    # "bs-automatic.f90",
    "bs-manual.f90",]

optimizations = [
    "",
    "-O2",
    "-O3",
    "-Os",]

num_trials = 20
import sys

print("Compiler,Version,Optimization,Multiplicity,Element,Num Points,Error,Time")
for compiler in compilers:
    for v in versions:
        for opt in optimizations:
            command = comp_comm.format(compiler=compiler, version=v, opt=opt)
            # Clean up old files.
            code, out, err = run(clean_comm.split())
            # Compile new code.
            code, out, err = run(command.split())
            for test in range(num_trials):
                print(command, sys.stderr)
                # Run new code.
                code, out, err = run([run_comm])
                # Print output in CSV format.
                for row in out:
                    if len("".join(row).strip()) == 0: continue
                    values = [compiler, v, opt] + row.split(',')
                    values = [val.strip() for val in values]
                    print(",".join(values))

