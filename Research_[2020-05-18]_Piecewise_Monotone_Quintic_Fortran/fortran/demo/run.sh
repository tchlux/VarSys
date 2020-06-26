#!/bin/sh

SRC_DIR=".."
CUR_DIR="$(pwd)"
# Name the executable 'cli' for "command line interface".
EXEC="cli"

# Declare the variables for compilation.
F03=gfortran
OPTS="-O3"
MAIN="cli.f90"
EXEC_NAME="$CUR_DIR/$EXEC"
LIB="blas.f lapack.f"

# Move to the source code directory before compiling.
echo ""
echo "Moving to source code directory.."
cd "$SRC_DIR"

# Compile the command line interface executable.
echo ""
echo "Compiling the command line interface.."
$F03 $OPTS REAL_PRECISION.f90 EVAL_BSPLINE.f90 SPLINE.f90 MQSI.f90 $MAIN -o "$EXEC_NAME" $LIB

# Move back to the starting directory.
cd "$CUR_DIR"

# Run the executable, this will show the usage information.
echo ""
echo "Executing the CLI with no arguments to show usage information.."
echo "  ./$EXEC"
echo "-----------------------------------------------------------------"
./$EXEC
echo "-----------------------------------------------------------------"
echo "Using the command line interface to model 'demo.pts'.."
echo "  ./$EXEC demo.data demo.pts demo.out"
# Remove any existing output file.
rm -f "$CUR_DIR/demo.out"
# Run the CLI to get the output.
./$EXEC "demo.data" "demo.pts" "demo.out"
echo ""
echo "Done."
echo "  Check 'demo.out' for results of MQSI model."
