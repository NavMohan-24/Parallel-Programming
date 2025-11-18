# Add header paths
export CPLUS_INCLUDE_PATH="$CPLUS_INCLUDE_PATH:/opt/homebrew/opt/openblas/include:/opt/homebrew/opt/lapack/include:/opt/homebrew/include/eigen3:/opt/homebrew/include/eigen3/Eigen/:/opt/homebrew/include"

# Add library paths
export LIBRARY_PATH="$LIBRARY_PATH:/opt/homebrew/opt/openblas/lib:/opt/homebrew/opt/lapack/lib"

# Optionally also extend runtime library search path (for dynamic linking)
export LD_LIBRARY_PATH="$LD_LIBRARY_PATH:/opt/homebrew/opt/openblas/lib:/opt/homebrew/opt/lapack/lib"