
set(MCCD_SRCS
  bvh_front.cpp
  DeformBVH.cpp
  orphan.cpp
  timing.cpp
  deformingNow.cpp
  DeformBVH-flat.cpp
  DeformModel.cpp
  ccd_vf.cpp
  intersect.cpp
  UpperBVH-flat.cpp
  Decompose.cpp
  UpperBVH.cpp
  memory.cpp
)

add_library(mccd ${MCCD_SRCS})

