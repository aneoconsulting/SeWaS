#ifndef __CONSTANTS_H__
#define __CONSTANTS_H__

/* Spatial dimensions */
enum
{
  _X,
  _Y,
  _Z,
  _DIM
};

/* Stress field components */
enum
{
  _XX,
  _YY,
  _ZZ,
  _XY,
  _XZ,
  _YZ,
  _NB_STRESS_FIELD_COMPONENTS
};

/* Halo locations */
enum
{
  _LEFT,
  _RIGHT,
  _BACKWARD,
  _FORWARD,
  _BOTTOM,
  _TOP,
  _NB_LOCATIONS
};

/* Task types */
enum
{
  COMPUTE_VELOCITY,
  EXTRACT_VELOCITY_HALO,
  UPDATE_VELOCITY,
  DUMP_VELOCITY,
  DISPLAY_VELOCITY,
  COMPUTE_STRESS,
  EXTRACT_STRESS_HALO,
  UPDATE_STRESS,
  DUMP_STRESS,
  DISPLAY_STRESS,
  INITIALIZE_FIELDS,
  NB_TASK_TYPES
};

#endif
