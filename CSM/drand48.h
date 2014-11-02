/*
 * The drand48 function - which exists on POSIX but not on Windows
 *
 * Extracted from mainRot.cpp during the 2014 reorganization
 * Created by Itay Zandbank
 */

#ifndef CSM_DRAND48_H  // Use the CSM_ prefix to make sure we're not confused with the standard C headers
#define CSM_DRAND48_H

double drand48(void);

#endif