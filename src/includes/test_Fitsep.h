#pragma once
// Needed by Recipe being tested
#include "marquardt.h"
#include "sep_marquardt.h"
#include "fgauss_sep.h"

// Needed only by test routine
#ifndef _NR_RAN
#include "ran.h"
#define _NR_RAN
#endif

#ifndef _NR_GAMMA
#include "gamma.h"
#define _NR_GAMMA
#endif

#ifndef _NR_DEVIATES
#include "deviates.h"
#define _NR_DEVIATES
#endif

#ifndef _NR_FIT_EXAMPLES
#include "fit_examples.h"
#define _NR_FIT_EXAMPLES
#endif
