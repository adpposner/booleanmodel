/**
 * @file generatevars.h
 * @author Russell Posner (rposner@uchc.edu)
 * @brief Used to iterate over parameters using Sobol sequence
 * @version 0.1
 * @date 2019-06-08
 * 
 * @copyright Copyright (c) 2019
 * 
 */
#ifndef GENERATE_VARS_RP__
#define GENERATE_VARS_RP__

#include "f_interface.h"

void init_param_selector();

void end_param_selector();

modelparams getNextParams();
#endif