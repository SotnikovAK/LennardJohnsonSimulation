#pragma once
#pragma once
#include "Libraries.h"

const int fov = 90;
const double ratio = std::tan(fov / 2.0);
const int screenwidth = 1000;
const int screenheight = 1000;

const double alpha = 1;


const double MassParticale = 1.67 * pow(10, -26);
const double KT = 2.67 * pow(10, -21);


const double delta_t = 0.001;
const double Collisation = 1.4 * delta_t ;

const int max_n = 10000;

const double FinishForceRo = 6.25 * pow(alpha,2);

const int Temperature = 200;
const double bolzman = 1.38 * pow(10, -23);


const int endtime = 1000;

const double relaxTime = 400;


