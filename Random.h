#pragma once

#include <random>
std::seed_seq seed_{ time(NULL) };
std::mt19937 mt(seed_);

