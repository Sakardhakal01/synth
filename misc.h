#pragma once

#ifndef MISC_H
#define MISC_H

#include <iostream>
#include <cmath>

using namespace std;

enum { SAW_DOWN, SAW_UP };

struct OSCIL
{
	double twopiovrsr;
	double curfreq;
	double curphase;
	double incr;
};

struct GTABLE
{
	double* table;
	unsigned long length;
};

GTABLE* new_gtable(unsigned long);
GTABLE* new_sine(unsigned long);
GTABLE* new_triangle(unsigned long, unsigned long);
GTABLE* new_square(unsigned long, unsigned long);
GTABLE* new_saw(unsigned long, unsigned long, int);
void gtable_free(GTABLE**);

struct OSCILT
{
	OSCIL osc;
	const GTABLE* gtable;
	double dtablen;
	double sizeovrsr;
};

OSCILT* new_oscilt(double, const GTABLE*, double);
double tick(OSCILT*, double);
double getphase(OSCILT*);

#endif
