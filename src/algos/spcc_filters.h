/*
 * This file is part of Siril, an astronomy image processor.
 * Copyright (C) 2005-2011 Francois Meyer (dulle at free.fr)
 * Copyright (C) 2012-2023 team free-astro (see more in AUTHORS file)
 * Reference site is https://free-astro.org/index.php/Siril
 *
 * Siril is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * Siril is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with Siril. If not, see <http://www.gnu.org/licenses/>.
 */

#ifndef SPCC_FILTERS_H
#define SPCC_FILTERS_H

// Define the filters
spectral_intensity Johnson_B, Johnson_V, Optolong_Blue, Optolong_Green, Optolong_Red;

// Define the sensors
spectral_intensity Sony_IMX571M;

// Sensor transmittance data
float Johnson_UBVI_wl[] = {	360.f, 370.f, 380.f, 390.f, 400.f, 410.f, 420.f, 430.f, 440.f, 450.f, 460.f, 470.f,
									480.f, 490.f, 500.f, 510.f, 520.f, 530.f, 540.f, 550.f, 560.f, 570.f, 580.f, 590.f,
									600.f, 610.f, 620.f, 630.f, 640.f, 650.f, 660.f, 670.f, 680.f, 690.f, 700.f};

float Johnson_B_sr[] = { 		0.f, 0.03f, 0.134f, 0.567f, 0.92f, 0.978f, 1.f, 0.978f, 0.935f, 0.853f, 0.74f, 0.64f,
									0.536f, 0.424f, 0.325f, 0.235f, 0.15f, 0.095f, 0.043f, 0.09f, 0.f, 0.f, 0.f, 0.f, 0.f,
									0.f, 0.f, 0.f, 0.f, 0.f, 0.f, 0.f, 0.f, 0.f, 0.f};

float Johnson_V_sr[] = {		0.f, 0.f, 0.f, 0.f, 0.f, 0.f, 0.f, 0.f, 0.f, 0.f, 0.f, 0.f, 0.03f, 0.163f, 0.458f,
									0.78f, 0.967f, 1.f, 0.973f, 0.898f, 0.792f, 0.684f, 0.574f, 0.461f, 0.359f, 0.27f,
									0.197f, 0.135f, 0.081f, 0.045f, 0.025f, 0.017f, 0.013f, 0.009f, 0.f};


float Optolong_Blue_wl[] = {  380,382,384,386,388,390,392,394,396,398,400,402,404,406,408,410,412,414,416,418,
                                    420,422,424,426,428,430,432,434,436,438,440,442,444,446,448,450,452,454,456,458,
                                    460,462,464,466,468,470,472,474,476,478,480,482,484,486,488,490,492,494,496,498,
                                    500,502,504,506,508,510,512,514,516,518,520};

float Optolong_Blue_sr[] = {    0.487103,0.641832,0.702067,0.747306,0.792559,0.822237,0.84495,0.871089,0.889818,
                                    0.90034,0.900255,0.911452,0.932061,0.943272,0.942841,0.931497,0.932145,0.949223,
                                    0.958034,0.958862,0.96115,0.95695,0.950615,0.958319,0.964633,0.962691,0.963774,
                                    0.964245,0.963791,0.961815,0.961697,0.968066,0.969004,0.96647,0.968639,0.965631,
                                    0.961919,0.965998,0.968974,0.961068,0.960942,0.962811,0.959943,0.961922,0.9626,
                                    0.964021,0.9676,0.963836,0.960556,0.966053,0.964567,0.954384,0.956556,0.963683,
                                    0.966899,0.959314,0.952424,0.958815,0.95878,0.950698,0.9508,0.932266,0.852827,
                                    0.716573,0.544945,0.326982,0.138395,0.044799,0.014681,0.004195,0};

float Optolong_Green_wl[] = { 472,474,476,478,480,482,484,486,488,490,492,494,496,498,500,502,504,506,508,510,
                                    512,514,516,518,520,522,524,526,528,530,532,534,536,538,540,542,544,546,548,550,
                                    552,554,556,558,560,562,564,566,568,570};

float Optolong_Green_sr[] = { 0,0.00211,0.031464,0.098462,0.282547,0.428682,0.579431,0.70095,0.783954,0.824455,
                                    0.860767,0.897699,0.913057,0.922797,0.94429,0.958425,0.966152,0.966793,0.966125,
                                    0.968089,0.96738,0.968383,0.965945,0.967077,0.965016,0.964923,0.967084,0.964595,
                                    0.964797,0.966811,0.966968,0.965919,0.969215,0.971862,0.969491,0.966879,0.968857,
                                    0.965309,0.965128,0.962766,0.956522,0.950539,0.947976,0.9423,0.890753,0.672564,
                                    0.274385,0.067188,0.010192,0};

float Optolong_Red_wl[] = {   572,574,576,578,580,582,584,586,588,590,592,594,596,598,600,602,604,606,608,610,612,
                                    614,616,618,620,622,624,626,628,630,632,634,636,638,640,642,644,646,648,650,652,654,
                                    656,658,660,662,664,666,668,670,672,674,676,678,680,682,684,686,688,690,692,694,696,
                                    698,700};

float Optolong_Red_sr[] = {   0,0.02392,0.064858,0.119023,0.272804,0.372528,0.630103,0.817813,0.941956,0.971376,
                                    0.967548,0.972185,0.973656,0.972037,0.972571,0.971647,0.967951,0.970224,0.968007,
                                    0.96926,0.970465,0.966566,0.966465,0.967856,0.969271,0.969676,0.966155,0.967553,
                                    0.967322,0.968854,0.966806,0.967249,0.968454,0.968244,0.966644,0.966334,0.969328,
                                    0.969359,0.968264,0.966955,0.965487,0.968318,0.968753,0.970422,0.972236,0.973576,
                                    0.973862,0.973514,0.974338,0.975678,0.976927,0.974341,0.964257,0.96429,0.960805,
                                    0.956705,0.901116,0.672745,0.371317,0.191692,0.107868,0.025254,0.014725,0.006082,0};

// Sensor QE data

float Sony_IMX571_wl[] = {	420,430,440,450,460,470,480,490,500,510,520,530,540,550,560,570,580,590,600,610,620,
									630,640,650,660,670,680,690,700};

float Sony_IMX571_qe[] = {	0.789702,0.845261,0.875913,0.902740,0.914637,0.914116,0.897399,0.899673,0.897235,
									0.893329,0.883355,0.876612,0.863033,0.841178,0.825449,0.801004,0.777065,0.757080,
									0.723538,0.696523,0.676412,0.656132,0.625822,0.602794,0.577713,0.542536,0.530722,
									0.496716,0.471553};

// SPCC White Points

const cmsCIExyY Whitepoint_D58 = {0.32598, 0.33532, 1.0}; // Sun as white point, modelled as D58 black body
// TODO: could model this better based on a real solar spectrum

const cmsCIExyY Whitepoint_average_galaxy = {0.345702915, 0.358538597, 1.0}; // D50
// TODO: for testing, D50 is used. Needs replacing with a value computed from real galactic spectra.

#endif
