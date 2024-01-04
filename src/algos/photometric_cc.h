#ifndef SRC_ALGOS_PHOTOMETRIC_CC_H_
#define SRC_ALGOS_PHOTOMETRIC_CC_H_

#include <stdio.h>
#include <glib.h>

#include "core/siril.h"
#include "core/proto.h"
#include "algos/PSF.h"
#include "algos/photometry.h"
#include "algos/astrometry_solver.h"

__inline __attribute__((always_inline)) int xisnanf(float x) { return x != x; }

typedef struct struct_coeff {
	float value;
	int channel;
} coeff;

typedef enum {
	CHANNEL_HIGHEST,
	CHANNEL_MIDDLE,
	CHANNEL_LOWEST
} normalization_channel;

struct photometric_cc_data {
	fits *fit;			// the image to process
	gboolean bg_auto;		// automatically select an area for bkg neutralization
	rectangle bg_area;		// the area for background if not bg_auto
	siril_catalogue *ref_stars;
	siril_cat_index catalog;		// catalog used for photometry
	limit_mag_mode mag_mode;	// automatically limit magnitude of the catalog
	double magnitude_arg;		// if not automatic, use this limit magnitude

	pcc_star *stars;		// the list of stars with BV index in the image
	int nb_stars;			// the number of stars in the array
	float fwhm;			// representative FWHM for stars
	gchar *datalink_path;	// to hold the datalink path for SPCC
	gboolean spcc;			// set if doing SPCC
	gboolean spcc_mono_sensor; // for SPCC
	int selected_sensor_m; // for SPCC
	int selected_sensor_rgb; // for SPCC
	int selected_filters; // for SPCC
	cmsCIExyYTRIPLE primaries; // used for SPCC source profile
	cmsCIExyY whitepoint; // used for SPCC source profile
};

int apply_photometric_color_correction(fits *fit, const float *kw, const coeff *bg, const float *mins, const float *maxs, int norm_channel);
int get_stats_coefficients(fits *fit, rectangle *area, coeff *bg, float *mins, float *maxs, int *norm_channel);
int photometric_cc(struct photometric_cc_data *args);
gpointer photometric_cc_standalone(gpointer p);
pcc_star *convert_siril_cat_to_pcc_stars(siril_catalogue *siril_cat, int *nbstars);
/*
#define TK_TABULATION_MAX 391
static const float tK[] = { 1000.0,1100.0,1200.0,1300.0,1400.0,1500.0,1600.0,1700.0,1800.0,1900.0,2000.0,2100.0,2200.0,2300.0,2400.0,2500.0,2600.0,2700.0,2800.0,2900.0,3000.0,3100.0,3200.0,3300.0,3400.0,3500.0,3600.0,3700.0,3800.0,3900.0,4000.0,4100.0,4200.0,4300.0,4400.0,4500.0,4600.0,4700.0,4800.0,4900.0,5000.0,5100.0,5200.0,5300.0,5400.0,5500.0,5600.0,5700.0,5800.0,5900.0,6000.0,6100.0,6200.0,6300.0,6400.0,6500.0,6600.0,6700.0,6800.0,6900.0,7000.0,7100.0,7200.0,7300.0,7400.0,7500.0,7600.0,7700.0,7800.0,7900.0,8000.0,8100.0,8200.0,8300.0,8400.0,8500.0,8600.0,8700.0,8800.0,8900.0,9000.0,9100.0,9200.0,9300.0,9400.0,9500.0,9600.0,9700.0,9800.0,9900.0,10000.0,10100.0,10200.0,10300.0,10400.0,10500.0,10600.0,10700.0,10800.0,10900.0,11000.0,11100.0,11200.0,11300.0,11400.0,11500.0,11600.0,11700.0,11800.0,11900.0,12000.0,12100.0,12200.0,12300.0,12400.0,12500.0,12600.0,12700.0,12800.0,12900.0,13000.0,13100.0,13200.0,13300.0,13400.0,13500.0,13600.0,13700.0,13800.0,13900.0,14000.0,14100.0,14200.0,14300.0,14400.0,14500.0,14600.0,14700.0,14800.0,14900.0,15000.0,15100.0,15200.0,15300.0,15400.0,15500.0,15600.0,15700.0,15800.0,15900.0,16000.0,16100.0,16200.0,16300.0,16400.0,16500.0,16600.0,16700.0,16800.0,16900.0,17000.0,17100.0,17200.0,17300.0,17400.0,17500.0,17600.0,17700.0,17800.0,17900.0,18000.0,18100.0,18200.0,18300.0,18400.0,18500.0,18600.0,18700.0,18800.0,18900.0,19000.0,19100.0,19200.0,19300.0,19400.0,19500.0,19600.0,19700.0,19800.0,19900.0,20000.0,20100.0,20200.0,20300.0,20400.0,20500.0,20600.0,20700.0,20800.0,20900.0,21000.0,21100.0,21200.0,21300.0,21400.0,21500.0,21600.0,21700.0,21800.0,21900.0,22000.0,22100.0,22200.0,22300.0,22400.0,22500.0,22600.0,22700.0,22800.0,22900.0,23000.0,23100.0,23200.0,23300.0,23400.0,23500.0,23600.0,23700.0,23800.0,23900.0,24000.0,24100.0,24200.0,24300.0,24400.0,24500.0,24600.0,24700.0,24800.0,24900.0,25000.0,25100.0,25200.0,25300.0,25400.0,25500.0,25600.0,25700.0,25800.0,25900.0,26000.0,26100.0,26200.0,26300.0,26400.0,26500.0,26600.0,26700.0,26800.0,26900.0,27000.0,27100.0,27200.0,27300.0,27400.0,27500.0,27600.0,27700.0,27800.0,27900.0,28000.0,28100.0,28200.0,28300.0,28400.0,28500.0,28600.0,28700.0,28800.0,28900.0,29000.0,29100.0,29200.0,29300.0,29400.0,29500.0,29600.0,29700.0,29800.0,29900.0,30000.0,30100.0,30200.0,30300.0,30400.0,30500.0,30600.0,30700.0,30800.0,30900.0,31000.0,31100.0,31200.0,31300.0,31400.0,31500.0,31600.0,31700.0,31800.0,31900.0,32000.0,32100.0,32200.0,32300.0,32400.0,32500.0,32600.0,32700.0,32800.0,32900.0,33000.0,33100.0,33200.0,33300.0,33400.0,33500.0,33600.0,33700.0,33800.0,33900.0,34000.0,34100.0,34200.0,34300.0,34400.0,34500.0,34600.0,34700.0,34800.0,34900.0,35000.0,35100.0,35200.0,35300.0,35400.0,35500.0,35600.0,35700.0,35800.0,35900.0,36000.0,36100.0,36200.0,36300.0,36400.0,36500.0,36600.0,36700.0,36800.0,36900.0,37000.0,37100.0,37200.0,37300.0,37400.0,37500.0,37600.0,37700.0,37800.0,37900.0,38000.0,38100.0,38200.0,38300.0,38400.0,38500.0,38600.0,38700.0,38800.0,38900.0,39000.0,39100.0,39200.0,39300.0,39400.0,39500.0,39600.0,39700.0,39800.0,39900.0,40000.0 };

static const float x_1931_2deg_jv[] = { 0.6499,0.6361,0.6226,0.6095,0.5966,0.5841,0.572,0.5601,0.5486,0.5375,0.5267,0.5162,0.5062,0.4965,0.4872,0.4782,0.4696,0.4614,0.4535,0.446,0.4388,0.432,0.4254,0.4192,0.4132,0.4075,0.4021,0.3969,0.3919,0.3872,0.3827,0.3784,0.3743,0.3704,0.3666,0.3631,0.3596,0.3563,0.3532,0.3502,0.3473,0.3446,0.3419,0.3394,0.3369,0.3346,0.3323,0.3302,0.3281,0.3261,0.3242,0.3223,0.3205,0.3188,0.3171,0.3155,0.314,0.3125,0.311,0.3097,0.3083,0.307,0.3058,0.3045,0.3034,0.3022,0.3011,0.3,0.299,0.298,0.297,0.2961,0.2952,0.2943,0.2934,0.2926,0.2917,0.291,0.2902,0.2894,0.2887,0.288,0.2873,0.2866,0.286,0.2853,0.2847,0.2841,0.2835,0.2829,0.2824,0.2818,0.2813,0.2807,0.2802,0.2797,0.2792,0.2788,0.2783,0.2778,0.2774,0.277,0.2765,0.2761,0.2757,0.2753,0.2749,0.2745,0.2742,0.2738,0.2734,0.2731,0.2727,0.2724,0.2721,0.2717,0.2714,0.2711,0.2708,0.2705,0.2702,0.2699,0.2696,0.2694,0.2691,0.2688,0.2686,0.2683,0.268,0.2678,0.2675,0.2673,0.2671,0.2668,0.2666,0.2664,0.2662,0.2659,0.2657,0.2655,0.2653,0.2651,0.2649,0.2647,0.2645,0.2643,0.2641,0.2639,0.2638,0.2636,0.2634,0.2632,0.2631,0.2629,0.2627,0.2626,0.2624,0.2622,0.2621,0.2619,0.2618,0.2616,0.2615,0.2613,0.2612,0.261,0.2609,0.2608,0.2606,0.2605,0.2604,0.2602,0.2601,0.26,0.2598,0.2597,0.2596,0.2595,0.2593,0.2592,0.2591,0.259,0.2589,0.2588,0.2587,0.2586,0.2584,0.2583,0.2582,0.2581,0.258,0.2579,0.2578,0.2577,0.2576,0.2575,0.2574,0.2573,0.2572,0.2571,0.2571,0.257,0.2569,0.2568,0.2567,0.2566,0.2565,0.2564,0.2564,0.2563,0.2562,0.2561,0.256,0.2559,0.2559,0.2558,0.2557,0.2556,0.2556,0.2555,0.2554,0.2553,0.2553,0.2552,0.2551,0.2551,0.255,0.2549,0.2548,0.2548,0.2547,0.2546,0.2546,0.2545,0.2544,0.2544,0.2543,0.2543,0.2542,0.2541,0.2541,0.254,0.254,0.2539,0.2538,0.2538,0.2537,0.2537,0.2536,0.2535,0.2535,0.2534,0.2534,0.2533,0.2533,0.2532,0.2532,0.2531,0.2531,0.253,0.253,0.2529,0.2529,0.2528,0.2528,0.2527,0.2527,0.2526,0.2526,0.2525,0.2525,0.2524,0.2524,0.2523,0.2523,0.2523,0.2522,0.2522,0.2521,0.2521,0.252,0.252,0.2519,0.2519,0.2519,0.2518,0.2518,0.2517,0.2517,0.2517,0.2516,0.2516,0.2515,0.2515,0.2515,0.2514,0.2514,0.2513,0.2513,0.2513,0.2512,0.2512,0.2512,0.2511,0.2511,0.2511,0.251,0.251,0.2509,0.2509,0.2509,0.2508,0.2508,0.2508,0.2507,0.2507,0.2507,0.2506,0.2506,0.2506,0.2505,0.2505,0.2505,0.2505,0.2504,0.2504,0.2504,0.2503,0.2503,0.2503,0.2502,0.2502,0.2502,0.2502,0.2501,0.2501,0.2501,0.25,0.25,0.25,0.25,0.2499,0.2499,0.2499,0.2498,0.2498,0.2498,0.2498,0.2497,0.2497,0.2497,0.2497,0.2496,0.2496,0.2496,0.2496,0.2495,0.2495,0.2495,0.2495,0.2494,0.2494,0.2494,0.2494,0.2493,0.2493,0.2493,0.2493,0.2492,0.2492,0.2492,0.2492,0.2491,0.2491,0.2491,0.2491,0.2491,0.249,0.249,0.249,0.249,0.2489,0.2489,0.2489,0.2489,0.2489,0.2488,0.2488,0.2488,0.2488,0.2487 };

static const float y_1931_2deg_jv[] = { 0.3474,0.3594,0.3703,0.3801,0.3887,0.3962,0.4025,0.4076,0.4118,0.415,0.4173,0.4188,0.4196,0.4198,0.4194,0.4186,0.4173,0.4158,0.4139,0.4118,0.4095,0.407,0.4044,0.4018,0.399,0.3962,0.3934,0.3905,0.3877,0.3849,0.382,0.3793,0.3765,0.3738,0.3711,0.3685,0.3659,0.3634,0.3609,0.3585,0.3561,0.3538,0.3516,0.3494,0.3472,0.3451,0.3431,0.3411,0.3392,0.3373,0.3355,0.3337,0.3319,0.3302,0.3286,0.327,0.3254,0.3238,0.3224,0.3209,0.3195,0.3181,0.3168,0.3154,0.3142,0.3129,0.3117,0.3105,0.3094,0.3082,0.3071,0.3061,0.305,0.304,0.303,0.302,0.3011,0.3001,0.2992,0.2983,0.2975,0.2966,0.2958,0.295,0.2942,0.2934,0.2927,0.2919,0.2912,0.2905,0.2898,0.2891,0.2884,0.2878,0.2871,0.2865,0.2859,0.2853,0.2847,0.2841,0.2836,0.283,0.2825,0.2819,0.2814,0.2809,0.2804,0.2799,0.2794,0.2789,0.2785,0.278,0.2776,0.2771,0.2767,0.2763,0.2758,0.2754,0.275,0.2746,0.2742,0.2738,0.2735,0.2731,0.2727,0.2724,0.272,0.2717,0.2713,0.271,0.2707,0.2703,0.27,0.2697,0.2694,0.2691,0.2688,0.2685,0.2682,0.2679,0.2676,0.2673,0.2671,0.2668,0.2665,0.2663,0.266,0.2657,0.2655,0.2652,0.265,0.2648,0.2645,0.2643,0.2641,0.2638,0.2636,0.2634,0.2632,0.2629,0.2627,0.2625,0.2623,0.2621,0.2619,0.2617,0.2615,0.2613,0.2611,0.2609,0.2607,0.2606,0.2604,0.2602,0.26,0.2598,0.2597,0.2595,0.2593,0.2592,0.259,0.2588,0.2587,0.2585,0.2584,0.2582,0.258,0.2579,0.2577,0.2576,0.2574,0.2573,0.2572,0.257,0.2569,0.2567,0.2566,0.2565,0.2563,0.2562,0.2561,0.2559,0.2558,0.2557,0.2555,0.2554,0.2553,0.2552,0.255,0.2549,0.2548,0.2547,0.2546,0.2545,0.2543,0.2542,0.2541,0.254,0.2539,0.2538,0.2537,0.2536,0.2535,0.2534,0.2533,0.2532,0.2531,0.253,0.2529,0.2528,0.2527,0.2526,0.2525,0.2524,0.2523,0.2522,0.2521,0.252,0.2519,0.2518,0.2517,0.2516,0.2516,0.2515,0.2514,0.2513,0.2512,0.2511,0.2511,0.251,0.2509,0.2508,0.2507,0.2507,0.2506,0.2505,0.2504,0.2503,0.2503,0.2502,0.2501,0.25,0.25,0.2499,0.2498,0.2497,0.2497,0.2496,0.2495,0.2495,0.2494,0.2493,0.2493,0.2492,0.2491,0.2491,0.249,0.2489,0.2489,0.2488,0.2487,0.2487,0.2486,0.2485,0.2485,0.2484,0.2484,0.2483,0.2482,0.2482,0.2481,0.2481,0.248,0.248,0.2479,0.2478,0.2478,0.2477,0.2477,0.2476,0.2476,0.2475,0.2474,0.2474,0.2473,0.2473,0.2472,0.2472,0.2471,0.2471,0.247,0.247,0.2469,0.2469,0.2468,0.2468,0.2467,0.2467,0.2466,0.2466,0.2465,0.2465,0.2464,0.2464,0.2463,0.2463,0.2463,0.2462,0.2462,0.2461,0.2461,0.246,0.246,0.2459,0.2459,0.2459,0.2458,0.2458,0.2457,0.2457,0.2456,0.2456,0.2456,0.2455,0.2455,0.2454,0.2454,0.2454,0.2453,0.2453,0.2452,0.2452,0.2452,0.2451,0.2451,0.245,0.245,0.245,0.2449,0.2449,0.2449,0.2448,0.2448,0.2447,0.2447,0.2447,0.2446,0.2446,0.2446,0.2445,0.2445,0.2445,0.2444,0.2444,0.2444,0.2443,0.2443,0.2443,0.2442,0.2442,0.2442,0.2441,0.2441,0.2441,0.244,0.244,0.244,0.2439,0.2439,0.2439,0.2438 };
*/
#endif /* SRC_GUI_PHOTOMETRIC_CC_H_ */

