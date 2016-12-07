/*
 * sampleVonMisesMex.c - This function generates samples from a unimodal 
 * 		von Mises distribution. The sampling scheme used by this 
 * 		function was adopted from [1].
 *
 * The calling syntax is:
 *
 *		angles = sampleVonMises(mu, kappa, nSamples)
 * 
 * References:
 * 		[1] L. Barabesi (1995): "Generating von Mises Variates by the 
 * 			Ratio-of-Uniforms Method"
 *
 * Author:
 *      Copyright (c) 2016
 *      Christopher Schymura (christopher.schymura@rub.de)
 * 		Cognitive Signal Processing Group
 * 		Ruhr-Universitaet Bochum
 * 		Universitaetsstr. 150, 44801 Bochum
 * 
 * License information:
 * 		This program is free software: you can redistribute it and/or 
 *		modify it under the terms of the GNU General Public License as 
 *		published by the Free Software Foundation, either version 3 of 
 * 		the License, or (at your option) any later version.
 *
 * 		This material is distributed in the hope that it will be useful,
 *		but WITHOUT ANY WARRANTY; without even the implied warranty of
 *		MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
 *		GNU General Public License for more details.
 *
 *		You should have received a copy of the GNU General Public License
 *		along with this program. If not, see <http://www.gnu.org/licenses/>
 * 
 * This is a MEX file for MATLAB.
*/

#include "mex.h"
#include <stdlib.h>
#include <time.h>
#include <math.h>

/* Random number generation for the von Mises distribution */
void sample(double mu, double kappa, double *angles, mwSize numSamples)
{
	/* Initialize variables */
	mwSize i;					/* Counter variable */
	double samplingParameter;	/* Sampling parameter (see [1]) */
	double angle;				/* Angular value during sampling */
	double randomSample1;		/* Random values (uniform) */
	double randomSample2;
	
	/* Initialize random number generator */
	srand( time(NULL) + clock() + random() );
	
	/* Compute sampling parameter */
	if ( kappa > 1.3 ) {
		samplingParameter = 1 / ( sqrt(kappa) );
	}
	else {
		samplingParameter = M_PI * exp( -kappa );
	}
	
	for ( i = 0; i < numSamples; i++ ) {
		while ( true ) {
			/* Generate two random samples from the uniform distribution */
			randomSample1 = (double) rand() / (double) RAND_MAX;
			randomSample2 = (double) rand() / (double) RAND_MAX;
			
			/* Compute angular value */
			angle = samplingParameter * ( 2 * randomSample2 - 1 ) / randomSample1;
			
			/* Check if angle is in valid range - if not, start over */
			if ( fabs(angle) > M_PI ) {
				continue;
			}
			
			/* Check first stopping condition */
			if ( kappa * pow(angle, 2) < 4 - 4 * randomSample1 ) {
				break;
			}
			
			/* Check second stopping condition */
			if ( kappa * cos(angle) < 2 * log(randomSample1) + kappa ) {
				continue;
			}
			else {
				break;
			}
		}
		/* Assign sample to output vector */
		angles[i] = atan2( sin(angle + mu), cos(angle + mu) );		
	}	
}

/* Main function */
void mexFunction( int nlhs, mxArray *plhs[],
                  int nrhs, const mxArray *prhs[] )
{
	/* Initialize variables */
	double mu;					/* Circular mean */
	double kappa;				/* Concentration parameter */
	mwSize numSamples; 			/* Number of samples to be generated */
	double *angles;				/* Vector of generated samples */
	
	/* Check number of inputs */
	if ( nrhs != 2 && nrhs != 3 ) {
		mexErrMsgIdAndTxt( "VonMisesMixture:sampleVonMises:nrhs",
						  "Wrong number of input arguments." );
	}
	
	/* Check number of outputs */
	if ( nlhs != 1 ) {
		mexErrMsgIdAndTxt( "VonMisesMixture:sampleVonMises:nlhs",
						  "One output required." );
	}
	
	/* Check input types. */
	if (!mxIsDouble( prhs[0] ) ||
		!mxIsScalar( prhs[0] ) ||
        mxGetScalar( prhs[0] ) < -M_PI ||
        mxGetScalar( prhs[0] ) > M_PI) {
		mexErrMsgIdAndTxt( "VonMisesMixture:sampleVonMises:notScalar",
						   "The distribution mean must be a real scalar value between -pi and pi." );
	}
	
	if (!mxIsDouble( prhs[1] ) ||
		!mxIsScalar( prhs[1] ) ||
		mxGetScalar( prhs[1] ) < 0) {
		mexErrMsgIdAndTxt( "VonMisesMixture:sampleVonMises:notRealPositiveScalar",
						   "The concentration parameter must be a real, positive scalar value." );
	}
	
	/* Check optional third input parameter */
	if ( nrhs != 3 ) {
		/* Set number of samples to default value if not specified */
		numSamples = 1;
	}
	else {
		/* Check type of optional input */
		if (!mxIsDouble( prhs[2] ) ||
			!mxIsScalar( prhs[2] ) ||
			mxGetScalar( prhs[2] ) <= 0) {
			mexErrMsgIdAndTxt( "VonMisesMixture:sampleVonMises:notIntPositiveScalar",
							   "Number of samples must be an integer, positive scalar value." );
		}		
		/* Assign desired number of samples and cast to integer */
		numSamples = (int) mxGetScalar( prhs[2] );
	}
	
	/* Get distribution parameters */
	mu = mxGetScalar( prhs[0] );
	kappa = mxGetScalar( prhs[1] );
	
	/* Allocate output vector */
	plhs[0] = mxCreateDoubleMatrix( numSamples, 1, mxREAL );
	
	/* Get pointer to output vector */
	angles = mxGetPr( plhs[0] );
	
	/* Generate samples */
	sample( mu, kappa, angles, numSamples );
}