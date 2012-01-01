#include <R.h>
#include <Rmath.h>
#include <math.h>
#include <stdlib.h>

void cpmMLEBartlett(double *S, double *W, int *n, double *Ds)
{
	int i,j;
	double mu1,mu2,sigma1,sigma2,sigma,V1,V2,C,G;

    //Rprintf("Called\n");
	for (i = 1; i < (*n-2); i++) {
		j = i+1; //we have j obserations in the left sample
		
        mu1 = S[j-1]/(double)j; 
		mu2 = (S[*n-1]-S[j-1])/(double)(*n-j);
		V1 = W[j-1];
		V2 = ( W[*n-1] - W[j-1] - (j*(*n-j)*(mu1-mu2)*(mu1-mu2))/ (double)(*n) );
		sigma1 = V1/(double)(j-1);
		sigma2 = V2/(double)(*n-j-1);
		sigma = (V1+V2)/(double)(*n-2);	
		C = 1 + ( 1/(double)(j-1) + 1/(double)(*n-j-1) - 1/(double)(*n-2))/(double)3;
        
		G = ((j-1)*log(sigma/sigma1) + (*n-j-1)*log(sigma/sigma2))/C;
		Ds[i] = G;
	}
}

void cpmMLEJoint(double *S, int *nS, double *W, int *nW, int *n, int *nb, double *Ds)
{
	int i,n1,n2;
	double Sok,Son,Skn,mu1,mu2,C,G;


	for (i = 1; i < (*nS-2); i++) {
		n1 = i+1; n2 = *n - n1;
		mu1 = S[n1-1]/n1; 
		mu2 = (S[*n-1]-S[n1-1])/n2;
		
		Sok = W[n1-1]/(double)n1;
		Son = W[*n-1]/(double)*n;
		Skn = W[*n-1] - W[n1-1] - (n1*(*n-n1)*(mu1-mu2)*(mu1-mu2))/ *n;
		Skn = Skn/(double)n2;
		C = 1 + 11/(double)12 * (1/(double)n1 + 1/(double)n2 - 1/(double)*n) + (1/(double)(n1*n1) + 1/(double)(n2*n2) - 1/(double)(*n * *n));
		G = (n1*log(Son/Sok) + n2*log(Son/Skn))/C;

		Ds[i] = G;
		//Rprintf("%d: %f \n ",i,Ts[i]);	
	}
}


void cpmMLEJointAdjusted(double *S, int *nS, double *W, int *nW, int *n, int *nb, double *Ds)
{
    int i,n1,n2;
	double Sok,Son,Skn,mu1,mu2,C,G;


	for (i = 1; i < (*nS-2); i++) {
		n1 = i+1; n2 = *n - n1;
		mu1 = S[n1-1]/n1; 
		mu2 = (S[*n-1]-S[n1-1])/n2;
		
		Sok = W[n1-1]/(double)n1;
		Son = W[*n-1]/(double)*n;
		Skn = W[*n-1] - W[n1-1] - (n1*(*n-n1)*(mu1-mu2)*(mu1-mu2))/ *n;
		Skn = Skn/(double)n2;
		C = 1 + 11/(double)12 * (1/(double)n1 + 1/(double)n2 - 1/(double)*n) + (1/(double)(n1*n1) + 1/(double)(n2*n2) - 1/(double)(*n * *n));
		G = (n1*log(Son/Sok) + n2*log(Son/Skn))/C;

		Ds[i] = G;
		//Rprintf("%d: %f \n ",i,Ts[i]);	
	}

    if (*nS < 10) {return;}

    //now do the adjustment
    int len = *nS;
    double meanAdjustments[3] = {2.2989, 2.0814, 2.0335};
    double sdAdjustments[3] = {2.3151,2.0871, 2.0368};

    Ds[1] = (Ds[1] - meanAdjustments[0])/sdAdjustments[0];
    Ds[len-3] = (Ds[len-3] - meanAdjustments[0])/sdAdjustments[0];
    Ds[2] = (Ds[2] - meanAdjustments[1])/sdAdjustments[1];
    Ds[len-4] = (Ds[len-4] - meanAdjustments[1])/sdAdjustments[1];
    Ds[3] = (Ds[3] - meanAdjustments[2])/sdAdjustments[2];
    Ds[len-5] = (Ds[len-5] - meanAdjustments[2])/sdAdjustments[2];

    //now convert back to chi-square
    double truemean = 2;
    double truesd = 2;

    Ds[1] = (Ds[1] * truesd) + truemean;
    Ds[2] = (Ds[2] * truesd) + truemean;
    Ds[3] = (Ds[3] * truesd) + truemean;
    Ds[len-3] = (Ds[len-3] * truesd) + truemean;
    Ds[len-4] = (Ds[len-4] * truesd) + truemean;
    Ds[len-5] = (Ds[len-5] * truesd) + truemean;
}

//the W here is the V in hte paper i think 




void cpmMLECVM(double *X, int *nX, int *ranks, double *Ds)
{
	int i,j;
	double a,b,statistic,mu,sigma,N,prod,n0,n1;

	double *cumsums;
	cumsums = malloc(*nX * sizeof(double));
	
	N = *nX;

	mu = (double)1/6 + 1/(6*N);

	for (i = 1; i < (*nX-2) ;i++) {
		n0 = i+1;
		n1 = *nX-i-1;
		//Rprintf("%d: %f %f %d\n ",i,n0,n1,*nX);
		a = 1/n0;
		b = -1/n1;

		
		for (j = 0 ; j < *nX ; j++) {
			if (ranks[j] <= n0) {
				cumsums[j] = a;
			}
			else {
				cumsums[j] = b;
			}			
		}

		for (j = 1 ; j < *nX; j++) {
			cumsums[j] = cumsums[j-1] + cumsums[j]; 
		}
		statistic = 0;
		for (j = 0; j < *nX; j++) {
			statistic+= cumsums[j]*cumsums[j];
		}
		prod = n0*n1;
		sigma = sqrt(  (double) 1/45 * (N+1)/(N*N)  * (4*prod*N - 3*(n1*n1 + n0*n0) - 2*prod)/(4*prod));

		//Us[i] = (statistic  - mu)/sigma;
		Ds[i] = (statistic * (prod)/(N*N) - mu)/sigma; 
	}
	free(cumsums);
}


void cpmMLEKS(double *X, int *nX, int *orders, int *pvalues, int *adjustment, double *Ds)
{
	int i,j;
	double n0,n1,a,b,statistic, temp,z,correction;

	double *cumsums;
	cumsums = malloc(*nX * sizeof(double));
	
	
	for (i = 1; i < (*nX-2) ;i++) {
		n0 = i+1;
		n1 = *nX-i;
		//Rprintf("%d: %f %f %d\n ",i,n0,n1,*nX);
		a = 1/n0;
		b = -1/n1;

		
		for (j = 0 ; j < *nX ; j++) {
			if (orders[j] <= n0) {
				cumsums[j] = a;
			}
			else {
				cumsums[j] = b;
			}			
		}

		for (j = 1 ; j < *nX; j++) {
			cumsums[j] = cumsums[j-1] + cumsums[j]; 
		}
		statistic = 0;
		for (j = 0; j < *nX; j++) {
			temp = fabs(cumsums[j]);
			if(temp > statistic) {
				statistic=temp;
			}
		}
		
		//return only the statistic, not the pvalue
		if (*pvalues==0) {
			Ds[i] = statistic;
			continue;
		}	

		//now we compute p value
		//this is a continuity correction for the KS statistic so I can more accurately use the asymptotic distribution. See Kim 1969
		correction = 0;
		if (*adjustment > 0) {
			//rearrange so n0 >= n1
			if (n1 > n0) {
				a = n1;
				n1 = n0;
				n0 = a;
			} 
		
			if (n0 > 2*n1) {
				correction=1/(2*sqrt(n0));
			}
			else {
				if ((int) n0 % (int) n1 == 0) {
					correction=2/(3*sqrt(n0));
				}
				else {
					correction=2/(5*sqrt(n0));
				}
			}
		}
		z = statistic*sqrt((n0*n1)/(n0+n1)) + correction;

		z = z*z;
		Ds[i] = 2*(exp(-2*z) - exp(-8 * z));
	}
	free(cumsums);
}





double square(double x) {
        return(x*x);
}

void cpmMLEMood(double *X, int *nX, int *N, int *nN, int *ranks, int *nranks, double *Ds)
{
        int i;
        double n,n0,n1,M,mu,sd,med;
        double *cumsumsSquare;

        n=N[*nN-1];
        med = (n+1)/2;

        /*cumsumsSquare[i] is the sum of squared ranks (-median) for i,i+1,i+2,...,n*/
        cumsumsSquare = malloc(*nranks * sizeof(double));

       cumsumsSquare[0] =  (ranks[0] - med)*(ranks[0] - med);

        for (i = 1 ; i < *nranks; i++) {
                cumsumsSquare[i] = cumsumsSquare[i-1] + (ranks[i]-med)*(ranks[i]-med);
        }

        for (i = 1; i < (*nX-2) ;i++) {
                //n0 = N[i]-1;
                //n1 = N[*nN-1]-n0;
                n0 = i+1;
                n1 = n - n0;
                M = cumsumsSquare[i];

                mu= n0*(n*n - 1)/12;
                sd = sqrt(n0*n1*(n+1)*(n*n-4)/180);

                Ds[i] = (M-mu)/sd;

        }
        free(cumsumsSquare);
}




void cpmMLEMW(double *X, int *nX, int *N, int *nN, int *ranks, int *nranks, double *Ds)
{
	int i;
	double n0,n1,R1,U,mu,sd;
    int n = N[*nN-1];
    
	double *cumsums;
	cumsums = malloc(*nranks * sizeof(double));
	
	cumsums[0] =  ranks[0];
	
	for (i = 1 ; i < *nranks; i++) {
		cumsums[i] = cumsums[i-1] + ranks[i];
	}

	for (i = 1; i < (*nX-2) ;i++) {
		n0 = i+1;
        n1 = n-n0;
		R1 = cumsums[i];
		
		U = R1 - ( n0*(n0+1)/2) ;		
				
		mu = n0*n1/2;			
		sd = sqrt(n0*n1*(n0+n1+1)/12);
		//Rprintf("%d: %f %f %f\n ",i,U,mu,sd);
		Ds[i] = (U-mu)/sd;		
	}
	free(cumsums);
}


void cpmMLELepage(double *X, int *nX, int *N, int *nN, int *ranks, int *nranks, double *Ds) {
	double *DsMood = malloc(*nX * sizeof(double));
	
	cpmMLEMW(X,nX,N,nN,ranks,nranks,Ds);
	cpmMLEMood(X,nX,N,nN,ranks,nranks,DsMood);
	
	
	for (int i = 1; i < *nX-2; i++) {
		Ds[i] = Ds[i]*Ds[i] + DsMood[i]*DsMood[i];
	}
	
	free(DsMood);
}




void cpmMLEStudent(double *S, int *nS, double *W, int *nW, int *n, int *nb, double *Ds)
{
	int i,j;
	double temp,E,sigma;	//sigma is standard deiation of test stasitic
	sigma = (*nb + *nS) -2;	//degrees of freedom
	sigma = (sqrt(sigma/(sigma-2)));


	for (i = 1; i < *nS-2; i++) {
		j = i+1; 		
		temp = *n*S[i] - j*S[*nS-1];
		E = temp*temp / (*n*j*(*n-j));
		Ds[i] = sqrt( (*n-2)*E / (W[*nW-1]-E)) / sigma;
		//Rprintf("%d: %f \n ",i,Ts[i]);	
	}
}

void cpmMLEFET(double *S, int *nS, double *N, int *nN, int *n, double *lambda, double *Ds)
{
	int i,n0,n1,s0,s1;


	for (i = 1; i < *nS; i++) {
		n0 = N[i-1];
		n1 = *n-n0;
		s0 = S[i-1];
		s1 = S[*nS-1] - s0;
		//Ds[i-1] =  1.0 - my_gsl_cdf_hypergeometric_P(s0,s0+s1,n0+n1-s0-s1,n0);
        Ds[i-1] =  1.0 - phyper(s0,s0+s1,n0+n1-s0-s1,n0,1,0);
	}

	if (*nS > 3 && *lambda > 0) {
		for (i = 2; i < *nS; i++) { 
			Ds[i] = (1 - *lambda)*Ds[i-1] + *lambda * Ds[i];
		}
	}
		

}


