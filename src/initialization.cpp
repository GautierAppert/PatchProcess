#include "initialization.h" 

// [[Rcpp::export]]

void computeLoop (	
	NumericVector &m, 
	NumericVector &v, 
	NumericVector &w, 
	NumericVector &C, 
	IntegerVector &A, 
	IntegerVector &B, 
	IntegerVector &Wsx,
	List &paramList) 
{
	int i, j;

	// creates a Patches object to pass parameters from R
	Patches parameters(m, v, w, C, A, B, Wsx, paramList); 

	parameters.computeBeta();
	paramList["thresholdValue"] = parameters.thresholdValue;
	// loop on patch index
	for ( i = parameters.firstNewLabel- 1; i < parameters.lastNewLabel; ++i ) 
	{
		parameters.split(i); 
		Rcout << "patch number " << i+1-parameters.shift << " computed \n";
		Rcout << "at iteration number " << i+1 << "\n";
	}
	parameters.output();
	Rcout << "betaCoeff = " << parameters.betaCoeff << "\n";
}

void Patches::output(void) {
	int i,j,s;
	double area;

	#pragma omp parallel for private(i,j,s,area)
	for (j=0;j<lastNewLabel;++j) 
	{ 
		for (i=0; i<n; ++i) 
		{ 
			if (A[i+n*j]>0) { 
				for (s=0; s<d; ++s) { 
					if (B[s+d*j]>0) {
						Wsx[s+d*i] = j;
					}
				}
			}
		}
	} 
}


// compute the threshold from the mean variance of the pixels.
void Patches::computeBeta(void) {
	double m_buf, v_buf, square_buf;
	int s, i;

	#pragma omp parallel for private(i)
	for (i=0; i<n*n; ++i) { 
		A[i] = 0;
	}	
	#pragma omp parallel for private(i)
	for (i=0; i< n; ++i) { 
		A[i+n*i] = 1;
	}
	#pragma omp parallel for private(i)
	for (i=0;i<d*n;++i) { 
		B[i] = 1;
	}
	
	v_buf = 0;
	#pragma omp parallel for private(i, m_buf, square_buf) reduction(+:v_buf)
	for (s=0; s<d; ++s)
	{
		m_buf = 0;
		for (i=0; i<n; ++i)
		{		
			m_buf += m[s+d*i]; 
		}
		m_buf /= n;
		for (i=0; i<n; ++i)
		{
			square_buf = m[s+d*i] - m_buf;
			square_buf *= square_buf;
			v_buf += square_buf; 
		}
	}
	v_buf /= d*n;
	thresholdValue = thresholdCoeff * v_buf;
	Rcout << "sqrt(mean variance) = " << sqrt(v_buf) << "\n";
	Rcout << "sqrt(threshold) = " << sqrt(thresholdValue) << "\n";
}

void Patches::split(int iter) { 
	int i,j,s;
	double m_buf, v_buf, weight_buf, square_buf, a, double_buf;
	double crit, max_crit, Cmax_i, Cmax_j;
	int max_i, count, max_j, int_buf;

	// computes arg max C
	max_crit = -1;
	#pragma omp parallel private(i)  
	{
		double my_max_crit = -1;
		int my_max_i;
		#pragma omp for
		for (i=shift; i<iter-shift; ++i) 
		{ 
			if ( my_max_crit < C[i] ) {
				my_max_crit = C[i]; // selection according to the criterion.
				my_max_i = i;
			}
		}
		#pragma omp critical
		{
			if ( max_crit < my_max_crit ) {
				max_crit = my_max_crit;
				max_i = my_max_i;
			}
		}
	}

	if (max_crit < 0) {
		Rcout << "Nothing to split any more.\nLast split = " << iter-shift << "\n"; 
		++shift;
		return;
	}

	Rcout << "max_i = " << max_i + 1 << "\n";
	// computes C for each j != max_i
	max_crit = -1;
	#pragma omp parallel private(j, count, s, weight_buf, crit, m_buf, \
		 square_buf, v_buf)
	{
		double my_max_crit = -1;
		int my_max_j; 
		#pragma omp for  
		for (j=shift; j<iter-shift; ++j)
		{
			if (j == max_i) continue; // we are looking for j != max_i 
			count = 0;
			// compute the intersection between B_max_i and B_j.
			for (s=0; s<d; ++s) 
			{
				count += B[s+d*max_i] * B[s+d*j];	
			}
			if (count == 0) continue;
			
			// compute the sum of weights.
			weight_buf = w[max_i] + w[j];
			
			//  compute the criterion for j.
			crit = 0;
			
			// for each pixel do
			// compute the criterion using the law of total variance.
			for (s=0; s<d; ++s) 
			{
			  // Test if there is an intersection.
				if (B[s+d*max_i]*B[s+d*j] == 0) continue;	
				
				// compute the mean of the mean.
				m_buf = (w[max_i]*m[s+d*max_i] 
					+ w[j]*m[s+d*j])/weight_buf; 
				square_buf = m[s+d*max_i] - m_buf;
				
				// compute the variance of the means.
				square_buf *= square_buf; 
				
				// compute the means of the variance.
				v_buf = w[max_i] * ( v[s+d*max_i] + square_buf );
				square_buf = m_buf - m[s+d*j];
				square_buf *= square_buf; 
				v_buf += w[j] * ( v[s+d*j] + square_buf); 
				v_buf /= weight_buf;
				if (v_buf < thresholdValue) { 
					crit += 1; 
				} 
			}
			crit *= betaCoeff;
			crit -= weight_buf;
			if (my_max_crit < crit) {
				my_max_crit = crit;
				my_max_j = j;
			}
		}
		
		// get our max_j.
		#pragma omp critical
		{
			if (max_crit < my_max_crit) { 
				max_crit = my_max_crit;
				max_j = my_max_j;
			}
		}
	}
		 
		 
	if (max_crit < 0) {
		#pragma omp parallel for private(s,double_buf,int_buf)
		for (s=0;s<d;++s) 
		{
			double_buf = m[s+d*shift]; 
			m[s+d*shift] = m[s+d*max_i]; 
			m[s+d*max_i] = double_buf; 
			double_buf = v[s+d*shift]; 
			v[s+d*shift] = v[s+d*max_i]; 
			v[s+d*max_i] = double_buf; 
			int_buf = B[s+d*shift]; 
			B[s+d*shift] = B[s+d*max_i]; 
			B[s+d*max_i] = int_buf; 
		}
		double_buf = w[shift];
		w[shift] = w[max_i];
		w[max_i] = double_buf; 
		double_buf = C[shift];
		C[shift] = C[max_i];
		C[max_i] = double_buf;
		#pragma omp parallel for private(i,int_buf)
		for (i=0;i<n;++i)
		{
			int_buf = A[i+n*shift];
			A[i+n*shift] = A[i+n*max_i];
			A[i+n*max_i] = int_buf;
		}
		++shift;
		return;
	}

	Rcout << "max_j = " << max_j + 1 << "\n";
	C[iter-shift] = max_crit;
	Cmax_i = Cmax_j = 0;
	weight_buf = w[max_i] + w[max_j];
	
	#pragma omp parallel for private(s, square_buf) reduction(+:Cmax_i, Cmax_j) 
	for (s=0; s<d; ++s) // computes patch number iter  
	{
	  
	  // compute patch mean using mean of m[s,i] and m[s,j].
		m[s + d*(iter-shift)] = (w[max_i]*m[s+d*max_i] 
			+ w[max_j]*m[s+d*max_j])/weight_buf; 
	  
		square_buf = m[s + d*max_i] - m[s + d*(iter-shift)];
		square_buf *= square_buf; 
		v[s+d*(iter-shift)] = w[max_i] * ( v[s+d*max_i] + square_buf ); 
		square_buf = m[s+d*max_j] - m[s+d*(iter-shift)];
		square_buf *= square_buf; 
		v[s+d*(iter-shift)] += w[max_j] * ( v[s+d*max_j] + square_buf); 
		v[s+d*(iter-shift)] /= weight_buf;
		
		// Update A and B 
		// we need to compute the intersection
		// and check that the variance is under thresholdValue.
		if ( (B[s+d*max_i]*B[s+d*max_j] > 0) && (v[s+d*(iter-shift)] 
			< thresholdValue) ) 
		{
			B[s+d*(iter-shift)] = 1;
			B[s+d*max_i] = 0;
			B[s+d*max_j] = 0;
		} else {  
			B[s+d*(iter-shift)] = 0;
		}
		if (B[s+d*max_i] > 0) {
			Cmax_i += 1; 
		}
		if (B[s+d*max_j] > 0) {
			Cmax_j += 1;
		}
	}
	
	C[max_i] = Cmax_i*betaCoeff - w[max_i];
	C[max_j] = Cmax_j*betaCoeff - w[max_j];
	w[iter-shift] = w[max_i] + w[max_j];
	#pragma omp parallel for private(i)
	for (i=0; i<n; ++i)
	{
		A[i+n*(iter-shift)] = A[i+n*max_i] + A[i+n*max_j];
	}
	return;
}
