//****************************************************//
//         Olivier Catoni & Gautier Appert           //
//         GNU LESSER GENERAL PUBLIC LICENSE        // 
//*************************************************//



#include "initialization.h" 


// [[Rcpp::export]]    

void combineLabels( // Swx[w,x] = Qwx[Rwx[w,x],x]
	IntegerVector &Swx, 
	IntegerVector &Qwx, 
	IntegerVector &Rwx, 
	int S_height, int Q_height, int width) 
{
	int i,j;

	#pragma omp parallel for private(i)
	for (i=0; i<width*S_height; ++i)
	{
		Swx[i] = -1;
	}

	#pragma omp parallel for private(i,j)
	for (i=0;i<width;++i)
	{
		for (j=0;j<S_height;++j)
		{
			if (Rwx[j+S_height*i] >= 0) {
				Swx[j+S_height*i] = Qwx[Rwx[j+S_height*i]+Q_height*i];
			}
		}
	}
}

void mergeParam::merge(int iter)
{
  
  int i, j, int_buf;
  double crit, max_crit, Cmax_i, Cmax_j, double_buf;
  int max_i, max_j;
  max_crit = -1;
	int height = n;
	int width = lastNewLabel;
	int initial_width = firstNewLabel-1; 
  
  #pragma omp parallel private(i)
	{
		double my_max_crit=-1;
    int my_max_i;
    
		#pragma omp for 
		for (i=shift; i < iter-shift; ++i) 
		{
			if (my_max_crit < Crit[i]) { 
				my_max_crit = Crit[i];
				my_max_i = i;
			}
		}
   
    
		#pragma omp critical
		{
			if (max_crit < my_max_crit) { 
				max_crit = my_max_crit;
				max_i = my_max_i;
			}
		}
  
  }

	if (max_crit <= 0) {
		Rcout << "max_crit <= 0 for max_i !\n";
		++shift;
		return;
	}
	Rcout << "max_i = " << max_i + 1 << "\n";
    
	max_crit = -1;
	#pragma omp parallel private(j, crit, i) 
	{
		double my_max_crit = -1;
		int my_max_j;
		#pragma omp for
		for (j=shift; j < iter-shift; ++j) 
		{
			// j!= max_i.
			if (j == max_i) { 
				continue;
			}
			crit = 0;
			for (i=0; i<height; ++i)
			{
				if ( A[i+height*max_i]*A[i+height*j] > 0 ) { 
					crit += 1; 
				}
			}
			crit -= 2;
			if (my_max_crit < crit) { 
				my_max_crit = crit;
				my_max_j = j;
			}
		}
		#pragma omp critical
		{
			if (max_crit < my_max_crit) { 
				max_crit = my_max_crit;
				max_j = my_max_j;
			}
		}
	}

	if (max_crit < 0) { // We are more tolerant for max_j than for max_i
// where we ask for crit > 0. 
		#pragma omp parallel for private(i, int_buf, double_buf)
		for (i=0;i<height;++i)
		{
			int_buf = A[i+height*shift];
			A[i+height*shift] = A[i+height*max_i]; 
			A[i+height*max_i] = int_buf; 
		}
		double_buf = Crit[shift];
		Crit[shift] = Crit[max_i];
		Crit[max_i] = double_buf;
		#pragma omp parallel for private(i, int_buf)
		for (i=0;i<initial_width;++i) {
			int_buf = B[i+initial_width*shift];
			B[i+initial_width*shift] = B[i+initial_width*max_i]; 
			B[i+initial_width*max_i] = int_buf; 
		}
		Rcout << "Merged label " << max_i + 1 << " was moved to " 
			<< shift + 1 << "\n";
		++shift;
		return;
	}

	// display max_j on the terminal
	Rcout << "max_j = " << max_j + 1 << "\n";
	if (max_i > max_j) { 
		Jt[2*(iter-shift)] = max_i;
		Jt[2*(iter-shift)+1] = max_j; 
	} else {
		Jt[2*(iter-shift)] = max_j;
		Jt[2*(iter-shift)+1] = max_i; 
	}
		
	Cmax_i = Cmax_j = 0;
	#pragma omp parallel for private(i) reduction(+:Cmax_i,Cmax_j)
	for (i=0; i<height; ++i)
	{
		// intersection between A_max_i and A_max_j
		if ( A[i+height*max_i]*A[i+height*max_j] > 0 ) { 
			A[i+height*(iter-shift)] = 1; 
                  
			A[i+height*max_i] = 0; 
			A[i+height*max_j] = 0;
		}
		if ( A[i+height*max_i] > 0 ) {
			Cmax_i += 1;
		}
		if ( A[i+height*max_j] > 0 ) {
			Cmax_j += 1; 
		}
	}
	Cmax_i -= 2;
	Cmax_j -= 2;

           //////////////////////////////
          // Update the criterion.    //
         //////////////////////////////
             
	Crit[iter-shift] = max_crit;
	Crit[max_i] = Cmax_i;
	Crit[max_j] = Cmax_j;
              
        
         //////////////////////////////////////
        //       Update matrix B             //
       ////////////////////////////////////////
      
	// merge max_i and max_j.
	#pragma omp parallel for 
	for (i=0; i<initial_width; ++i) // computes patch number iter  
	{
		// UNION between B_max_i and B_max_j
		B[i+initial_width*(iter-shift)] = B[i+initial_width*max_i] 
			+ B[i+initial_width*max_j];
	}
        
	return;
}
         
// [[Rcpp::export]]    

void computeMerge( 
	IntegerVector &A, 
	IntegerVector &Asizes,
	IntegerVector &B,
	IntegerVector &Bsizes,
	IntegerVector &Twx,
	IntegerVector &Jt,
	List &paramList)
{
	int i;
	mergeParam parameters(A, Asizes, B, Bsizes, Twx, Jt, paramList); 
	parameters.initialization();
	for (i=parameters.firstNewLabel-1; i< parameters.lastNewLabel; ++i)
	{
		parameters.merge(i); 
		Rcout << "Merged patch number " << i+1 << " computed \n";
	}
	parameters.output();
}

void mergeParam::output(void) 
{
	int i,j,w;

	#pragma omp parallel for private(i)
	for (i=0;i<n*(firstNewLabel-1);++i)
	{ 
		Twx[i] = -1;
	}
	
	#pragma omp parallel for private(i,j,w)
	for (i=0;i<n;++i)
	{ 
		for (j=0;j<lastNewLabel;++j)
		{
			if (A[i+n*j]>0) { 
				for (w=0;w<(firstNewLabel-1);++w) 
				{
					if (B[w+(firstNewLabel-1)*j] > 0) {
						Twx[w+(firstNewLabel-1)*i] = j;	
					}
				}
			}
		}
	}
}

void mergeParam::initialization(void) 
{
	int i, j;
	#pragma omp parallel for private(i,j)
	for (j=0; j<firstNewLabel-1;++j) 
	{ 
		for (i=0; i<n; ++i) 
		{ 
			Asizes[j] += A[i+n*j];
		}
	}
	#pragma omp parallel for private(i)
	for (i=0;i<(firstNewLabel-1)*lastNewLabel;++i)
	{
		B[i] = 0;
	}
	#pragma omp parallel for private(i)
	for (i=0; i < (firstNewLabel-1); ++i)
	{
		B[i+(firstNewLabel-1)*i] = 1;
	}
	#pragma omp parallel for private(i)
	for (i=0; i<(firstNewLabel-1);++i)
	{
		Bsizes[i] = 1;
	}
	#pragma omp parallel for private(j) 
	for (j=0; j<(firstNewLabel-1); ++j) 
	{
		Crit[j] = Asizes[j] - 2;
	} 
	#pragma omp parallel for private(i)
	for (i=0; i< 2*lastNewLabel; ++i)
	{
		Jt[i] = -1;
	}
}


// [[Rcpp::export]]

void computeContext(IntegerVector &Jt_v, IntegerVector &Cont_v, 
	IntegerVector &ContOrder_v, 
	int number_of_labels , int number_of_merged_labels)
{
	int i, j, k;
	int *Jt = &Jt_v[0];
	int *Cont = &Cont_v[0];
	int *ContOrder = &ContOrder_v[0];
	int ContWeights[number_of_merged_labels];
	int contMax, contMaxIndex;

	#pragma omp parallel for private(i)
	for (i=0; i<number_of_merged_labels;++i)
	{
		ContOrder[i] = -1;
	}

	#pragma omp parallel for private(i)
	for (i=0; i<number_of_merged_labels*number_of_merged_labels; ++i)
	{
		Cont[i] = 0; 
	}
	#pragma omp parallel for private(i)
	for (i=number_of_labels; i<number_of_merged_labels; ++i)
	{
		if (Jt[2*i] >= 0) {
			Cont[Jt[2*i] + number_of_merged_labels*Jt[1+2*i]] = 1;
		// symmetric pairs : 
			Cont[Jt[1+2*i] + number_of_merged_labels*Jt[2*i]] = 1;
		}
	}

	contMax = -1;
	#pragma omp parallel private(i,j)
	{
		int my_max = -1;
		int my_max_index;
		#pragma omp for 
		for (i=0; i<number_of_merged_labels; ++i)
		{
			ContWeights[i] = 0;
			for (j=0; j<number_of_merged_labels; ++j)
			{
				ContWeights[i] += Cont[i+number_of_merged_labels*j];		
			}
			if (ContWeights[i] < my_max) {
				my_max = ContWeights[i];
				my_max_index = i;
			}
		}
		#pragma omp critical
		{
			if (contMax < my_max) {
				contMax = my_max;
				contMaxIndex = my_max_index;
			}
		}
	}
	for (k=0; k<number_of_merged_labels; ++k) {
		ContOrder[k] = contMaxIndex;
		#pragma omp parallel for private(i)
		for (i=0; i<number_of_merged_labels; ++i)
		{
			ContWeights[i] -= 
				Cont[i + number_of_merged_labels*contMaxIndex];	
			Cont[i + number_of_merged_labels*contMaxIndex] = 0;	
		}
		ContWeights[contMaxIndex] = 0;
		contMax = -1;
		#pragma omp parallel private(i)
		{
			int my_max = -1;
			int my_max_index;
			#pragma omp for
			for (i=0; i<number_of_merged_labels; ++i)
			{
				if (my_max < ContWeights[i]){
					my_max = ContWeights[i];
					my_max_index;
				}
			}
			#pragma omp critical
			{
				if (contMax < my_max) { 
					contMax = my_max;
					contMaxIndex = my_max_index;
				}
			}
		}
		if (contMax <= 0) { 
			break;
		}
	}
}

// [[Rcpp::export]]

void applySyntax(IntegerVector &Twx_v, 
	IntegerVector &Gtt_v, IntegerVector &Jt_v, 
	IntegerVector &Gt_v, IntegerVector &Uwx_v, 
	IntegerVector &Axw_v,
	IntegerVector &Axu_v, 
	int number_of_syntax_labels,
	int number_of_labels, 
	int number_of_patches,
	int number_of_images)
{
	int *Twx = &Twx_v[0];
	int *Gtt = &Gtt_v[0];
	int *Gt = &Gt_v[0];
	int *Jt = &Jt_v[0];
	int *Uwx = &Uwx_v[0];
	int *Axw = &Axw_v[0];
	int *Axu = &Axu_v[0];
	int i, j, buf1, buf2;

	#pragma omp parallel for private(i)
	for (i=0; i<number_of_patches; ++i)
	{
		Gt[i] = i;
	}
	#pragma omp parallel for private(i,buf1,buf2)
	for (i=number_of_patches; i<number_of_labels; ++i)
	{
		if (Jt[2*i] >= 0) { 
// reversed order since
// Gtt is the analogous of Twx (and not Txw)
			buf1 = Gtt[Jt[1+2*i] + number_of_labels*Jt[2*i]]; 
			buf2 = Gtt[Jt[2*i] + number_of_labels*Jt[1+2*i]]; 
			Gt[i] = (buf2 < buf1) ? buf1 : buf2; 
		} else {
			Gt[i] = i;
		}
	}
	#pragma omp parallel for private(i)
	for (i=0; i<number_of_patches*number_of_images; ++i) 
	{
			Uwx[i] = Twx[i] < 0 ? -1 : Gt[Twx[i]];
	}
	#pragma omp parallel for private(i)
	for (i=0; i<number_of_images*number_of_syntax_labels; ++i)
	{
		Axu[i] = 0;
	}
	#pragma omp parallel for private(i,j)
	for (i=0; i<number_of_images; ++i)
	{
		for (j=0; j<number_of_labels; ++j)
		{
			if (Axw[i + number_of_images * j] == 1) {
				Axu[i+number_of_images*Gt[j]] = 1;
			}
		}
	}
}
