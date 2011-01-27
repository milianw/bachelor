#include <iostream>
#include <string>
#include <cmath>
#include <bitset>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_const_mksa.h>
#include <gsl/gsl_complex.h>
#include <gsl/gsl_complex_math.h>
#include <gsl/gsl_eigen.h>
#include <gsl/gsl_sort_vector.h>


using namespace std;

//Physical constants, all in MKS units
const double hbar = GSL_CONST_MKSA_PLANCKS_CONSTANT_HBAR;
const double h= 6.62606896E-34;
const double NUC_MAGNETON = GSL_CONST_MKSA_NUCLEAR_MAGNETON; // (J/T)
const double g_1H = 5.585694701 ;            // proton g factor
const double g_14N = 0.403761;               // N14 g factor
const double GAMMA_1H = 267.513E6 ;          // (rad/s.T)
const double GAMMA_14N =  19.331E6 ;         // (rad/s.T)
const double Bohrm = 9.27400949E-24;         // Bohr magneton in J/T
// Reminder: GAMMA_1H * hbar =  NUC_MAGNETON * g_1H
//const double B = 1;                          // Static magnetic field (T)
//const double B2 = ;                        // RF field
//const double LARMOR_1H = B * GAMMA_1H /2/M_PI / 1.0E6;
const gsl_complex half_x[4] = { gsl_complex_rect(0,0), gsl_complex_rect(0.5,0),
				gsl_complex_rect(0.5,0), gsl_complex_rect(0,0)};
const gsl_complex half_y[4] = { gsl_complex_rect(0,0), gsl_complex_rect(0,-0.5),
				gsl_complex_rect(0,0.5), gsl_complex_rect(0,0)};
const gsl_complex half_z[4] = { gsl_complex_rect(0.5,0), gsl_complex_rect(0,0),
				gsl_complex_rect(0,0), gsl_complex_rect(-0.5,0)};
const int nprotons = 2;

void hamiltonian(const double B);
void printMatrix(const gsl_matrix_complex*, const string, const int);
void printMatrix(const gsl_matrix*, const string , const int);
void printVector(const gsl_vector_complex *vector, const string name, const int dim);
void printVector(const gsl_vector *vector, const string name, const int dim);



int main(int argc, char* argv[])
{
  //cout << 2.023 * Bohrm << endl;
  //cout << NUC_MAGNETON << endl;
  //cout << g_1H * NUC_MAGNETON << endl;
  //cout << GAMMA_1H * hbar << endl;
  
  /////////////////////////////////////////////////////////////////////////
  // The input file should be structured as follows:	      		 //
  // natoms n                     (number of coupled nuclei)    	 //
  // field x y z                  (The static field direction)  	 //
  // g                            (g tensor)                   		 //
  //   xx xy xz								 //
  //   yx yy yz								 //
  //   zx zy zz								 //
  // N                            (A tensor)				 //
  //   xx xy xz								 //
  //   yx yy yz								 //
  //   zx zy zz								 //
  // H                            (A tensor)				 //
  //   .........							 //
  //   .........							 //
  //   .........							 //
  // And so on...					     		 //
  /////////////////////////////////////////////////////////////////////////
  

  //Parsing input file
  //ifstream inputfile;
  //inputfile.open(argv[1]);
  // if(!inputfile)
  //   {
  //     cerr << "This program expects a single input file\n";
  //     cerr << "Check the comments in the source code for details\n"
  //     return(1);
  //   }


  //  for(int i=0;i<100;i++)
  //{

//   hamiltonian(0.362562);
  hamiltonian(0.3417757);
  //}

  
  return 0;
}



void hamiltonian(const double B)
{
  
  
  //Temporary values for testing===============================================
  //gtensor=============================
  double g[3][3]={2.0022838, 0.0, 0.0,
		  0.0, 2.0022838, 0.0,
		  0.0, 0.0, 2.0022838};
  gsl_matrix_complex* gtensor = gsl_matrix_complex_calloc(3,3);
  for(int i=0;i<3;i++)
    {
      for(int j=0; j<3; j++)
	gsl_matrix_complex_set(gtensor,i, j, gsl_complex_rect(g[i][j],0));
    }
  //static field========================
  gsl_vector_complex* Bstatic = gsl_vector_complex_calloc(3);
  {
    double B0_norm = B; //0.2839; //field in Tesla
    double B0[3]={0,0,1}; //field direction
    double B_temp=sqrt(B0[0]*B0[0]+B0[1]*B0[1]+B0[2]*B0[2]);
    for(int i=0;i<3;i++)
      B0[i] = B0[i]*B0_norm/B_temp;
    //cout << " " << B0[0] << '\t' << B0[1] << '\t' << B0[2] << endl;
    //cout.precision(5);
    //cout << "\nStatic Field (T): " << sqrt(B0[0]*B0[0]+B0[1]*B0[1]+B0[2]*B0[2]) << endl;
    for(int i=0;i<3; i++)
      gsl_vector_complex_set(Bstatic, i, gsl_complex_rect(B0[i],0));
  }
  //arbitrary A-tensor==================
  double a1[3][3]={1420, 0.0, 0.0,
		   0.0, 1420, 0.0,
		   0.0, 0.0, 1420}; //proton hyperfine coupling in T
  gsl_matrix_complex* atensor = gsl_matrix_complex_calloc(3,3);
  for(int i=0;i<3;i++)
    {
      for(int j=0; j<3; j++)
	gsl_matrix_complex_set(atensor,i, j, gsl_complex_rect(a1[i][j],0));
    }

  
  //define Pauli matrices======================================================
  gsl_matrix_complex* half_xv = gsl_matrix_complex_alloc(2,2);
  gsl_matrix_complex* half_yv = gsl_matrix_complex_alloc(2,2);
  gsl_matrix_complex* half_zv = gsl_matrix_complex_alloc(2,2);

  
  for(int i=0;i<2;i++)
    {
      for(int j=0;j<2;j++)
  	gsl_matrix_complex_set(half_xv,i,j,half_x[2*i+j]);
    }
  for(int i=0;i<2;i++)
    {
      for(int j=0;j<2;j++)
  	gsl_matrix_complex_set(half_yv,i,j,half_y[2*i+j]);
    }
  for(int i=0;i<2;i++)
    {
      for(int j=0;j<2;j++)
  	gsl_matrix_complex_set(half_zv,i,j,half_z[2*i+j]);
    }
  
  

  //Generate the direct product eigenstates====================================

  const int dim=pow(2,nprotons+1);
  //An array of binary bits to represent the spin states
  bitset<(nprotons+1)> states[dim];  //A binary representation of the states
  for(int i=0;i<dim;i++)           //This includes all spin half species (protons+1 electron)
    states[i]=i;                   //for each bit: 0=+0.5, 1=-0.5
  bitset<2> nitrogen[3];           //The -1,0,+1 hydrogen states
  for(int i=0;i<3;i++)             //00=-1, 01=0, 10=+1
    nitrogen[i]=i;
  //  cout << nitrogen[0] << '\t' << nitrogen[1] << '\t' << nitrogen[2] << endl;
  /*cout << "\nSpin states: \n";
  for(int i=0;i<dim;i++)
    cout << i+1 << '\t' << states[i] << endl;
  //Debug: Output the states===================================================
  cout << "\nProduct spin states: \n";
  cout << "======================\n";
  for(int i=0;i<dim;i++)
    cout << '\t' << states[i];
  cout << endl;
  for(int i=0;i<dim;i++)
    cout << "----------";
  cout << endl;
   */
  
  //cout << "!!!!!!!!!!!!!!!!!!!!\n";
  //for (int i =0; i<nprotons+1;i++)
  //cout << states[6][i] << endl;
  //cout << "!!!!!!!!!!!!!!!!!!!!\n";
  
  
  //Allocate Hamiltonian matrix
  gsl_matrix_complex* Ham = gsl_matrix_complex_calloc(dim,dim);
  gsl_vector_complex* s = gsl_vector_complex_calloc(3);
  short int a=0,b=0, flag=0;






  
  //Compute nZeeman============================================================  
  gsl_matrix_complex* nZeeman = gsl_matrix_complex_calloc(dim,dim);
  gsl_vector_complex* I = gsl_vector_complex_calloc(3);
  gsl_complex BI; //dot product of I and B
  for(int i=0; i<dim; i++)
    {
      for(int j=0; j<dim; j++)
  	{
	  if(states[i][nprotons]!=states[j][nprotons]) continue;  //matrix elements between different e states are zero
	  for(int k=0; k<nprotons;k++)    //nprotons is always the index of the electronic spin state
	    {
	      flag=0;
	      for(int l=0;l<nprotons;l++)
		{
		  if(l==k) continue;
		  if(states[i][l]!=states[j][l])
		    {
		      flag =1;
		    }
		}
	      if(flag==1) continue;
	      a= states[i][k];  //spin state of state k in row i
	      b= states[j][k];  //spin state of state k in column j
	      gsl_vector_complex_set(I, 0, gsl_matrix_complex_get(half_xv,a,b));
	      gsl_vector_complex_set(I, 1, gsl_matrix_complex_get(half_yv,a,b));
	      gsl_vector_complex_set(I, 2, gsl_matrix_complex_get(half_zv,a,b));
	      gsl_blas_zdotu(Bstatic, I, &BI);
	      //BI = gsl_complex_mul_real(BI, -1.0);
	      gsl_matrix_complex_set(nZeeman,i,j,gsl_complex_add(gsl_matrix_complex_get(nZeeman,i,j),BI));
	      BI = gsl_complex_rect(0,0);
	    }
	}
    }
  flag=0;
  gsl_matrix_complex_scale(nZeeman, gsl_complex_rect(-1.0*g_1H*NUC_MAGNETON,0));
  //To turn of the nZeeman term
  //gsl_matrix_complex_scale(nZeeman, gsl_complex_rect(0,0));
  gsl_matrix_complex_add(Ham, nZeeman);
  //Debug: output nZeeman=====================================================
  //printMatrix(nZeeman, "nZeeman", dim);
  
  

  
  //Compute Hyperfine couplings matrix=========================================  
  gsl_matrix_complex* hyperfine = gsl_matrix_complex_calloc(dim,dim);
  gsl_vector_complex* atensor_I= gsl_vector_complex_calloc(3);
  gsl_complex sAi;
  /*///
  gsl_vector_complex_set(s, 0, gsl_matrix_complex_get(half_xv,0,0));
  gsl_vector_complex_set(s, 1, gsl_matrix_complex_get(half_yv,0,1));
  gsl_vector_complex_set(s, 2, gsl_matrix_complex_get(half_zv,0,0));
  printVector(s, "a", 3);
  gsl_blas_zdotu(s, s, &sAi);
  cout << "zdotu(s, s): " << GSL_REAL(sAi) << ',' << GSL_IMAG(sAi) << endl;
  gsl_blas_zdotu(s, s, &sAi);
  cout << "zdotc(s, s): " << GSL_REAL(sAi) << ',' << GSL_IMAG(sAi) << endl;
  return;
  ///*/
  for(int i=0; i<dim; i++)
    {
      for(int j=0; j<dim; j++)
  	{
    	  a= states[i][nprotons];   //compute elements of s vector
  	  b= states[j][nprotons];
  	  gsl_vector_complex_set(s, 0, gsl_matrix_complex_get(half_xv,a,b));
  	  gsl_vector_complex_set(s, 1, gsl_matrix_complex_get(half_yv,a,b));
  	  gsl_vector_complex_set(s, 2, gsl_matrix_complex_get(half_zv,a,b));
  	  for(int k=0; k<nprotons; k++)    //loop over nuclei
  	    {
	      flag=0;
	      for(int l=0;l<nprotons;l++)
		{
		  if(l==k) continue;
		  if(states[i][l]!=states[j][l])
		    {
		      flag =1;
		    }
		}
	      if(flag==1) continue;
  	      a= states[i][k];    //compute elements of I vector
  	      b= states[j][k];
  	      gsl_vector_complex_set(I, 0, gsl_matrix_complex_get(half_xv,a,b));
  	      gsl_vector_complex_set(I, 1, gsl_matrix_complex_get(half_yv,a,b));
  	      gsl_vector_complex_set(I, 2, gsl_matrix_complex_get(half_zv,a,b));
	      //skip this for now, all nuclei assigned the free proton isotropic coupling
  	      //multiply atensor by I
//               cout << "\n\n~~~\n\n" << endl;
//               printMatrix(atensor, "atensor", 3);
//               printVector(I, "I", 3);
  	      gsl_blas_zgemv(CblasTrans, gsl_complex_rect(1,0), atensor, I, gsl_complex_rect(0,0), atensor_I);
  	      //multiply s by atensor_I
	      gsl_blas_zdotu(atensor_I, s, &sAi);
	      gsl_matrix_complex_set(hyperfine,i,j,gsl_complex_add(gsl_matrix_complex_get(hyperfine,i,j),sAi));
	      gsl_vector_complex_scale(atensor_I, gsl_complex_rect(0,0));
	    }
  	}
    }
  gsl_matrix_complex_scale(hyperfine, gsl_complex_rect(h*1.0E6,0));
  //To Turn off the hyperfine term
  gsl_matrix_complex_add(Ham, hyperfine);
  //Debug: Output Hyperfine matrix================================================
  //printMatrix(hyperfine, "Hyperfine", dim);





  
  
  
  //Compute eZeeman============================================================  
  gsl_matrix_complex* eZeeman = gsl_matrix_complex_calloc(dim,dim);
  gsl_vector_complex* Bstaticg= gsl_vector_complex_calloc(3); //first multiply the g tensor by the field
  gsl_blas_zgemv(CblasTrans, gsl_complex_rect(1.0,0), gtensor, Bstatic, gsl_complex_rect(0.0,0.0), Bstaticg);
  //depending on the convention, i might have to tranpose the gtensor here 
  //cout << "Bstaticg: \n";
  //for(int i=0; i<3; i++)
  //cout << GSL_REAL(gsl_vector_complex_get(Bstaticg, i)) << '\t';
  //cout << endl;
  gsl_complex Bgs;
  for(int i=0; i<dim; i++)
    {
      //cout << states[i] << " |";
      for(int j=0; j<dim; j++)
	{
	  a= states[i][nprotons];  //nprotons is always the index of the electron spin
	  b= states[j][nprotons];
	  //cout << '\t' << a << b;
	  flag=0;
	  for(int k=0; k<nprotons; k++)
	    {
	      if(states[i][k]!=states[j][k]) 
		{
		  flag =1;
		  //break;
		}
	    }
	  if(flag==1) continue;
	  gsl_vector_complex_set(s, 0, gsl_matrix_complex_get(half_xv,a,b));
	  gsl_vector_complex_set(s, 1, gsl_matrix_complex_get(half_yv,a,b));
	  gsl_vector_complex_set(s, 2, gsl_matrix_complex_get(half_zv,a,b));
	  gsl_blas_zdotu(Bstaticg, s, &Bgs);
	  gsl_matrix_complex_set(eZeeman,i,j,gsl_complex_add(gsl_matrix_complex_get(eZeeman,i,j),Bgs));
	}
      //cout << endl;
    }
  //a=0; b=0; flag=0;
  gsl_matrix_complex_scale(eZeeman, gsl_complex_rect(Bohrm,0));
  gsl_matrix_complex_add(Ham, eZeeman);
  //Debug: Output eZeeman matrix================================================
  //cout.precision(2);
//   printMatrix(eZeeman, "eZeeman", dim);
  //printMatrix(Ham, "Ham", dim);

  // have to check if g tensor is correct or transposed??
  // same for a tensor!!

  //check why with two atoms orientation makes a difference
  //compare the values in the hyperfine matrix, if they do not change with orientation then the probelm is somewhere in the 
  //							   zeeman matrices, or in numerical accuracy
  





  



  //gsl_matrix_complex* quadrupolar = gsl_matrix_complex_calloc(dim,dim);




  




  //Diagonalize the total Hamiltonian matrix===================================
  gsl_eigen_hermv_workspace* w = gsl_eigen_hermv_alloc(dim);
  gsl_vector* eigenvalues = gsl_vector_alloc(dim);
  gsl_matrix_complex* eigenvectors = gsl_matrix_complex_alloc(dim,dim);
  gsl_eigen_hermv(Ham, eigenvalues, eigenvectors,  w);
//   gsl_sort_vector(eigenvalues);
  
  //int ntransitions = (dim*(dim-1))/2;   //pairwise energy differences
  //gsl_vector* transitions = gsl_vector_alloc(ntransitions);  //vector to store transition frequencies
  //double swap;
  //short int sortflag;
  //cout.precision(5);
  //printMatrix(eigenvectors, "eigen", dim);


  //Compute Transition moments=================================================
  // gsl_complex product = gsl_complex_rect(0,0);
  // for(int i =0; i<dim; i++)
  //   {
  //     product = gsl_complex_add(product,gsl_complex_mul
  // 				(gsl_complex_conjugate(gsl_matrix_complex_get(eigenvectors,0,i))
  // 				 ,gsl_matrix_complex_get(eigenvectors,1,i)));
  //   }

  // cout << GSL_REAL(product) << endl;
  // cout << GSL_IMAG(product) << endl;
  // cout << gsl_complex_abs(product) << endl;

  //  gsl_vector_complex* Bperturb = gsl_vector_complex_calloc(3);
  //  {
  // double B1[3]={1,0,0}; //field direction
  // for(int i=0;i<3; i++)
  //  gsl_vector_complex_set(Bstatic, i, gsl_complex_rect(B1[i],0));
  // }
  
  
  gsl_complex xmoment;
  gsl_matrix_complex* moments = gsl_matrix_complex_calloc(dim,dim);
  for(int i=0; i<dim; i++)
    {
      for(int j=0; j<dim; j++)
  	{
	  for(int k=0; k<nprotons+1;k++)    //nprotons is always the index of the electronic spin state
	    {
	      flag=0;
	      for(int l=0;l<nprotons+1;l++)
		{
		  if(l==k) continue;
		  if(states[i][l]!=states[j][l])
		    {
		      flag =1;
		      break;
		    }
		}
	      if(flag==1) continue;
	      a= states[i][k];  //spin state of state k in row i
	      b= states[j][k];  //spin state of state k in column j
	      xmoment = gsl_matrix_complex_get(half_xv,a,b);
	      //gsl_vector_complex_set(I, 0, gsl_matrix_complex_get(half_xv,a,b)); //I here is either I or S
	      //gsl_vector_complex_set(I, 1, gsl_matrix_complex_get(half_yv,a,b));
	      //gsl_vector_complex_set(I, 2, gsl_matrix_complex_get(half_zv,a,b));
	      //gsl_blas_zdotu(Bstatic, I, &BI);
	      //BI = gsl_complex_mul_real(BI, -1.0);
	      if(k!=nprotons)
		xmoment = gsl_complex_mul_real(xmoment, -1.0*g_1H*NUC_MAGNETON);
	      if(k==nprotons)
		xmoment = gsl_complex_mul_real(xmoment, 2.023*Bohrm);
//               cout << i << '\t' << j << '\t' << k << '\t' << GSL_REAL(xmoment) << '\t' << GSL_IMAG(xmoment) << endl;
	      gsl_matrix_complex_set(moments,i,j,gsl_complex_add(gsl_matrix_complex_get(moments,i,j),xmoment));
	    }
//               cout << i << '\t' << j << '\t' << "FINAL:" << '\t' << GSL_REAL(gsl_matrix_complex_get(moments,i,j)) << '\t' << GSL_IMAG(gsl_matrix_complex_get(moments,i,j)) << endl;
	}
    }
  //cout.precision(5);
  //Transform the moments matrix==============================================
  gsl_matrix_complex* transformed_moments = gsl_matrix_complex_calloc(dim,dim);
  gsl_matrix_complex* intermediate = gsl_matrix_complex_calloc(dim,dim);
  gsl_blas_zgemm(CblasConjTrans, CblasNoTrans, gsl_complex_rect(1,0),
		 eigenvectors, moments, gsl_complex_rect(1,0), intermediate);
  gsl_blas_zgemm(CblasNoTrans, CblasNoTrans, gsl_complex_rect(1,0),
		 intermediate, eigenvectors, gsl_complex_rect(1,0), transformed_moments);
  //printMatrix(transformed_moments, "trans", dim);
  //Debug: output moments=====================================================
  //printMatrix(moments, "Moments", dim);
  gsl_matrix* moments_sq = gsl_matrix_calloc(dim,dim);
  int counter =0;
  for(int i=0;i<dim;i++)
    {
      for(int j=0;j<dim;j++)
	{
	  //gsl_vector_set(moments_sq, counter, gsl_complex_abs2(gsl_matrix_complex_get(moments, i, j)));
	  //counter++;
	  gsl_matrix_set(moments_sq, i, j, gsl_complex_abs2(gsl_matrix_complex_get(moments, i, j)));
	}
    }
  double max_moment = gsl_matrix_max(moments_sq);
  gsl_matrix_scale(moments_sq, 1.0/max_moment);
  
  //for(int i=0;i<dim*(dim-1)/2;i++)
  //cout << i+1 << ": " << gsl_vector_get(moments_sq, i) << endl;
  //printMatrix(moments_sq, "moments_sq", dim);
  
  
  
 
  // counter =0;
  // for(int i=0; i<dim; i++)
  //   {
  //     for(int j=i+1; j<dim; j++)
  // 	{
  // 	  gsl_vector_set(transitions,counter, abs(gsl_vector_get(eigenvalues, i)-gsl_vector_get(eigenvalues, j)));
  // 	  counter++;
  // 	}
  //   }
  
  
  
  
  
  //Output=====================================================================
  //cout << "\nEigenvalues: " << endl;
  //cout << "============\n";
  //cout.precision(5);
  //for(int i=0; i<dim; i++)
    //cout << scientific << gsl_vector_get(eigenvalues, i) << endl;
  //cout << "\nTransition Frequencies (GHz): " << endl;
  //cout << "================================\n";
  //cout.precision(5);
  gsl_matrix *transition_frequencies = gsl_matrix_calloc(dim,dim);
  for(int i=0;i<dim;i++)
    {
      for(int j=0;j<dim;j++)
	{
	  gsl_matrix_set(transition_frequencies, i, j,
			 abs(gsl_vector_get(eigenvalues, i)-gsl_vector_get(eigenvalues, j)));
	}
    }
  gsl_matrix_scale(transition_frequencies, 1.0/h/1.0E9);
//   printMatrix(transition_frequencies, "frequencies", dim);
  printVector(eigenvalues, "eigenvalues", dim);

  

  cout << "\n\nAllowed Transitions: \t\t\tB= " << fixed << B << endl;
  cout << "-------------------- " << endl;
  cout << "Transition\tFrequency (GHz)\t\tProbability" << endl;
  cout << "---------------------------------------------------\n";
  int transitions = 0;
  for(int i=0;i<dim;i++)
    {
      for(int j=i+1;j<dim;j++)
	{
	  if(gsl_matrix_get(moments_sq, i,j) >1.0E-6)
	    {
	      
	      cout << fixed << i+1 << " -> " << j+1 << "\t\t";
	      cout.precision(5);cout.width(10);
	      cout << right << gsl_matrix_get(transition_frequencies,i,j) << "\t\t";
	      cout.precision(8);
	      cout << gsl_matrix_get(moments_sq, i, j)
		   << endl;//<< "\t\t" << 1e24*GSL_REAL(gsl_matrix_complex_get(transformed_moments, i, j)) << endl;
              ++transitions;
	    }
	}
    }
  cout << transitions << " transitions in total" << endl;
  
  // for(int i=0; i<ntransitions; i++)
  //   {
  //     if(gsl_vector_get(transitions,i)>0)
  //       cout << i+1 << ": " << fixed << gsl_vector_get(transitions,i)/h/1.0E9 << " GHz" << endl;  
  //   }
  // cout << endl;
  // double difference=0;
  // for(int i=0; i<dim; i++)
  //   {
  //     for(int j=dim-1; j>i; j--)
  // 	{
  // 	  difference = abs(gsl_vector_get(eigenvalues, i)-gsl_vector_get(eigenvalues, j));
  // 	  if(difference>0)
  // 	    cout << i+1 << "->" << j+1 << '\t' << fixed << difference/h/1.0E9 << " GHz" << endl;
  // 	}
  //   }
  cout.clear();
  //consistency check
  // cout << "\nA0: " << (gsl_vector_get(eigenvalues, 0)-gsl_vector_get(eigenvalues, 1) 
  // 		       - gsl_vector_get(eigenvalues, 2) + gsl_vector_get(eigenvalues, 3))/h/1E6
  //      << " MHz " <<  endl;
  // cout << "A0: " << (gsl_vector_get(eigenvalues, 0)-gsl_vector_get(eigenvalues, 1) 
  // 		     - gsl_vector_get(eigenvalues, 2) + gsl_vector_get(eigenvalues, 3))/h/28.02495/1E6
  //      << " mT " <<endl;
  // cout << "A0: " << 1420*1E9*h/2.0022838/Bohrm       << " mT " <<endl;


  //Free memory================================================================
  gsl_matrix_complex_free(hyperfine);
  gsl_vector_complex_free(atensor_I);
  gsl_vector_complex_free(I);
  gsl_matrix_complex_free(nZeeman);  
  gsl_vector_complex_free(Bstaticg);
  gsl_vector_complex_free(s);
  gsl_matrix_complex_free(gtensor);
  gsl_vector_complex_free(Bstatic);
  gsl_matrix_complex_free(half_xv);
  gsl_matrix_complex_free(half_yv);
  gsl_matrix_complex_free(half_zv);
  gsl_vector_free(eigenvalues);
  gsl_eigen_hermv_free(w);
  gsl_matrix_complex_free(Ham);
  gsl_matrix_complex_free(eZeeman);
  gsl_matrix_complex_free(eigenvectors);
  //gsl_matrix_complex_free(quadrupolar);
  gsl_matrix_complex_free(atensor);
  gsl_matrix_free(moments_sq);
  gsl_matrix_complex_free(moments);
  gsl_matrix_complex_free(transformed_moments);
  gsl_matrix_complex_free(intermediate);
  gsl_matrix_free(transition_frequencies);
}


void printMatrix(const gsl_matrix_complex *matrix, const string name, const int dim)
{
  cout << "\n\n----------------------------------------------------\n";
  cout << name << ":REAL: \n";
  for(int i=0; i<dim; i++)
    {
      for(int j=0; j<dim; j++)
  	{
  	  cout << scientific << GSL_REAL(gsl_matrix_complex_get(matrix, i, j)) << '\t';
  	}
      cout << endl;
    }
  cout << name << ":IMAG: \n";
  for(int i=0; i<dim; i++)
    {
      for(int j=0; j<dim; j++)
  	{
  	  cout << GSL_IMAG(gsl_matrix_complex_get(matrix, i, j)) << '\t';
  	}
      cout << endl;
    }
  cout << "----------------------------------------------------\n";
  
}


void printMatrix(const gsl_matrix *matrix, const string name, const int dim)
{
  cout << "\n\n----------------------------------------------------\n";
  cout << name << ":REAL: \n";
  for(int i=0; i<dim; i++)
    {
      for(int j=0; j<dim; j++)
  	{
  	  cout << scientific << gsl_matrix_get(matrix, i, j) << "\t";
  	}
      cout << endl;
    }
  cout << "----------------------------------------------------\n";  
}

void printVector(const gsl_vector_complex *vector, const string name, const int dim)
{
  cout << "\n\n----------------------------------------------------\n";
  cout << name << ": \n";
  for(int i=0; i<dim; i++) {
    cout << scientific << '(' << GSL_REAL(gsl_vector_complex_get(vector, i)) << ',' << GSL_IMAG(gsl_vector_complex_get(vector, i)) << ')' << '\t';
  }
  cout << endl;
  cout << "----------------------------------------------------\n";
  
}

void printVector(const gsl_vector *vector, const string name, const int dim)
{
  cout << "\n\n----------------------------------------------------\n";
  cout << name << ": \n";
  for(int i=0; i<dim; i++) {
    cout << scientific << gsl_vector_get(vector, i) << '\t';
  }
  cout << endl;
  cout << "----------------------------------------------------\n";
  
}