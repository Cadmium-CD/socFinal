
#include <assert.h>
#include <ap_axi_sdata.h>

//#define DB_DEBUG


//#define MCR_SIZE 1024
#define NUM_PARA 3
typedef ap_axiu<32,4,5,5> AXI_VAL;

// function prototypes
/*void standalone_mmult (float A[32][32], float B[32][32], float C[32][32]);*/
void HLS_accel (AXI_VAL in_stream[27*50 + 27*500 +50*500 + NUM_PARA], AXI_VAL out_stream[50*500]);



/* ****************************** C++ TEMPLATES ************************************** */

// reference function
template <typename T,int DIMA,int DIMB,int DIMC>
void matrix_multiply_ref(T a[DIMA], T b[DIMB], T c[DIMC], T out[DIMC],T M, T N, T K)
{

	// matrix multiplication of a A*B matrix
	/*for (int ia = 0; ia < DIM; ++ia)
		for (int ib = 0; ib < DIM; ++ib)
		{

			float sum = 0;

			for (int id = 0; id < DIM; ++id)

				sum += a[ia][id] * b[id][ib];

			out[ia][ib] = sum;
		}*/
	int i,j,k;

	L1: for (i = 0; i < M; ++i) {
        L2: for (k = 0; k < K; ++k) {
    		int A_PART = a[i*K + k];
            L3: for (j = 0; j < N; ++j) {
				out[i*N + j] = c[i*N + j] + A_PART * b[k*N + j];
			}
		}
	}

		return;
}

// --------------------------------------------------------
// function to be accelerated in HW
template <typename T,int DIMA,int DIMB,int DIMC>
void mmult_hw(T a[DIMA], T b[DIMB], T c[DIMC], T out[DIMC],T M, T N, T K)
{

	//int const FACTOR = DIM/2;
	#pragma HLS INLINE
	//#pragma HLS array_partition variable=a block factor=FACTOR dim=2
	//#pragma HLS array_partition variable=b block factor=FACTOR dim=1

	// matrix multiplication of a A*B matrix
	/*L1:for (int ia = 0; ia < DIM; ++ia)
		L2:for (int ib = 0; ib < DIM; ++ib)
		{
			//#pragma HLS PIPELINE II=1
			T sum = 0;
			L3:for (int id = 0; id < DIM; ++id)
				sum += a[ia][id] * b[id][ib];
			out[ia][ib] = sum;
		}*/
	//gemm
	//printf("k = %d",K);
	int i,j,k;

	L1: for (i = 0; i < M; ++i) {
        L2: for (k = 0; k < K; ++k) {
    		int A_PART = a[i*K + k];
            L3: for (j = 0; j < N; ++j) {
				out[i*N + j] = c[i*N + j] + A_PART * b[k*N + j];
				//printf("out = %d c = %d a = %d b = %d \r\n",out[i*N + j],c[i*N + j],a[i*K + k],b[k*N + j]);
			}
		}
	}

	//printf("out = %d c = %d a = %d b = %d",out[0],c[0],a[0],b[0]);
		return;
}

// --------------------------------------------------------
// functions to insert and extract elements from an axi stream
// includes conversion to correct data type

template <typename T, int U, int TI, int TD>
T pop_stream(ap_axiu <sizeof(T)*8,U,TI,TD> const &e)
{
#pragma HLS INLINE

	assert(sizeof(T) == sizeof(int));
	/*union
	{
		int ival;
		T oval;
	} converter;
	converter.ival = e.data;
	T ret = converter.oval;*/
	T ret = e.data;

	volatile ap_uint<sizeof(T)> strb = e.strb;
	volatile ap_uint<sizeof(T)> keep = e.keep;
	volatile ap_uint<U> user = e.user;
	volatile ap_uint<1> last = e.last;
	volatile ap_uint<TI> id = e.id;
	volatile ap_uint<TD> dest = e.dest;

	return ret;
}

template <typename T, int U, int TI, int TD> ap_axiu <sizeof(T)*8,U,TI,TD> push_stream(T const &v, bool last = false)
{
#pragma HLS INLINE
	ap_axiu<sizeof(T)*8,U,TI,TD> e;

	assert(sizeof(T) == sizeof(int));
	/*union
	{
		int oval;
		T ival;
	} converter;
	converter.ival = v;
	e.data = converter.oval;*/
	e.data = v;

	// set it to sizeof(T) ones
	e.strb = -1;
	e.keep = 15; //e.strb;
	e.user = 0;
	e.last = last ? 1 : 0;
	e.id = 0;
	e.dest = 0;
	return e;
}

// --------------------------------------------------------------------
// function to be accelerated in HW wrapped with AXI4-Stream interface

template <typename T, int DIMA, int DIMB, int DIMC,  int U, int TI, int TD>
void wrapped_mmult_hw (
	AXI_VAL in_stream[NUM_PARA + DIMA + DIMB + DIMC],
	AXI_VAL out_stream[DIMC])
{

#pragma HLS INLINE

	//T lda;
	//T ldb;
	//T ldc;
	T M;
	T N;
	T K;
	T a[DIMA];
	T b[DIMB];
	T c[DIMC];
	T out[DIMC];

	assert(sizeof(T)*8 == 32);

	//lda = pop_stream<T,U,TI,TD>(in_stream[0]);
	//ldb = pop_stream<T,U,TI,TD>(in_stream[1]);
	//ldc = pop_stream<T,U,TI,TD>(in_stream[2]);
	M = pop_stream<T,U,TI,TD>(in_stream[0]);
	N = pop_stream<T,U,TI,TD>(in_stream[1]);
	K = pop_stream<T,U,TI,TD>(in_stream[2]);
	//printf("M = %d N = %d K = %d",M,N,K);
	// stream in first matrix
	for(int i=0; i<DIMA; i++){
#pragma HLS PIPELINE II=1
			a[i] = pop_stream<T,U,TI,TD>(in_stream[i + NUM_PARA]);
		}

		// stream in second matrix
		for(int i=0; i<DIMB; i++){
#pragma HLS PIPELINE II=1
				b[i] = pop_stream<T,U,TI,TD>(in_stream[i + NUM_PARA + DIMA]);
				//printf("b = %d \r\n",b[i]);
			}

			for(int i=0; i<DIMC; i++){
#pragma HLS PIPELINE II=1
				c[i] = pop_stream<T,U,TI,TD>(in_stream[i + NUM_PARA + DIMA + DIMB]);
				//printf("c = %d \r\n",c[i]);
			}
			//printf("begin");
			// do HW multiplication
			mmult_hw<T,DIMA,DIMB,DIMC>(a,b,c,out,M,N,K);

			// stream out result matrix
			for(int i=0; i<DIMC; i++){
					#pragma HLS PIPELINE II=1
					out_stream[i] = push_stream<T,U,TI,TD>(out[i],i == (DIMC-1));
					//printf("out = %d \r\n",out[i]);
				}
				return;

}



// test the functions
template <typename T, int DIMA, int DIMB, int DIMC,  int U, int TI, int TD>
int test_matrix_mult(void)
{
	int i,j, err;
	T M = 50;
	T N = 500;
	T K = 27;
	T lda = K;
	T ldb = N;
	T ldc = N;
	T matOp1[DIMA];
	T matOp2[DIMB];
	T matOp3[DIMC];
	T matMult_sw[DIMC];
	T matMult_hw[DIMC];

	/** Matrix Initiation */
	for(i = 0; i<DIMA; i++){
			matOp1[i] = (int)(i);
	}
    i = 0;
	for(i = 0; i<DIMB; i++){
			matOp2[i] = (int)(i * 2);
	}
	i = 0;
	for(i = 0; i<DIMC; i++){
				matOp3[i] = 0;
		}
	i = 0;
	/** End of Initiation */


	printf("DEBUGGING AXI4 STREAMING DATA TYPES!\r\n");

	// prepare data for the DUT
	AXI_VAL inp_stream[NUM_PARA + DIMA + DIMB + DIMC];
	AXI_VAL out_stream[DIMC];

	assert(sizeof(T)*8 == 32);

	//inp_stream[0] = push_stream<T,U,TI,TD>(lda,0);
	//inp_stream[1] = push_stream<T,U,TI,TD>(ldb,0);
	//inp_stream[2] = push_stream<T,U,TI,TD>(ldc,0);
	inp_stream[0] = push_stream<T,U,TI,TD>(M,0);
	inp_stream[1] = push_stream<T,U,TI,TD>(N,0);
	inp_stream[2] = push_stream<T,U,TI,TD>(K,0);

	// stream in the first input  matrix
	for(int i=0; i<DIMA; i++){
			inp_stream[i + NUM_PARA] = push_stream<T,U,TI,TD>(matOp1[i], i == DIMA - 1);
			//printf("i = %d\r\n",i);
			//printf("a = %d\r\n",matOp1[i]);
		}
	    i = 0;
		// stream in the second input  matrix
		for(int i=0; i<DIMB; i++){
				inp_stream[i + NUM_PARA + DIMA] = push_stream<T,U,TI,TD>(matOp2[i],i == (DIMB-1));
			}
		    i = 0;
			for(int i=0; i<DIMC; i++){
				inp_stream[i + NUM_PARA + DIMA + DIMB] = push_stream<T,U,TI,TD>(matOp3[i],i == (DIMC-1));
			}

			i = 0;
			//call the DUT
			HLS_accel (inp_stream, out_stream);
			//wrapped_mmult_hw<T, DIMA, DIMB, DIMC,  U, TI, TD>(inp_stream, out_stream);

			// extract the output matrix from the out stream 
			for(int i=0; i<DIMC; i++){
					matMult_hw[i] = pop_stream<T,U,TI,TD>(out_stream[i]);
				}
			   i = 0;

	/* reference Matrix Multiplication */
	matrix_multiply_ref<T, DIMA, DIMB, DIMC>(matOp1, matOp2, matOp3, matMult_sw ,M, N, K);

	for (i = 0; (i<DIMC); i++){
				//printf("sw = %d hw = %d\r\n",matMult_sw[i],matMult_hw[i]);
		}
	/** Matrix comparison */
	err = 0;
	for (i = 0; (i<DIMC && !err); i++){
			if (matMult_sw[i] != matMult_hw[i]) 
				err++;
			//printf("sw = %d hw = %d",matMult_sw[i],matMult_hw[i]);
	}
	if (err == 0)
		printf("Matrixes identical ... Test successful!\r\n");
	else
		printf("Test failed!\r\n");

	return err;
}


