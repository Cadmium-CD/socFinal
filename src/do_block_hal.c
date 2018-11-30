#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <sys/mman.h>
#include <fcntl.h>
#include "do_block_hal.h"
#define MAP_SIZE 524288UL //consider block size 100 27 1000
#define MAP_MASK (MAP_SIZE -1)
  volatile void *base;
  volatile int *flag_address;  //store special address to control run/stop
  //volatile int *c_address;
  unsigned long int addr0;


int do_block_hal (int lda,int ldb, int ldc, int M, int N, int K, int* A, int* B, int* C)
{

  int finish_flag; //flag 1 means computation done
  addr0 = 0x40000000ul;
  finish_flag = 0;
  //predefine address


  //Memory map
  int fd = open("/dev/mem", O_RDWR|O_SYNC);
  if(!fd)
    {
      printf("Unable to open /dev/mem.  Ensure it exists (major=1, minor=1)\n");
      return -1;
    }
  Drv_initialize(fd);

  //Input data
  Drv_inputdata1(lda);
  Drv_inputdata1(ldb);
  Drv_inputdata1(ldc);
  Drv_inputdata1(M);
  Drv_inputdata1(N);
  Drv_inputdata1(K);
  Drv_inputdata2(A, M, K, lda);
  Drv_inputdata2(B, K, N, ldb);
  Drv_inputdata2(C, M, N, ldc);

  //Run accelerator
  Drv_startfunction();

  //Wait until finish
  finish_flag = Drv_waitfunction();

  //outdata data and ummap
  //if (1)
  if (finish_flag)
  {
    Drv_outputdata(C, M, N, ldc);
    Drv_stopfunction(fd);
  }

  //return finish_flag;

}

void Drv_initialize(int fd)
{
  //Open memory as a file
  base = (unsigned int *)mmap(NULL, MAP_SIZE, PROT_READ|PROT_WRITE, MAP_SHARED, fd, addr0 & ~MAP_MASK);
    if((base == MAP_FAILED))
    {
      printf("mapping failed\n");
      fflush(stdout);
    }

  //synchronize the virtual memory and physical memory before ummap
  msync((void*)base, MAP_SIZE, MS_SYNC);
  flag_address = base;
}

void Drv_inputdata1(int data1)
{
   base = base + 4;
   *((int*) base) = data1;
}

void Drv_inputdata2(int *data2, int height, int width, int ld)
{
  int i,j;
  for (i = 0; i<height; ++i)
    for (j = 0; j<width; ++j)
    {
        base = base + 4;
	//printf("%d \n", data2[i*ld + j]);
        *((int*) base) = data2[i*ld + j];
    }
}

void Drv_startfunction()
{
  *(int*) flag_address = 12;
}

int Drv_waitfunction()
{
  volatile int result;
  do{
    //printf("waiting...\n");
    result = *((int*) flag_address);
  } while (result==12); //continue waiting until result is 0
  printf("result = %d \n",result);
  return 1;
}

void Drv_outputdata(int *data, int height, int width, int ld)
{
  //int* temp_data;
  //temp_data = (int*) malloc(height*width*sizeof(int));
  int i,j;
  //printf("base = %d \n",*((int*)(base+4-(height*width*4))));
  /*memcpy(temp_data, (int*)(base+4-(height*width*4)),height*width*sizeof(int)); //copy data back
  // cut data to subblock
  for (i = 0; i<height; ++i)
    for (j = 0; j<width; ++j)
    {
      data[i*ld + j]=temp_data[i*width+j];
      //printf("temp_data = %d \n",temp_data[i*width+j]);
    }
  printf("C[0] = %d \n",data[0]);
  free(temp_data);*/
  //c_address = base+4-(height*width*4);
  for (i = 0; i<height; ++i)
    for (j = 0; j<width; ++j)
    {
       data[i*ld + j] = *((int*) (base+4-(height*width*4)+i*width*4+j*4) );
    }
  //printf("C[0] = %d \n",data[0]);
  *(int*) flag_address = 3;
}

void Drv_stopfunction(int fd)
{
    munmap((void*)flag_address, MAP_SIZE);
    close(fd);
}
