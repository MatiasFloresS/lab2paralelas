#include <stdio.h>
#include <pmmintrin.h>
#include <emmintrin.h>
#include <fcntl.h>
#include <smmintrin.h>
#include <string.h>
#include <getopt.h>
#include <omp.h>
#include <math.h>
#include <iostream>
#include <unistd.h>
#include <limits.h>
using namespace std;

#define n 16


struct my_struct{
	__m128 reg,reg2;
};

struct minheap{
    float a; 
    int b,c; 
};

void intercambio_elementos(minheap *elem1, minheap *elem2){
    minheap aux = *elem1;
    *elem1 = *elem2;
    *elem2 = aux;
}


class MinHeap{
    minheap *arrayheap;
    int sizeheap; 

public:
    MinHeap(minheap array[], int size);
    void MinHeapify(int );

    int left(int i){
    	return (2*i + 1);
    }

    int right(int i){
    	return (2*i + 2);
    }
 
    minheap raiz(){
    	return arrayheap[0];
    }
 
    void cambiar(minheap aux){
    	arrayheap[0] = aux;
    	MinHeapify(0);
    }
};

float *merge(float **arr, int filas, int columnas){
    float *salida = new float[filas*columnas];

    minheap *arrayheap = new minheap[filas];
    for (int i = 0; i < filas; i++){
        arrayheap[i].a = arr[i][0]; 
        arrayheap[i].b = i;  
        arrayheap[i].c = 1;  
    }

    MinHeap heap(arrayheap, filas); 
 	for (int aux2 = 0; aux2 < columnas*filas; aux2++){
        minheap aux = heap.raiz();
        salida[aux2] = aux.a;
 
        if (aux.c < n){
            aux.a = arr[aux.b][aux.c];
            aux.c += 1;
        }
        else {
        	aux.a =  INT_MAX;
 		}
        heap.cambiar(aux);
    }
 
    return salida;
}

MinHeap::MinHeap(minheap array[], int size){
    sizeheap = size;
    arrayheap = array; 
    int i = (sizeheap - 1)/2;
    while (i >= 0){
        MinHeapify(i);
        i--;
    }
}

void MinHeap::MinHeapify(int i){
	int l,r,menor;
	l=left(i);
	r=right(i);
	menor = i;
    if (l < sizeheap && arrayheap[l].a < arrayheap[i].a){
        menor = l;
    }
    if (r < sizeheap && arrayheap[r].a < arrayheap[menor].a){
        menor = r;
    }
    if (menor != i){
        intercambio_elementos(&arrayheap[i], &arrayheap[menor]);
        MinHeapify(menor);
    }
}

struct my_struct bmn(__m128 a, __m128 b) {

 	__m128 auxiliar,auxiliar2;
 	float array[4] __attribute__((aligned(n))) ;
 	float array1[4] __attribute__((aligned(n))) ;

	a = _mm_shuffle_ps(a,a,_MM_SHUFFLE(3,1,2,0));
	b = _mm_shuffle_ps(b,b,_MM_SHUFFLE(3,1,2,0));

	auxiliar2=_mm_min_ps(a,b);
	auxiliar=_mm_max_ps(a,b);

    a = _mm_shuffle_ps(auxiliar,auxiliar,_MM_SHUFFLE(3,1,2,0));
	b = _mm_shuffle_ps(auxiliar2,auxiliar2,_MM_SHUFFLE(3,1,2,0)); 

	auxiliar = _mm_shuffle_ps(a,b, _MM_SHUFFLE(1,0,1,0));
	auxiliar2 = _mm_shuffle_ps(a,b, _MM_SHUFFLE(3,2,3,2));
	a=auxiliar;
	b=auxiliar2;

	a =_mm_shuffle_ps(a,a,_MM_SHUFFLE(3,1,2,0));
	b =_mm_shuffle_ps(b,b,_MM_SHUFFLE(3,1,2,0));

	auxiliar2 = _mm_min_ps(a,b);
	auxiliar = _mm_max_ps(a,b);

	a = _mm_shuffle_ps(auxiliar,auxiliar2,_MM_SHUFFLE(1,0,1,0));
	b = _mm_shuffle_ps(auxiliar,auxiliar2,_MM_SHUFFLE(3,2,3,2));

	a = _mm_shuffle_ps(a,a,_MM_SHUFFLE(3,1,2,0));
	b = _mm_shuffle_ps(b,b,_MM_SHUFFLE(3,1,2,0));

	auxiliar2 = _mm_min_ps(a,b);
	auxiliar = _mm_max_ps(a,b);

	a = _mm_shuffle_ps(auxiliar,auxiliar2,_MM_SHUFFLE(1,0,1,0));
	b = _mm_shuffle_ps(auxiliar,auxiliar2,_MM_SHUFFLE(3,2,3,2));

	auxiliar = _mm_shuffle_ps(a,a,_MM_SHUFFLE(0,2,1,3));
	auxiliar2 = _mm_shuffle_ps(b,b,_MM_SHUFFLE(0,2,1,3));

	// merge simd

	_mm_store_ps(array, auxiliar2);
    _mm_store_ps(array1, auxiliar);    


    struct my_struct r;
    r.reg = auxiliar2;
    r.reg2 = auxiliar;
	
    return r;
}

void MostrarNumeros(float arreglo[], int size){
	for (int i=0; i < size; i++){
       printf("%f\n",arreglo[i]);
    }

}

void funcion (float z[],float **v,int var1){
	float array[4]__attribute__((aligned(n)))= {z[0],z[1],z[2],z[3]};
	float array1[4] __attribute__((aligned(n)))= {z[4],z[5],z[6],z[7]};
	float array2[4] __attribute__((aligned(n)))= {z[8],z[9],z[10],z[11]};
	float array3[4] __attribute__((aligned(n)))= {z[12],z[13],z[14],z[15]};
    __m128 a,b,c,d,auxiliar,auxiliar2;
    __m128 temp, temp2, temp3, tempz, tempz2, fila1, fila2, fila3, fila4;
	__m128 A,B,C,D;
	__m128 aux1,aux2,Ef,Ff;
	__m128 or2,or3,or4,or5;
	int i;
	struct my_struct reg;
	struct my_struct reg2;
	struct my_struct reg3;
	struct my_struct reg4;
	struct my_struct reg5;
  

    a= _mm_load_ps(array);
    b= _mm_load_ps(array1);
    c= _mm_load_ps(array2);
    d= _mm_load_ps(array3);

	// red min max
	auxiliar=_mm_min_ps(a,c);
	auxiliar2=_mm_max_ps(a,c);
	a= auxiliar;
	c=auxiliar2;

	auxiliar=_mm_min_ps(b,d);
	auxiliar2=_mm_max_ps(b,d);
	b=auxiliar;
	d=auxiliar2;

	auxiliar=_mm_min_ps(a,b);
	auxiliar2=_mm_max_ps(a,b);
	a=auxiliar;
	b=auxiliar2;

	auxiliar=_mm_min_ps(c,d);
	auxiliar2=_mm_max_ps(c,d);
	c=auxiliar;
	d=auxiliar2;

	auxiliar=_mm_min_ps(b,c);
	auxiliar2=_mm_max_ps(b,c);
	b=auxiliar;
	c=auxiliar2;

	// bitonica


	temp = _mm_shuffle_ps(a,b, _MM_SHUFFLE(1,0,1,0));
	fila2 = _mm_shuffle_ps(c, d, _MM_SHUFFLE(1,0,1,0));

	fila1 = _mm_shuffle_ps(temp,fila2,_MM_SHUFFLE(2,0,2,0));
	fila2 = _mm_shuffle_ps(temp,fila2,_MM_SHUFFLE(3,1,3,1));

	temp = _mm_shuffle_ps(a,b, _MM_SHUFFLE(3,2,3,2));
	temp2 = fila1;
	temp3 = fila2;
	fila2 = _mm_shuffle_ps(c,d, _MM_SHUFFLE(3,2,3,2));
     
	fila1= _mm_shuffle_ps(temp,fila2,_MM_SHUFFLE(2,0,2,0));
	fila2 = _mm_shuffle_ps(temp,fila2,_MM_SHUFFLE(3,1,3,1));
	a = temp2;
	b = _mm_shuffle_ps(temp3,temp3,_MM_SHUFFLE(0,1,2,3));
	
	c = fila1;
	d = _mm_shuffle_ps(fila2,fila2,_MM_SHUFFLE(0,1,2,3));
	
	reg=bmn(a,b);
	reg2=bmn(c,d);
	A = reg.reg;
	B = reg.reg2;
    C = reg2.reg;
    D = reg2.reg2;

    C = _mm_shuffle_ps(C,C,_MM_SHUFFLE(0,1,2,3));
    reg3 = bmn(A,C);
	_mm_store_ps(array, reg3.reg);


	Ef = _mm_shuffle_ps(B,B,_MM_SHUFFLE(0,0,0,0));
	Ff = _mm_shuffle_ps(D,D,_MM_SHUFFLE(0,0,0,0));
	aux1 = _mm_cmpgt_ps(Ef, Ff); 
	aux2 = _mm_cmplt_ps(Ef,Ff);
	or2 = _mm_blendv_ps(B, D, aux1);
	or3 = _mm_blendv_ps(B,D,aux2);

	or2 = _mm_shuffle_ps(or2,or2,_MM_SHUFFLE(0,1,2,3));
	reg4=bmn(reg3.reg2,or2);
    _mm_store_ps(array1, reg4.reg);


    or3 = _mm_shuffle_ps(or3,or3,_MM_SHUFFLE(0,1,2,3));
	reg5=bmn(reg4.reg2,or3);
	_mm_store_ps(array2, reg5.reg);

    _mm_store_ps(array3, reg5.reg2);
    
    int p;
	
	#pragma omp critic(v)
	for (p = 0; p < 4; p++)
	{
		v[var1][p] = array[p];
		v[var1][p+4] = array1[p];
		v[var1][p+8] = array2[p];
		v[var1][p+12] = array3[p];
	}

    return;
}



int main(int argc, char *argv[]){
	int numero=0;
	char archivo[128];
	char archivo_salida[128];
	int N;
	int debug;
	int c;
	int L;
	int index=0;
	int tid;
	
	static struct option long_options[] = {
		{"Input",1,NULL,'i'},
		{"Output",1,NULL,'o'},
		{"Numbers",1,NULL,'N'},
		{"Debug",1,NULL,'d'},
		{"NiveldeRecursividad",1,NULL,'L'},
		{0,0,0,0}
	};

	while(1){
		c = getopt_long(argc, argv, "i:o:N:d:L:",long_options,&index);

		if(c==-1){
			break;
		}

		switch(c){
			case 'i':
				strcpy(archivo,optarg);
				break;
			case 'o':
				strcpy(archivo_salida,optarg);
				break;
			case 'N':
				N = atoi(optarg);
				break;
			case 'd':
				debug = atoi(optarg);
				break;
			case 'L':
				L = atoi(optarg);
				break;

			default:
				
				break;
		}
	}
	int cant_hebras = pow (2,L) ;
	omp_set_num_threads(cant_hebras);
	int fd = open(archivo, O_RDONLY); 
	if(fd == -1){ 
		printf("%s\n", "No existe el fichero");
		
		return 0; 
	} 
	float *buff;
	size_t size = sizeof(float) * N; 
	posix_memalign((void **) &buff, n, size);
	
	int i,q=0,w=0,e=0,r=0,var=0;
	int cont =0;
	
	float *arreglo;
	float matriz[N/n][n]__attribute__((aligned(n)));
	
	for( i = 0; i < N; i++){
		 	read(fd, buff, size);
		 	matriz[var][q]= *buff;
		 	q++;
		 	buff++;

		 	if(q % n == 0){
		 		var++;
		 		q=0;
		 	}
		}
	close(fd); 
    int var1;
	int aux = N/n;
	int numero2= aux/cant_hebras;
	float matriz2[aux][n]__attribute__((aligned(n)));
	float **matriz3 = (float **) malloc(sizeof(float *)*aux);
	int s;
    for(s=0; s<aux; s++){
    /* Allocate array, store pointer  */
    matriz3[s] = (float *) malloc(sizeof(float)*n);
	}
	#pragma omp parallel private(tid) shared(matriz2)
	{
		#pragma omp for schedule(static, numero2) ordered
		for (var1= 0; var1 < N/n; var1++){
		 	funcion(matriz[var1],matriz3,var1);
		}

		#pragma omp barrier
	}

	float tamanio = sizeof(matriz2)/sizeof(matriz2[0]);
	printf("%f\n",tamanio );
	int filas = aux;
	int columnas=16;
	float *salidaFinal = merge(matriz3,filas,columnas);	

	
	if(debug == 1){
		MostrarNumeros(salidaFinal, n*filas);
	}	

	int filesdesc2 = open(archivo_salida, O_WRONLY | O_CREAT, 0644);
	float nwrite;
	int buffHelp = n*tamanio;
	int wHelp = 0;

	while(buffHelp > 0){
		nwrite = salidaFinal[wHelp];
		write(filesdesc2, &nwrite, sizeof(float));
		wHelp++;
		buffHelp--;
	}
	close(filesdesc2);	
		
}