#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <stdlib.h>

//NUMERO BLOQUES = INTEGERMASALTO(NUMVERTICES/MAXTRHEADSXBLOCK)
//VARIABLES GLOBALES
#define NUMVERTICES 300      //MAXIMUM 300
#define MAXTRHEADSXBLOCK 32 //[1 - 32]

//ID gpudevice that is used
int gpudev = 0;

//Graph representation with
int *EG;    //Double array of edges |NUMVERTICES|
int *VG;    //Double array of vertices

int C;      //Current vertex INDEX
int NUMBEREDGES;

//MST edge list: Shows the path that is followed.
int *R1source;
int *R2destination;
int *R3weigth;

//Temporal arrays used for reduction results
int *T1weights;
int *T2indexes;

//------- FUNCIONES --------
void printDoubleArray(int *VX)
{
    int lengthArray = NUMVERTICES;

    for(int i = 0; i <lengthArray; i++)
    {
        printf("%i ", VX[i]);

    }//End for 2

    printf("\n");

    for(int i = lengthArray; i < (lengthArray*2); i++)
    {
         printf("%i ",VX[i]);

    }//End for 2

    printf("\n");

}//Fin funcion printDoubleArray
//------------------------
void printArrayRange(int *VX,int start,int end)
{

    for(int i = start; i <= end; i++)
    {
        printf("%2i ", VX[i]);

    }//End for 2

    printf("\n");

}//Fin funcion printArrayRange
//------------------------
/*
*Function that creates the vaules of the graph
*in the strucutre.
*/
void setGraph()
{
    //Asign and initialize memory for the double array of integers
    //Because is a doble array it is multiplied by 2 only
    VG = (int *) calloc(NUMVERTICES*2, NUMVERTICES*2*sizeof(int) );
    
    int numberEdges = 0;
    int randValue = 0;

    //Inicializacion valores VG
    for (int i = 0; i < NUMVERTICES; i++)
    {
        //Set the index of the Vertex
        VG[i] = i;
        //Set in random way the # of vertices to 
        //wich this vertex is connected

        //#Of vertices can not be 0,becuase all the veritces
        //have to be connected so ensure that at least all
        // the nodes are connected to at leas 2 other vertices
        randValue = rand() % (NUMVERTICES-2) +2;
        VG[i+NUMVERTICES] = randValue;
        //Keep track of the number edges
        numberEdges = numberEdges + randValue;

    }//Fin for 1

    //!!!SAVE IN GLOBAL VARIABLE!!!
    NUMBEREDGES = numberEdges;

    //----------------
    printf("-- Source Vertex --\n");
    printDoubleArray(VG);
    printf("------\n");
    printf("TOTAL EDGES: %i\n",numberEdges);
    //----------------

    //Asign and initialize memory for the double array of integers
    //Because is a doble array it is multiplied by 2 only
    EG = (int *) calloc(numberEdges*2,numberEdges*2*sizeof(int));

    //Initialize EDGE Double array values
    int indxEdges = 0;

    for(int i = 0; i < NUMVERTICES; i++)
    {
        //Num of vertices to wich vertex i has a path
        int numVerticesConn = VG[i+NUMVERTICES];

        //1)Set the destinatio id of the vertex, which can not
        //be repeated and can not be the same as the source
        //vertex i
        int indxDestination = 0;

        //2)Set randomly the value of the weight of edge 1)
        //values of weight from 1 - 100 
        for(int j = 0; j < numVerticesConn; j++)
        {
            //1)
            //Ojo: indxDestination = j a menos que se encuentr
            //que source = destino; en ese caso y por el resto
            //del for j, indxDestination ira uno arriba que j
            if( i == j )
            {
                indxDestination++;
            }//End if

            EG[indxEdges] = indxDestination;

            //2)
            EG[indxEdges+numberEdges] = rand() % (100) +1;

            indxEdges++;
            indxDestination++; 

        }//End for 3

    }//Fin for 2

    //----------------
    //printf("-- Destination vertex --\n");
    //printArrayRange(EG,0,numberEdges-1);
    //printf("-- Weigth of Edge --\n");
    //printArrayRange(EG,numberEdges,(numberEdges*2)-1);
    //----------------

}//Fin funcion setGraph
//--------------------------------
//Function that initializes values of R1,R2,R3 according to
//the Root vertex; and also define and initializes with 0s
//the temporal arrays
void setVariables()
{
    //Rs length = |NUMVERTICES|-1 because final path always
    //has one less than the #of vertices 
    R1source = (int *)calloc(NUMVERTICES-1,NUMVERTICES-1*sizeof(int));
    R2destination = (int *)calloc(NUMVERTICES-1,NUMVERTICES-1*sizeof(int));
    R3weigth = (int *)calloc(NUMVERTICES-1,NUMVERTICES-1*sizeof(int));

    //Look for the actual weights in VE and VG
    //for the source and destination and in case
    //of not being found asign 0 as the weight
    int numDestinations = VG[C+NUMVERTICES];

    int startIndex = 0;
    for(int k = 0; k < C; k++)
    {
        startIndex = startIndex+VG[k+NUMVERTICES];
    }//End for

    numDestinations = numDestinations+startIndex;

    //----------
    //printf("Range of values in EG(%i - %i)\n", startIndex, numDestinations);
    //----------
    
    //Set by default all the edges taking as the origin 
    //the root source, to all posible destinations
    int indxValidDestinations = -1;
    for(int i = 0; i < NUMVERTICES; i++)
    {

        //Only do not take as destination when source
        //and destination are equal
        if(C != i)
        {
            indxValidDestinations++;

            //Set source index
            R1source[indxValidDestinations] = C;
            R2destination[indxValidDestinations] = i;

            //Recorrer solamente los destinos para el source
            for(int j = startIndex; j < numDestinations; j++)
            {
                int idDestino = EG[j];

                //Se encontro el destino .:. poner peso correspondiente
                //----------
                //printf("%i == %i\n", idDestino, i);
                //----------
                if(idDestino == i)
                {
                    R3weigth[indxValidDestinations] = EG[j+NUMBEREDGES];
                }//End if

            }//End for 2

        }//End if


    }//Fin for 1

    //--------------
    //Recordar que para el print se considera un elemento menos
    //del limite superior ya que realmnete hace el print hasta
    //la posicion indicada
    //printf("R1: \n");
    //printArrayRange(R1source,0,NUMVERTICES-2);
    //printf("R2: \n");
    //printArrayRange(R2destination,0,NUMVERTICES-2);
    //printf("R3: \n");
    //printArrayRange(R3weigth,0,NUMVERTICES-2);
    //--------------


}//End fucntions setVariables
//--------------------------------
__global__ void kernel1(int *v, int *e, int *r1, int *r2, int *r3, int *c, int *t1, int *t2)
{
    //Define and construct T1 and T2
    //T1weights = (int *)calloc(MAXTRHEADSXBLOCK,MAXTRHEADSXBLOCK*sizeof(int));
    //T2indexes = (int *)calloc(MAXTRHEADSXBLOCK,MAXTRHEADSXBLOCK*sizeof(int));

    int idBloque = blockIdx.x;
    //ID de cada hilo (IDHILOBLOQUE+IDBLOQUE*HILOSPORBLOQUE)
    int i = threadIdx.x + idBloque*blockDim.x;

    //MIN REDUCTION AND WRITE RESULTS IN T1 AND T2

    //1)All threads in the grid make reduction operation on an array
    //of input data, ann obtain min weight and index of each thread
    
    //Solo trabajar |v|-1 hilos 
    //V-1 porque Rs son de size |V|-1
    if( i < NUMVERTICES-1 )
    {
        //----------------------
        //printf("| idh: %i | ", i);
        //printf("| %i %i | ", v[i],v[i+NUMVERTICES]);
        //----------------------

        //---------------------
        //printf("| %i  %i : %i | ", r1[i], r2[2], r3[2]);
        //printf(" %i < %i //", r3[i], t1[idBloque]);
        //---------------------
        //Con weiht mwnor al actual pero que sea un
        //weigth valido (diferente de 0)
        if(r3[i] < t1[idBloque] && r3[i] != 0 )
        {
            //Guardar Weight
            t1[idBloque] = r3[i];

            //Guardar Indice
            t2[idBloque] = r2[i];

        }//Nuevo menor encontrado

    }//End if

    //printf("| %i | ", r3[i]);

    //i < MAXNUMBEREDGES
    /*if(i < 15)
    {
                                 //[i+MAXNUMBEREDGES]
        printf("! %i %i ! ", e[i],e[i+15]);

    }//End if*/


    //2)All threads in every block make reduction of the result data in 1)
    //And obtain the minim value and index of every thread block

}//End function kernel1
//--------------------------------
__global__ void kernel2(int *numBlocks, int *weights, int *indxs)
{
    int N = numBlocks[0];

    //Reservar espacio en zona de memoria compartida
    __shared__ int temporal[MAXTRHEADSXBLOCK];
    __shared__ int tempids[MAXTRHEADSXBLOCK];

    //Indice de cada hilo en un solo bloque
    int i = threadIdx.x;
    
    if(i < N)
    {
        //Copiamos el vector de pesos en temporal y sincronizamos
        temporal[i] = weights[i];
        tempids[i] = indxs[i];
        __syncthreads();

        //---------------------
        //printf("|%i)  %i : %i | ", i ,weights[i], indxs[i]);
        //printf("| %i | ", temporal[i]);
        //----------------------

        //Inicio de reduccion paralela
        int salto = N/2;

        //log2(N) iteraciones
        while(salto)
        {
            //Solo trabajan la mitad de los hilos
            if(i < salto)
            {
                //Si se encuentra un vertex con peso menor se elige
                //como mejor candidato
                if( temporal[i+salto] < temporal[i] && temporal[i] != 0 )
                {
                    temporal[i] = temporal[i+salto];
                    tempids[i] = tempids[i+salto];

                }//End if

            }//End if 2
            __syncthreads();
            salto = salto/2;

        }//End while

        //Hilo 0 escibe el resultado final en memoria global
        if(i == 0)
        {
            weights[0] = temporal[0];
            indxs[0] = tempids[0];

        }//End if 2

    }//End if 1 


}//End function kernel2
//--------------------------------
//Comparing and update MST
__global__ void kernel3(int *v, int *e, int *r1, int *r3, int *c, int *numEdges)
{
    //1)Read current Vertex index

    //2)Fin the weight between current vertex and the
    //other vertices (n other Vertex in actual moment)

    //For every Vertex n if new weight with this C vertex is < old weight
    // Adjust corresponding values of R1 and R3 by :
    //if(W[n] < R3[n] )
    //R1[n] = C
    //R3[n] = W[n]

    //Look for the actual weights in VE and VG
    //for the source and destination and in case
    //of not being found asign 0 as the weight
    int NE = numEdges[0];
    int C = c[0];
    int numDestinations = v[C+NUMVERTICES];

    int startIndex = 0;
    for(int k = 0; k < C; k++)
    {
        startIndex = startIndex+v[k+NUMVERTICES];
    }//End for

    numDestinations = numDestinations+startIndex;

    //----------
    //printf("Range of values in EG(%i - %i)\n", startIndex, numDestinations);
    //----------

    int idBloque = blockIdx.x;
    //ID de cada hilo (IDHILOBLOQUE+IDBLOQUE*HILOSPORBLOQUE)
    int i = threadIdx.x + idBloque*blockDim.x;
    
    //Set by default all the edges taking as the origin 
    //the root source, to all posible destinations
    int indxValidDestinations = -1;
    if(i < NUMVERTICES)
    {

        //Only do not take as destination when source
        //and destination are equal
        if(C != i)
        {
            indxValidDestinations++;

            //Recorrer solamente los destinos para el source
            for(int j = startIndex; j < numDestinations; j++)
            {
                int idDestino = e[j];

                //Se encontro el destino .:. poner peso correspondiente
                //----------
                //printf("%i == %i\n", idDestino, i);
                //----------
                if(idDestino == i)
                {

                    if(e[j+NE] < r3[indxValidDestinations])
                    {
                        r3[indxValidDestinations] = e[j+NE];
                        r1[indxValidDestinations] = C;

                    }//End if

                }//End if

            }//End for 2

        }//End if


    }//Fin for 1
    

}//End function kernel3
//--------------------------------
void primMST(int *v, int *e, int *r1, int * r2, int *r3, int c)
{
    //Define size of CUDA grid
    int g_row = (int)ceil((float)NUMVERTICES/(float)MAXTRHEADSXBLOCK);
    int g_col = (int)ceil((float)NUMVERTICES/(float)MAXTRHEADSXBLOCK); 
    int numBloques = g_row;
    dim3 bloques(g_col,g_row);
    dim3 hilos(MAXTRHEADSXBLOCK,MAXTRHEADSXBLOCK);

    cudaEvent_t start, stop; 

    printf("Bloques: %i == %i \n", bloques, numBloques);
    printf("Hilos: %i \n", hilos);
    printf("Grid (%d,%d)\n", g_row, g_col); 

    //vARIABLES IN DEVICE
    int *VGD, *VED, *R1D, *R2D, *R3D;   //Arrays
    int *T1D, *T2D;
    int *CD, *NED;                      //Variables 

    //Define and construct T1 and T2? HERE
    T1weights = (int *)calloc(numBloques,numBloques*sizeof(int));
    T2indexes = (int *)calloc(numBloques,numBloques*sizeof(int));

    //Initialize temporal weights with a very high value
    //in order to make that any wieght is better than 
    //the init value
    for(int i = 0; i < numBloques; i++)
    {
        T1weights[i] = 99999;

    }//End for 1

    //--------------
    //Recordar que para el print se considera un elemento menos
    //del limite superior ya que realmnete hace el print hasta
    //la posicion indicada
    printf("R1: \n");
    printArrayRange(r1,0,NUMVERTICES-2);
    printf("R2: \n");
    printArrayRange(r2,0,NUMVERTICES-2);
    printf("R3: \n");
    printArrayRange(r3,0,NUMVERTICES-2);
    //--------------

    //TRANSFER FROM HOST (CPU) TO DEVICE GPU
    cudaSetDevice(gpudev);

    cudaEventCreate(&start); cudaEventCreate(&stop);
    cudaEventRecord(start,0); 

    //1)Asignar memoria para variables en GPU
    cudaMalloc(&VGD, NUMVERTICES*2*sizeof(int) );
    cudaMalloc(&VED, NUMBEREDGES*2*sizeof(int) );
    cudaMalloc(&R1D, (NUMVERTICES-1)*sizeof(int) );
    cudaMalloc(&R2D, (NUMVERTICES-1)*sizeof(int) );
    cudaMalloc(&R3D, (NUMVERTICES-1)*sizeof(int) );
    cudaMalloc(&T1D, (numBloques)*sizeof(int) );
    cudaMalloc(&T2D, (numBloques)*sizeof(int) );
    cudaMalloc(&CD, int(sizeof(int)) );
    cudaMalloc(&NED, int(sizeof(int)) );


    //INICIO LOOP |NUMVERTICES|-1 VECES
    int iCountVer = 0;

    while(iCountVer < NUMVERTICES-1)
    {
        //----
        printf("---- %i) ----- \n", iCountVer);
        printf("Current vertex: %i \n", c);
        //----

        //2)Copiar datos del host al device
        cudaMemcpy(VGD,v,NUMVERTICES*2*sizeof(int),cudaMemcpyDefault);
        cudaMemcpy(VED,e,NUMBEREDGES*2*sizeof(int),cudaMemcpyDefault);
        cudaMemcpy(R1D,r1,(NUMVERTICES-1)*sizeof(int),cudaMemcpyDefault);
        cudaMemcpy(R2D,r2,(NUMVERTICES-1)*sizeof(int),cudaMemcpyDefault);
        cudaMemcpy(R3D,r3,(NUMVERTICES-1)*sizeof(int),cudaMemcpyDefault);
        cudaMemcpy(T1D,T1weights,numBloques*sizeof(int),cudaMemcpyDefault);
        cudaMemcpy(T2D,T2indexes,numBloques*sizeof(int),cudaMemcpyDefault);
        cudaMemcpy(CD,&c,sizeof(int),cudaMemcpyDefault);
        cudaMemcpy(NED,&NUMBEREDGES,sizeof(int),cudaMemcpyDefault);

        //3)Ejecutar kernel
        //INVOQUE KERNEL 1 AND WRITE RESULTS IN T1 AND T2
        kernel1<<<bloques, hilos>>>(VGD,VED,R1D,R2D,R3D,CD,T1D,T2D);

        //4)Copiar datos del device al host
        //T1 Y T2
        
        // Valores de T1[0] y T2[0] son añadidos
        // a los correspondientes R1 Y R3
        //T2[0] sobreescribe a C
        cudaMemcpy(T1weights,T1D,numBloques*sizeof(int),cudaMemcpyDefault);
        cudaMemcpy(T2indexes,T2D,numBloques*sizeof(int),cudaMemcpyDefault);
        //---------------
        printf("\n Minimum weight found for each block (Global memory reduction) \n");
        printf("Id: \n");
        printArrayRange(T2indexes,0,numBloques-1);
        printf("Weight: \n");
        printArrayRange(T1weights,0,numBloques-1);
        //---------------

        //Verificar si se inicia al Kernel 2
        //MAXTRHEADSXBLOCK > numBloques > 1
        if(numBloques > 1)
        {
            //Definir variable en device
            int *NBD;

            //1)Asinar memoria para vairable en GPU/device
            cudaMalloc(&NBD, int(sizeof(int)) );

            //2)Copiar datos del host al device
            cudaMemcpy(NBD,&numBloques,sizeof(int),cudaMemcpyDefault);

            //3)ejecutar kermel2
            printf("Invoke Kernel2\n");
            kernel2<<<1,hilos>>>(NBD,T1D,T2D);

            //4)Copiar datos del device al host
            cudaMemcpy(T1weights,T1D,numBloques*sizeof(int),cudaMemcpyDefault);
            cudaMemcpy(T2indexes,T2D,numBloques*sizeof(int),cudaMemcpyDefault);

            //---------------
            printf("\n 2) Minimum weight found of each block (After shared memory reduction) \n");
            printf("Id: \n");
            printArrayRange(T2indexes,0,numBloques-1);
            printf("Weight: \n");
            printArrayRange(T1weights,0,numBloques-1);
            //---------------
            
            //5)liberar memoria NBD

        }//End if kernel2

        //---------------
        printf("Minimum weight found: %i for vertex with ID: %i \n", T1weights[0], T2indexes[0]);
        //---------------

        //ADDING NEW MST EDGE

        //1)Add previously found minimum weight Edge (T1[0] WEIGHT T2[0] EDGE)
        //to the MST by moving this Edge to the first position of R1 R2 R3
        r2[0] = T2indexes[0];  //R2[C] ? 
        r3[0] = T1weights[0];  //R3[C] ?

        //2)Saving current Vertex C= R2[T2[0]]
        c = r2[T2indexes[0]];

        //Copiar datos del host al device que han sido modificados
        cudaMemcpy(R2D,r2,(NUMVERTICES-1)*sizeof(int),cudaMemcpyDefault);
        cudaMemcpy(R3D,r3,(NUMVERTICES-1)*sizeof(int),cudaMemcpyDefault);
        cudaMemcpy(CD,&c,sizeof(int),cudaMemcpyDefault);

        //Kernel 3: Comparing and updating MST
        kernel3<<<bloques, hilos>>>(VGD,VED,R1D,R3D,CD,NED);

        //Copiar datos del device al host
        cudaMemcpy(r1,R1D,(NUMVERTICES-1)*sizeof(int),cudaMemcpyDefault);
        cudaMemcpy(r3,R3D,(NUMVERTICES-1)*sizeof(int),cudaMemcpyDefault);

        
        //--------------
        //Recordar que para el print se considera un elemento menos
        //del limite superior ya que realmnete hace el print hasta
        //la posicion indicada
        printf("--- MST ACTUALIZADO ---: \n");
        printf("R1: \n");
        printArrayRange(r1,0,NUMVERTICES-2);
        printf("R2: \n");
        printArrayRange(r2,0,NUMVERTICES-2);
        printf("R3: \n");
        printArrayRange(r3,0,NUMVERTICES-2);
        //--------------

        iCountVer++;

    }//End while
    
    //FIN LOOP |NUMVERTICES|-1 VECES


    //5) Liberar Memoria
    cudaFree(VGD);
    cudaFree(VED);
    cudaFree(R1D);
    cudaFree(R2D);
    cudaFree(R3D);
    cudaFree(T1D);
    cudaFree(T2D);
    cudaFree(CD);
    cudaFree(NED);

    cudaEventRecord(stop, 0); 
    cudaEventSynchronize(stop);

}//En function primMST
//---- FIN FUNCIONES -----
//Inicio del programa
int main(int argc, char **argv)
{
    setGraph();
    
    //Set root vertex of the MST
    C = 2;

    setVariables();

    printf("IDs threads: \n");
    primMST(VG,EG,R1source,R2destination,R3weigth,C);
    printf("\n");


    printf("Fin del programa V4\n");

}//Fin del main