#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <stdlib.h>

//VARIABLES GLOBALES
#define NUMVERTICES 10

//Graph representation with
int *EG;    //Double array of edges |NUMVERTICES|
int *VG;    //Double array of vertices

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
    printf("-- Destination vertex --\n");
    printArrayRange(EG,0,numberEdges-1);
    printf("-- Weigth of Edge --\n");
    printArrayRange(EG,numberEdges,(numberEdges*2)-1);
    //----------------

}//Fin funcion setGraph
//---- FIN FUNCIONES -----
//Inicio del programa
int main(int argc, char **argv)
{
    setGraph();
    printf("Fin del programa V1\n");

}//Fin del main